#!/usr/bin/env python3
"""Re-filter DupCaller SBS/indel calls by their own read family.

For every record in a DupCaller VCF, the read family that produced the call
is rebuilt from the BAM the same way DupCaller groups reads
(DupCaller_sub/funcs/misc.py:_compute_read_label + callBam's read filters):

    reads with reference_start == INFO/SP, orientation-normalised duplex
    barcode == INFO/TAG1+TAG2, and identical signed template_length
    (plus, for downstream-mate families, the "rugged" mates DupCaller
    reroutes into the family at SP -- mates whose own alignment start is a few
    bp earlier, e.g. 3' soft clip or an indel near the read start).

INFO/TL is used if present; otherwise the template length is picked among the
TL groups at SP whose F1R2/F2R1 read counts match INFO/F1R2,F2R1.

A record passes when
  1) every read of that family supports the variant (SBS/MNV: alt base at
     every position; indel: carries an indel that yields the same haplotype
     as the VCF allele), and
  2) the per-strand consensus base quality is > --min-cq on BOTH strands
     (F1R2 and F2R1). SBS: at the variant base (min over bases for MNVs).
     Indel: at anchor+1 -- the read base right after the VCF anchor base in
     read space (first inserted base / first base after the deletion).

Consensus quality of a strand = posterior Phred quality of the strand
consensus base under the standard independent-error model
(P(read base = b | true base = b) = 1-e, else e/3, e = 10^(-BQ/10), flat
prior). When all n reads agree it is ~ sum(BQ) + 10*log10(3)*(n-1), i.e. a
single Q37 read gives ~37, two agreeing Q37 reads ~79.

Usage:
  family_consensus_filter.py -i calls.vcf -b sample.bam -o out.vcf [-f ref.fa]
      [--min-cq 60] [--all] [--allow-noncovering]
"""
import argparse
import math
import sys
from collections import defaultdict

import pysam

DEFAULT_REF = "GRCh38.d1.vd1.fa"
INDEL_WINDOW = 300  # bp of reference on each side used for haplotype comparison
RUGGED_WINDOW = 300  # max bp a rerouted downstream mate may start before SP


# --------------------------------------------------------------------------
# family reconstruction (mirrors DupCaller callBam)
# --------------------------------------------------------------------------
def dupcaller_read_ok(rec):
    if rec.is_supplementary or rec.is_secondary or rec.is_qcfail:
        return False
    if not rec.is_proper_pair or rec.cigartuples is None:
        return False
    for t in ("DT", "dt"):
        if rec.has_tag(t) and rec.get_tag(t) == "SQ":
            return False
    # 5' soft clip -> skipped by DupCaller
    if rec.is_forward and rec.cigartuples[0][0] != 0:
        return False
    if rec.is_reverse and rec.cigartuples[-1][0] != 0:
        return False
    return True


def is_f1r2(rec):
    return (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2)


def family_label(rec, bc_tag, bc_sep):
    bc1, bc2 = rec.get_tag(bc_tag).split(bc_sep)
    if is_f1r2(rec):
        return bc1, bc2
    return bc2, bc1


def _label_ok(rec, tags, bc_tag, bc_sep):
    return (
        dupcaller_read_ok(rec)
        and rec.has_tag(bc_tag)
        and family_label(rec, bc_tag, bc_sep) == tags
    )


def _reads_starting_at(bam, chrom, pos, tags, bc_tag, bc_sep):
    return [
        r
        for r in bam.fetch(chrom, pos, pos + 1)
        if r.reference_start == pos and _label_ok(r, tags, bc_tag, bc_sep)
    ]


def dupcaller_families_at(bam, chrom, sp, tags, bc_tag, bc_sep):
    """{signed TL: reads} of the families DupCaller assembles at position sp
    for this barcode pair. Downstream mates (TL < 0) follow DupCaller's
    rugged-mate rerouting (misc.py:_index_rugged_mates): when the upstream
    family's members disagree on their mate start, every mate not at the
    largest mate start is processed in the family at that largest start."""
    upstream = {}

    def processed_at(rec):
        if rec.template_length >= 0:
            return rec.reference_start
        key = (rec.next_reference_start, -rec.template_length)
        if key not in upstream:
            ups = [
                u
                for u in _reads_starting_at(bam, chrom, key[0], tags, bc_tag, bc_sep)
                if u.template_length == key[1]
            ]
            maxm = max((u.next_reference_start for u in ups), default=None)
            rugged = {u.query_name for u in ups if u.next_reference_start != maxm}
            upstream[key] = (maxm, rugged)
        maxm, rugged = upstream[key]
        return maxm if rec.query_name in rugged else rec.reference_start

    # materialise first: processed_at() fetches again on the same handle,
    # which would disturb a live fetch iterator
    cands = [
        rec
        for rec in bam.fetch(chrom, max(0, sp - RUGGED_WINDOW), sp + 1)
        if sp - RUGGED_WINDOW <= rec.reference_start <= sp
        and _label_ok(rec, tags, bc_tag, bc_sep)
    ]
    groups = defaultdict(list)
    for rec in cands:
        if processed_at(rec) == sp:
            groups[rec.template_length].append(rec)
    return groups


def get_family(bam, chrom, sp, tag1, tag2, info, bc_tag, bc_sep):
    """Return (reads, tl, status)."""
    groups = dupcaller_families_at(bam, chrom, sp, (tag1, tag2), bc_tag, bc_sep)
    if not groups:
        return [], None, "family_not_found"

    if "TL" in info:
        tl = int(info["TL"])
        if tl not in groups:
            return [], tl, "family_not_found"
        return groups[tl], tl, "ok"

    if len(groups) == 1:
        tl = next(iter(groups))
        return groups[tl], tl, "ok"

    want = (int(info.get("F1R2", -1)), int(info.get("F2R1", -1)))
    matches = [
        tl
        for tl, rs in groups.items()
        if (sum(is_f1r2(r) for r in rs), sum(not is_f1r2(r) for r in rs)) == want
    ]
    if len(matches) == 1:
        return groups[matches[0]], matches[0], "ok"
    return [], None, "ambiguous_family"


# --------------------------------------------------------------------------
# consensus quality
# --------------------------------------------------------------------------
def consensus_quality(obs):
    """obs: list of (base, bq). Returns (consensus_base, phred_q)."""
    if not obs:
        return None, 0.0
    ll = {}
    for b in "ACGT":
        s = 0.0
        for base, q in obs:
            e = min(10 ** (-q / 10.0), 0.75)
            s += math.log10(1 - e) if base == b else math.log10(e / 3)
        ll[b] = s
    best = max(ll, key=ll.get)
    # log10 P(not best) - log10 P(best), computed stably
    others = [ll[b] - ll[best] for b in ll if b != best]
    m = max(others)
    log_err = m + math.log10(sum(10 ** (o - m) for o in others))
    log_err -= math.log10(1 + 10 ** log_err)  # normalise by total
    return best, -10 * log_err


def strand_cq(per_read, want_base):
    """per_read: list of (is_f1r2, base, bq). Returns {strand: q}, where q is
    the consensus quality if the strand consensus equals want_base else 0
    (want_base None -> accept whatever the consensus is)."""
    out = {}
    for strand, name in ((True, "F1R2"), (False, "F2R1")):
        obs = [(b, q) for s, b, q in per_read if s == strand and b is not None]
        cb, q = consensus_quality(obs)
        if cb is None or (want_base is not None and cb != want_base):
            q = 0.0
        out[name] = q
    return out


# --------------------------------------------------------------------------
# SBS / MNV
# --------------------------------------------------------------------------
def eval_sbs(reads, pos0, ref, alt):
    """Returns (n_support, n_reads, n_noncover, min_cq_by_strand)."""
    n_support = n_noncover = 0
    cq = {"F1R2": float("inf"), "F2R1": float("inf")}
    per_pos_reads = [[] for _ in alt]
    for rec in reads:
        r2q = {
            rp: qp
            for qp, rp in rec.get_aligned_pairs(matches_only=True)
            if pos0 <= rp < pos0 + len(alt)
        }
        seq, qual = rec.query_sequence, rec.query_qualities
        supports, covers = True, True
        for k, a in enumerate(alt):
            qp = r2q.get(pos0 + k)
            if qp is None:
                covers = supports = False
                per_pos_reads[k].append((is_f1r2(rec), None, 0))
                continue
            b = seq[qp]
            if b != a:
                supports = False
            per_pos_reads[k].append((is_f1r2(rec), b, qual[qp]))
        n_support += supports
        n_noncover += not covers
    for k, a in enumerate(alt):
        q = strand_cq(per_pos_reads[k], a)
        for s in cq:
            cq[s] = min(cq[s], q[s])
    return n_support, n_noncover, cq


# --------------------------------------------------------------------------
# indel
# --------------------------------------------------------------------------
def read_indel_events(rec):
    """Yield (op, ref_pos_after_anchor, length, inserted_seq, anchor_qpos)
    for every I/D op preceded by an aligned base."""
    rpos, qpos = rec.reference_start, 0
    prev_m = False
    for op, ln in rec.cigartuples:
        if op in (0, 7, 8):
            rpos += ln
            qpos += ln
            prev_m = True
            continue
        if op == 1:
            if prev_m:
                yield ("I", rpos, ln, rec.query_sequence[qpos : qpos + ln], qpos - 1)
            qpos += ln
        elif op in (2, 3):
            if prev_m and op == 2:
                yield ("D", rpos, ln, "", qpos - 1)
            rpos += ln
        elif op == 4:
            qpos += ln
        prev_m = False


def eval_indel(reads, fasta, chrom, pos0, ref, alt):
    """pos0 = 0-based anchor. Returns (n_support, n_noncover, cq)."""
    is_ins = len(alt) > len(ref)
    ln = abs(len(alt) - len(ref))
    w0 = max(0, pos0 - INDEL_WINDOW)
    w1 = pos0 + len(ref) + INDEL_WINDOW
    win = fasta.fetch(chrom, w0, w1).upper()
    vcf_hap = win[: pos0 - w0] + alt.upper() + win[pos0 - w0 + len(ref) :]

    n_support = n_noncover = 0
    per_read = []
    for rec in reads:
        seq, qual = rec.query_sequence, rec.query_qualities
        match_qpos = None
        for op, rp, l, ins, aq in read_indel_events(rec):
            if l != ln or (op == "I") != is_ins or not (w0 < rp < w1 - ln):
                continue
            o = rp - w0
            hap = win[:o] + ins + win[o:] if op == "I" else win[:o] + win[o + l :]
            if hap == vcf_hap:
                match_qpos = aq
                break
        if match_qpos is None:
            # does the read even span the anchor and the base after it?
            spans = rec.reference_start <= pos0 and rec.reference_end > pos0 + len(ref)
            n_noncover += not spans
            per_read.append((is_f1r2(rec), None, 0))
            continue
        n_support += 1
        # anchor+1 in read space: the read base after the VCF (left-aligned)
        # anchor; fall back to the read's own event anchor if the VCF anchor
        # is not aligned in this read.
        qa = None
        for qp, rp in rec.get_aligned_pairs(matches_only=True):
            if rp == pos0:
                qa = qp
                break
        if qa is None or qa > match_qpos:
            qa = match_qpos
        qn = qa + 1
        if qn >= len(seq):
            per_read.append((is_f1r2(rec), None, 0))
        else:
            per_read.append((is_f1r2(rec), seq[qn], qual[qn]))
    cq = strand_cq(per_read, None)
    return n_support, n_noncover, cq


# --------------------------------------------------------------------------
# VCF I/O (plain text: tolerant of trailing NUL padding)
# --------------------------------------------------------------------------
NEW_HEADER = [
    '##INFO=<ID=FAM_TL,Number=1,Type=Integer,Description="Signed template length of the read family used by family_consensus_filter">',
    '##INFO=<ID=FAM_N,Number=3,Type=Integer,Description="Family reads: total,supporting the variant,not covering the site">',
    '##INFO=<ID=FAM_CQ,Number=2,Type=Float,Description="Consensus base quality F1R2,F2R1 (SBS: variant base, min over MNV bases; indel: base at anchor+1)">',
    '##FILTER=<ID=family_not_found,Description="No reads with SP/TAG1/TAG2 (and TL) found in the BAM">',
    '##FILTER=<ID=ambiguous_family,Description="Several template-length groups at SP match TAG1/TAG2 and F1R2/F2R1; family not determinable">',
    '##FILTER=<ID=family_discordant,Description="Not every read of the read family (same TAG1,TAG2,TL,SP) supports the variant">',
    '##FILTER=<ID=low_consensus_bq,Description="Consensus base quality <= threshold on at least one strand">',
]


def parse_info(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-i", "--vcf", required=True)
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("-f", "--reference", default=DEFAULT_REF, help="needed for indels")
    ap.add_argument("--min-cq", type=float, default=60, help="strict >, default 60")
    ap.add_argument("--all", action="store_true",
                    help="write every record with FILTER annotated (default: passing only)")
    ap.add_argument("--allow-noncovering", action="store_true",
                    help="ignore family reads that do not cover the site in criterion 1 "
                         "(default: such reads count as not supporting)")
    ap.add_argument("--barcode-tag", default="DB")
    ap.add_argument("--barcode-sep", default="-")
    a = ap.parse_args()

    bam = pysam.AlignmentFile(a.bam)
    fasta = None
    n_in = n_pass = 0
    reasons = defaultdict(int)

    with open(a.vcf) as fin, open(a.out, "w") as fout:
        for line in fin:
            line = line.rstrip("\n").rstrip("\x00")
            if not line:
                continue
            if line.startswith("##"):
                fout.write(line + "\n")
                continue
            if line.startswith("#"):
                fout.write("\n".join(NEW_HEADER) + "\n" + line + "\n")
                continue
            f = line.split("\t")
            chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
            info = parse_info(f[7])
            n_in += 1

            reads, tl, status = get_family(
                bam, chrom, int(info["SP"]), info["TAG1"], info["TAG2"],
                info, a.barcode_tag, a.barcode_sep,
            )
            fails = []
            if status != "ok":
                fails.append(status)
                n_sup = n_nc = 0
                cq = {"F1R2": 0.0, "F2R1": 0.0}
            elif len(ref) == len(alt):
                n_sup, n_nc, cq = eval_sbs(reads, pos - 1, ref, alt)
            else:
                if fasta is None:
                    fasta = pysam.FastaFile(a.reference)
                n_sup, n_nc, cq = eval_indel(reads, fasta, chrom, pos - 1, ref, alt)

            if status == "ok":
                needed = len(reads) - (n_nc if a.allow_noncovering else 0)
                if n_sup < needed or needed == 0:
                    fails.append("family_discordant")
                if not (cq["F1R2"] > a.min_cq and cq["F2R1"] > a.min_cq):
                    fails.append("low_consensus_bq")

            for r in fails:
                reasons[r] += 1
            passed = not fails
            n_pass += passed
            if not passed and not a.all:
                continue

            extra = [
                f"FAM_N={len(reads)},{n_sup},{n_nc}",
                f"FAM_CQ={cq['F1R2']:.2f},{cq['F2R1']:.2f}",
            ]
            if tl is not None:
                extra.insert(0, f"FAM_TL={tl}")
            f[7] = f[7] + ";" + ";".join(extra)
            if not passed:
                f[6] = ";".join(fails) if f[6] in ("PASS", ".") else f[6] + ";" + ";".join(fails)
            fout.write("\t".join(f) + "\n")

    msg = f"{a.vcf}: {n_pass}/{n_in} pass"
    if reasons:
        msg += " | " + ", ".join(f"{k}={v}" for k, v in sorted(reasons.items()))
    print(msg, file=sys.stderr)


if __name__ == "__main__":
    main()
