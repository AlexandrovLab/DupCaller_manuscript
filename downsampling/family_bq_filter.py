#!/usr/bin/env python3
"""Re-filter DupCaller SBS/indel calls by their read family, BQ-aware.

The read family of each call (same TAG1/TAG2, template length and start
position SP, incl. DupCaller's rugged-mate rerouting) is rebuilt exactly as
in family_consensus_filter.py (imported from the same directory).

Each family read is classified at the call site as
  support    - carries the variant (SBS/MNV: alt base at every position;
               indel: an indel giving the same haplotype as the VCF allele)
  nonsupport - covers the site but does not carry the variant
  noncover   - does not cover the site (counted as nonsupport with unknown BQ
               unless --allow-noncovering, in which case it is ignored)
and given a BQ at the site:
  SBS:   base quality at the variant base (MNV: min over the bases for a
         supporting read; max over the non-alt bases for a nonsupporting read,
         i.e. it is "low BQ" only if every disagreeing base is low BQ)
  indel: base quality at anchor+1 in read space (supporting read: first
         inserted base / first base after the deletion; nonsupporting read:
         the read base after the base aligned to the VCF anchor)

A record passes when
  1) >= --min-support (3) supporting reads have BQ >= --min-bq (18), and
     every family read with BQ >= 18 supports the variant;
  2) at most --max-nonsupport (1) read does not support the variant, and each
     such read has BQ < 18;
  3) the consensus base quality is > --min-cq (60) on BOTH strands (F1R2 and
     F2R1): SBS at the variant base (min over MNV bases), indel at anchor+1.
     Consensus = posterior Phred quality under the independent-error model
     (see family_consensus_filter.consensus_quality), computed from the
     supporting reads only; with --cq-include-dissent the nonsupporting read
     enters the strand consensus as a disagreeing observation.

Usage:
  family_bq_filter.py -i calls.vcf -b sample.bam -o out.vcf [-f ref.fa]
      [--min-bq 18] [--min-support 3] [--max-nonsupport 1] [--min-cq 60]
      [--cq-include-dissent] [--all] [--allow-noncovering]
"""
import argparse
import os
import sys
from collections import defaultdict

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_consensus_filter import (  # noqa: E402
    DEFAULT_REF,
    INDEL_WINDOW,
    consensus_quality,
    get_family,
    is_f1r2,
    parse_info,
    read_indel_events,
)


# --------------------------------------------------------------------------
# per-read classification
# --------------------------------------------------------------------------
# each read -> (strand_is_f1r2, status, bq, [(base, bq) per site position])
# status in {"support", "nonsupport", "noncover"}


def classify_sbs(reads, pos0, alt):
    out = []
    for rec in reads:
        r2q = {
            rp: qp
            for qp, rp in rec.get_aligned_pairs(matches_only=True)
            if pos0 <= rp < pos0 + len(alt)
        }
        seq, qual = rec.query_sequence, rec.query_qualities
        obs = []
        for k in range(len(alt)):
            qp = r2q.get(pos0 + k)
            obs.append(None if qp is None else (seq[qp], qual[qp]))
        if any(o is None for o in obs):
            out.append((is_f1r2(rec), "noncover", None, obs))
        elif all(b == a for (b, _), a in zip(obs, alt)):
            out.append((is_f1r2(rec), "support", min(q for _, q in obs), obs))
        else:
            bq = max(q for (b, q), a in zip(obs, alt) if b != a)
            out.append((is_f1r2(rec), "nonsupport", bq, obs))
    return out


def classify_indel(reads, fasta, chrom, pos0, ref, alt):
    is_ins = len(alt) > len(ref)
    ln = abs(len(alt) - len(ref))
    w0 = max(0, pos0 - INDEL_WINDOW)
    w1 = pos0 + len(ref) + INDEL_WINDOW
    win = fasta.fetch(chrom, w0, w1).upper()
    vcf_hap = win[: pos0 - w0] + alt.upper() + win[pos0 - w0 + len(ref):]

    out = []
    for rec in reads:
        seq, qual = rec.query_sequence, rec.query_qualities
        match_qpos = None
        for op, rp, l, ins, aq in read_indel_events(rec):
            if l != ln or (op == "I") != is_ins or not (w0 < rp < w1 - ln):
                continue
            o = rp - w0
            hap = win[:o] + ins + win[o:] if op == "I" else win[:o] + win[o + l:]
            if hap == vcf_hap:
                match_qpos = aq
                break
        qa = None
        for qp, rp in rec.get_aligned_pairs(matches_only=True):
            if rp == pos0:
                qa = qp
                break
        if match_qpos is not None:
            # anchor+1 after the VCF (left-aligned) anchor; fall back to the
            # read's own event anchor if the VCF anchor is not aligned here
            if qa is None or qa > match_qpos:
                qa = match_qpos
            status = "support"
        else:
            spans = rec.reference_start <= pos0 and rec.reference_end > pos0 + len(ref)
            status = "nonsupport" if spans and qa is not None else "noncover"
        qn = None if qa is None else qa + 1
        if qn is None or qn >= len(seq):
            # no anchor+1 base in this read -> no BQ to judge it by
            out.append((is_f1r2(rec), "noncover", None, [None]))
        else:
            out.append((is_f1r2(rec), status, qual[qn], [(seq[qn], qual[qn])]))
    return out


# --------------------------------------------------------------------------
# consensus quality
# --------------------------------------------------------------------------
def strand_cq(classified, alt, include_dissent):
    """{strand: q} where q is the min over site positions of the strand
    consensus quality, 0 if the consensus is not the supporting base.
    alt None (indel) -> the supporting reads' majority base at anchor+1."""
    npos = len(alt) if alt is not None else 1
    cq = {}
    for strand, name in ((True, "F1R2"), (False, "F2R1")):
        q_min = float("inf")
        for k in range(npos):
            sup = [o[k] for s, st, _, o in classified
                   if s == strand and st == "support" and o[k] is not None]
            if alt is not None:
                want = alt[k]
            else:
                want, _ = consensus_quality(sup)
            if want is None:
                q_min = 0.0
                break
            obs = list(sup)
            if include_dissent:
                for s, st, _, o in classified:
                    if s != strand or st != "nonsupport" or o[k] is None:
                        continue
                    b, q = o[k]
                    if alt is None:
                        # indel dissent: its base is not on the alt haplotype,
                        # count it as an observation of a different base
                        b = next(x for x in "ACGT" if x != want)
                    obs.append((b, q))
            cb, q = consensus_quality(obs)
            q_min = min(q_min, q if cb == want else 0.0)
        cq[name] = q_min
    return cq


# --------------------------------------------------------------------------
# VCF I/O (plain text: tolerant of trailing NUL padding)
# --------------------------------------------------------------------------
def new_header(a):
    return [
        '##INFO=<ID=FAM_TL,Number=1,Type=Integer,Description="Signed template length of the read family used by family_bq_filter">',
        f'##INFO=<ID=FAM_N,Number=5,Type=Integer,Description="Family reads: total,supporting with BQ>={a.min_bq},supporting with BQ<{a.min_bq},nonsupporting,not covering the site">',
        '##INFO=<ID=FAM_NSBQ,Number=.,Type=Integer,Description="BQ of each nonsupporting read (. = not covering the site)">',
        '##INFO=<ID=FAM_CQ,Number=2,Type=Float,Description="Consensus base quality F1R2,F2R1 (SBS: variant base, min over MNV bases; indel: base at anchor+1)">',
        '##FILTER=<ID=family_not_found,Description="No reads with SP/TAG1/TAG2 (and TL) found in the BAM">',
        '##FILTER=<ID=ambiguous_family,Description="Several template-length groups at SP match TAG1/TAG2 and F1R2/F2R1; family not determinable">',
        f'##FILTER=<ID=low_hq_support,Description="Fewer than {a.min_support} supporting reads with BQ>={a.min_bq}">',
        f'##FILTER=<ID=hq_dissent,Description="A family read not supporting the variant has BQ>={a.min_bq} (or does not cover the site)">',
        f'##FILTER=<ID=too_many_dissent,Description="More than {a.max_nonsupport} family read(s) do not support the variant">',
        f'##FILTER=<ID=low_consensus_bq,Description="Consensus base quality <= {a.min_cq:g} on at least one strand">',
    ]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("-i", "--vcf", required=True)
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("-f", "--reference", default=DEFAULT_REF, help="needed for indels")
    ap.add_argument("--min-bq", type=int, default=18, help="high-quality read base BQ, >=, default 18")
    ap.add_argument("--min-support", type=int, default=3,
                    help="min supporting reads with BQ>=--min-bq, default 3")
    ap.add_argument("--max-nonsupport", type=int, default=1,
                    help="max nonsupporting reads (each must have BQ<--min-bq), default 1")
    ap.add_argument("--min-cq", type=float, default=60, help="per-strand consensus BQ, strict >, default 60")
    ap.add_argument("--cq-include-dissent", action="store_true",
                    help="include the nonsupporting read in the strand consensus (default: supporting reads only)")
    ap.add_argument("--all", action="store_true",
                    help="write every record with FILTER annotated (default: passing only)")
    ap.add_argument("--allow-noncovering", action="store_true",
                    help="ignore family reads that do not cover the site "
                         "(default: such reads count as nonsupporting with unknown BQ -> fail)")
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
                fout.write("\n".join(new_header(a)) + "\n" + line + "\n")
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
            classified = []
            cq = {"F1R2": 0.0, "F2R1": 0.0}
            if status != "ok":
                fails.append(status)
            else:
                if len(ref) == len(alt):
                    classified = classify_sbs(reads, pos - 1, alt)
                    cq = strand_cq(classified, alt, a.cq_include_dissent)
                else:
                    if fasta is None:
                        fasta = pysam.FastaFile(a.reference)
                    classified = classify_indel(reads, fasta, chrom, pos - 1, ref, alt)
                    cq = strand_cq(classified, None, a.cq_include_dissent)

            sup_hq = sum(st == "support" and bq >= a.min_bq for _, st, bq, _ in classified)
            sup_lq = sum(st == "support" and bq < a.min_bq for _, st, bq, _ in classified)
            noncover = [c for c in classified if c[1] == "noncover"]
            dissent = [c for c in classified if c[1] == "nonsupport"]
            if not a.allow_noncovering:
                dissent += noncover
            if status == "ok":
                if sup_hq < a.min_support:
                    fails.append("low_hq_support")
                if any(bq is None or bq >= a.min_bq for _, _, bq, _ in dissent):
                    fails.append("hq_dissent")
                if len(dissent) > a.max_nonsupport:
                    fails.append("too_many_dissent")
                if not (cq["F1R2"] > a.min_cq and cq["F2R1"] > a.min_cq):
                    fails.append("low_consensus_bq")

            for r in fails:
                reasons[r] += 1
            passed = not fails
            n_pass += passed
            if not passed and not a.all:
                continue

            n_ns = sum(c[1] == "nonsupport" for c in classified)
            extra = [
                f"FAM_N={len(reads)},{sup_hq},{sup_lq},{n_ns},{len(noncover)}",
                f"FAM_CQ={cq['F1R2']:.2f},{cq['F2R1']:.2f}",
            ]
            if dissent:
                extra.insert(1, "FAM_NSBQ=" + ",".join(
                    "." if bq is None else str(bq) for _, _, bq, _ in dissent))
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
