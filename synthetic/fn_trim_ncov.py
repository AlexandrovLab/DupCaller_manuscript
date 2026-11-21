#!/usr/bin/env python3
"""Diagnose why mock false negatives were missed by DupCaller.

For every record in a false-negative VCF (e.g. compare_*/DC_*_compare_truth/
<sample>/0001.vcf) this checks the two masks below, using the same rules as
DupCaller call (code_rerun15_20261010):

  trim   every alt-supporting duplex family has the variant inside the
         read/fragment end trim (--trimF/--trimR, default 7 bp), as computed
         by funcs/misc.py:determineTrimLength + trim_mask
  n_cov  matched-normal depth < --minNdepth (default 10), using the same
         `samtools mpileup -Q minBq -q mapq` depth as
         funcs/depth.py:extractDepthRegion (dups/anomalous pairs excluded,
         as mpileup does by default)

Trim rule per read, window = the read's aligned reference span
[reference_start, reference_end):
  tlen > 0 (left mate):  left  = trimF
                         right = max(trimR, trimF - min(trimF, |tlen| - reflen))
  otherwise (right mate, processed after its mate):
                         overlap = max(0, reflen + mate_reflen - |tlen| - trimR)
                         left  = max(trimR, overlap, trimF - min(trimF, |tlen| - reflen))
                         right = trimF
  a position at offset o is trimmed if o < left or (reflen - 1 - o) < right.
(DupCaller applies this to the family window using its longest read; the
per-read version here matches it whenever a family's reads share coordinates,
as they do in these mocks.)

Duplex family = reads with the same (chrom, fragment start, fragment end,
sorted DB UMI pair); strand = (is_read1 != is_reverse). A family supports
the variant (untrimmed) when both strands have >=1 alt read (untrimmed).
Reads must pass MAPQ >= --mapq and not be secondary/supplementary/QC-fail/
unmapped; duplicate-flagged reads are kept since DupCaller uses them as
family members. SNV alt bases need BQ >= --minBq. Indels match on type and
length within --indelWindow bp of the truth position (alignment can shift
them inside repeats).

reason column (';'-joined):
  n_cov        normal depth < minNdepth
  trim         alt duplex families exist, but every one is trimmed
  no_duplex    no duplex family carries the alt on both strands
  other        none of the above (the miss has another cause)

usage: fn_trim_ncov.py -b mock.bam -v 0001.vcf -n normal.bam [normal.bam ...] -o out.tsv
"""
import argparse
import sys
from collections import defaultdict

import pysam

CIGAR_M, CIGAR_I, CIGAR_D, CIGAR_EQ, CIGAR_X = 0, 1, 2, 7, 8
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x800


def cigar_ref_len(cigar):
    n, num = 0, ""
    for c in cigar:
        if c.isdigit():
            num += c
        else:
            if c in "MDN=X":
                n += int(num)
            num = ""
    return n


def trim_lengths(read, trimF, trimR):
    tlen = abs(read.template_length)
    reflen = read.reference_length
    frag_short = trimF - min(trimF, max(tlen - reflen, 0)) if tlen else 0
    if read.template_length > 0:
        return trimF, max(frag_short, trimR)
    overlap = 0
    if read.has_tag("MC") and tlen:
        overlap = max(0, reflen + cigar_ref_len(read.get_tag("MC")) - tlen - trimR)
    return max(frag_short, overlap, trimR), trimF


def fragment_coords(read):
    tlen = abs(read.template_length)
    if read.template_length > 0:
        return read.reference_start, read.reference_start + tlen
    if read.template_length < 0:
        return read.reference_end - tlen, read.reference_end
    return read.reference_start, read.reference_end


def family_key(read):
    fs, fe = fragment_coords(read)
    umi = read.get_tag("DB") if read.has_tag("DB") else ""
    return (fs, fe, tuple(sorted(umi.split("-"))))


def find_indel(read, kind, length, center, window):
    """Return the 0-based ref position of a matching indel in this read, or None.
    Insertions are placed at the base after the anchor, deletions at their
    first deleted base."""
    best = None
    rpos = read.reference_start
    for op, n in read.cigartuples:
        if op == CIGAR_I:
            if kind == "I" and n == length and abs(rpos - center) <= window:
                if best is None or abs(rpos - center) < abs(best - center):
                    best = rpos
        elif op == CIGAR_D:
            if kind == "D" and n == length and abs(rpos - center) <= window:
                if best is None or abs(rpos - center) < abs(best - center):
                    best = rpos
            rpos += n
        elif op in (CIGAR_M, CIGAR_EQ, CIGAR_X, 3):
            rpos += n
    return best


def snv_alt(read, pos0, alt, min_bq):
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if rpos == pos0:
            if read.query_qualities[qpos] < min_bq:
                return False
            return read.query_sequence[qpos].upper() == alt
        if rpos > pos0:
            break
    return False


def normal_depth(normals, chrom, start0, end0, min_bq, mapq, max_depth):
    """Min over [start0, end0] of the summed mpileup depth (DupCaller's
    extractDepthRegion). Positions with no pileup line count as 0."""
    depth = defaultdict(int)
    for nb in normals:
        out = pysam.mpileup("-Q", str(min_bq), "-q", str(mapq), "-d", str(max_depth),
                            "-r", f"{chrom}:{start0 + 1}-{end0 + 1}", nb)
        for line in out.split("\n"):
            if line:
                p = line.split("\t")
                depth[int(p[1]) - 1] += int(p[3])
    return min(depth.get(i, 0) for i in range(start0, end0 + 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-b", "--bam", required=True, help="mock tumor bam")
    ap.add_argument("-v", "--vcf", required=True, help="false-negative vcf (e.g. 0001.vcf)")
    ap.add_argument("-n", "--normalBams", nargs="+", required=True, help="matched-normal bam(s)")
    ap.add_argument("-o", "--out", default="-", help="per-site tsv (default stdout)")
    ap.add_argument("-s", "--summary", help="optional summary tsv of reason counts")
    ap.add_argument("--sample", default="", help="sample label for the summary")
    ap.add_argument("--trimF", type=int, default=7)
    ap.add_argument("--trimR", type=int, default=7)
    ap.add_argument("--minNdepth", type=int, default=10)
    ap.add_argument("--mapq", type=int, default=40)
    ap.add_argument("--minBq", type=int, default=18)
    ap.add_argument("--maxPileupDepth", type=int, default=1000000)
    ap.add_argument("--indelWindow", type=int, default=10)
    a = ap.parse_args()

    bam = pysam.AlignmentFile(a.bam)
    cols = ["chrom", "pos", "ref", "alt", "type",
            "tumor_reads", "alt_reads", "alt_reads_untrimmed",
            "alt_frags", "alt_frags_untrimmed",
            "alt_duplex_fam", "alt_duplex_fam_untrimmed",
            "min_dist_read_end", "min_dist_frag_end", "max_dist_read_end_alt",
            "normal_depth", "reason"]
    out = sys.stdout if a.out == "-" else open(a.out, "w")
    out.write("\t".join(cols) + "\n")
    counts = defaultdict(int)
    n_sites = 0

    for rec in pysam.VariantFile(a.vcf):
        for alt in rec.alts:
            n_sites += 1
            ref = rec.ref
            chrom = rec.chrom
            if len(ref) == len(alt) == 1:
                vtype, center = "SNV", rec.pos - 1
            elif len(alt) > len(ref):
                vtype, center = "INS", rec.pos      # base after the anchor
            else:
                vtype, center = "DEL", rec.pos      # first deleted base
            ilen = abs(len(alt) - len(ref))
            fetch_s = max(0, center - (a.indelWindow if vtype != "SNV" else 0))
            fetch_e = center + 1 + (a.indelWindow if vtype != "SNV" else 0)

            tumor_reads = alt_reads = alt_untrim = 0
            frag_any, frag_untrim = set(), set()
            fam_strands = defaultdict(set)
            fam_strands_untrim = defaultdict(set)
            dists_read, dists_frag = [], []
            seen = set()
            for read in bam.fetch(chrom, fetch_s, fetch_e):
                if read.flag & SKIP_FLAGS or read.mapping_quality < a.mapq or read.mate_is_unmapped:
                    continue
                key = (read.query_name, read.is_read1)
                if key in seen:
                    continue
                seen.add(key)
                if vtype == "SNV":
                    if not (read.reference_start <= center < read.reference_end):
                        continue
                    tumor_reads += 1
                    if not snv_alt(read, center, alt.upper(), a.minBq):
                        continue
                    p = center
                else:
                    if read.reference_start <= center < read.reference_end:
                        tumor_reads += 1
                    p = find_indel(read, vtype[0], ilen, center, a.indelWindow)
                    if p is None:
                        continue
                alt_reads += 1
                off = p - read.reference_start
                reflen = read.reference_length
                left, right = trim_lengths(read, a.trimF, a.trimR)
                trimmed = off < left or (reflen - 1 - off) < right
                fs, fe = fragment_coords(read)
                dists_read.append(min(off, reflen - 1 - off))
                dists_frag.append(min(p - fs, fe - 1 - p))
                fam = family_key(read)
                strand = read.is_read1 != read.is_reverse
                frag_any.add(read.query_name)
                fam_strands[fam].add(strand)
                if not trimmed:
                    alt_untrim += 1
                    frag_untrim.add(read.query_name)
                    fam_strands_untrim[fam].add(strand)

            dup = sum(len(s) == 2 for s in fam_strands.values())
            dup_un = sum(len(s) == 2 for s in fam_strands_untrim.values())
            nd = normal_depth(a.normalBams, chrom, rec.pos - 1, rec.pos - 1 + len(ref) - 1 + (1 if vtype == "INS" else 0),
                              a.minBq, a.mapq, a.maxPileupDepth)
            reasons = []
            if nd < a.minNdepth:
                reasons.append("n_cov")
            if dup == 0:
                reasons.append("no_duplex")
            elif dup_un == 0:
                reasons.append("trim")
            reason = ";".join(reasons) if reasons else "other"
            counts[reason] += 1
            row = [chrom, rec.pos, ref, alt, vtype, tumor_reads, alt_reads, alt_untrim,
                   len(frag_any), len(frag_untrim), dup, dup_un,
                   min(dists_read) if dists_read else "NA",
                   min(dists_frag) if dists_frag else "NA",
                   max(dists_read) if dists_read else "NA",
                   nd, reason]
            out.write("\t".join(map(str, row)) + "\n")
    if out is not sys.stdout:
        out.close()

    lines = [f"{r}\t{c}\t{c / n_sites:.3f}" for r, c in sorted(counts.items(), key=lambda x: -x[1])] if n_sites else []
    sys.stderr.write(f"{a.sample or a.vcf}: {n_sites} FN\n" + "\n".join("  " + l for l in lines) + "\n")
    if a.summary:
        with open(a.summary, "w") as f:
            f.write("sample\treason\tn\tfraction\n")
            for l in lines:
                f.write(f"{a.sample}\t{l}\n")


if __name__ == "__main__":
    main()
