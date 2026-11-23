#!/usr/bin/env python3
"""Build +-2000bp merged rescue region files per sample from compare_20261117_dmgfix's
NS_only SNV + INDEL sets (snv_intersect/indel_intersect/<sample>/0001.vcf),
plus up to N_CONTROL shared (0002.vcf) SNVs and indels per sample as positive
controls: a region-restricted --rescue run reusing the genome-wide -E/-mr
should re-call every shared site as PASS. Also writes sites.tsv listing every
site to classify (role = ns_only | control)."""
import csv, os, random, subprocess, sys

ROOT = "/tscc/lustre/restricted/alexandrov-ddn/users/yuc211/NanoSeq_Sanger_Analysis/runtime_bechmarking"
CMP = os.path.join(ROOT, "compare_20261117_dmgfix")
OUT_DIR = os.path.join(ROOT, "ns_only_rescue_investigation_20261117_dmgfix")
FLANK = 2000
N_CONTROL = 30
SEED = 20261010


def records(path):
    if not os.path.isfile(path):
        return []
    with open(path) as fh:
        return [l.rstrip("\n").split("\t") for l in fh if l[0] != "#"]


def build_for_sample(s, writer):
    rng = random.Random(f"{SEED}-{s}")
    sites = []
    for cat in ["snv", "indel"]:
        isec = os.path.join(CMP, f"{cat}_intersect", s)
        for f in records(os.path.join(isec, "0001.vcf")):
            sites.append(("ns_only", cat, f))
        shared = records(os.path.join(isec, "0002.vcf"))
        for f in rng.sample(shared, min(N_CONTROL, len(shared))):
            sites.append(("control", cat, f))
    if not any(r == "ns_only" for r, _, _ in sites):
        return None
    for role, cat, f in sites:
        writer.writerow([s, role, cat, f[0], f[1], f[3], f[4], f[7]])
    ivs = sorted((f[0], max(0, int(f[1]) - 1 - FLANK), int(f[1]) + FLANK) for _, _, f in sites)
    merged = []
    for c, st, en in ivs:
        if merged and merged[-1][0] == c and st <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], en)
        else:
            merged.append([c, st, en])
    bed = os.path.join(OUT_DIR, f"{s}_regions.bed")
    with open(bed, "w") as fh:
        for c, st, en in merged:
            fh.write(f"{c}\t{st}\t{en}\n")
    subprocess.run(f"bgzip -f {bed} && tabix -f -p bed {bed}.gz", shell=True, check=True)
    n_ns = sum(r == "ns_only" for r, _, _ in sites)
    return n_ns, len(sites) - n_ns, len(merged)


if __name__ == "__main__":
    samples = sys.argv[1:] or sorted(os.listdir(os.path.join(CMP, "snv_intersect")))
    with open(os.path.join(OUT_DIR, "sites.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "role", "category", "chrom", "pos", "ref", "alt", "nas_info"])
        for s in samples:
            r = build_for_sample(s, w)
            print(f"{s}: no NS_only variants, skipped" if r is None
                  else f"{s}: {r[0]} NS_only + {r[1]} control sites -> {r[2]} merged regions")
