#!/usr/bin/env python3
import gzip
import re
import subprocess
import sys

# usage: nanoseq_indel_build_chart.py <synthetic_benchmark_dir> > chart.tsv
ROOT = sys.argv[1]
CORR_RATE = 0.0000000118
CORR_FRAC = 0.76

SAMPLES = [f"mock_{g}indel_{i}" for g in [10, 50, 100, 200, 500] for i in range(1, 6)]

def dupcaller_tpfpfn(sample):
    log = f"{ROOT}/compare_20260910/logs/{sample}.recall.log"
    with open(log) as f:
        line = f.read().strip()
    # format: [c1 c2] prec rec tp fp fn
    nums = re.findall(r"[\d.eE+-]+", line.split("]")[1])
    prec, rec, tp, fp, fn = (float(x) for x in nums[:5])
    return tp, fp, fn, prec, rec

def dupcaller_coverage(sample):
    f = f"{ROOT}/DupCaller_new_indels/bams/mkdped/{sample}_latest/SBS/{sample}_latest_sbs_burden.txt"
    with open(f) as fh:
        for line in fh:
            if line.startswith("Duplex coverage"):
                return float(line.split("\t")[1])
    raise RuntimeError(f"no duplex coverage for {sample}")

def nanoseq_tpfpfn(sample):
    f = f"{ROOT}/NAS_indel_z10_rerun_20260911/compare/truth/{sample}/counts.tsv"
    with open(f) as fh:
        _, tp, fp, fn = fh.read().strip().split("\t")
    return int(tp), int(fp), int(fn)

def nanoseq_coverage(sample):
    f = f"{ROOT}/NAS/sams/{sample}/tmpNanoSeq/post/results.mut_burden.tsv"
    with open(f) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "observed":
                return float(parts[2])
    raise RuntimeError(f"no coverage for {sample}")

def mutation_number(sample):
    m = re.match(r"mock_(\d+)indel_\d+", sample)
    return int(m.group(1))

def f1(p, r):
    return 2 * p * r / (p + r) if (p + r) else 0.0

rows = []
for sample in SAMPLES:
    mutnum = mutation_number(sample)

    dtp, dfp, dfn, dprec, drec = dupcaller_tpfpfn(sample)
    dcov = dupcaller_coverage(sample)
    dcorr = dcov * CORR_RATE * CORR_FRAC
    dfp_c = dfp - dcorr
    dprec_c = dtp / (dtp + dfp_c)
    df1 = f1(dprec, drec)
    df1_c = f1(dprec_c, drec)

    ntp, nfp, nfn = nanoseq_tpfpfn(sample)
    ncov = nanoseq_coverage(sample)
    nprec = ntp / (ntp + nfp) if (ntp + nfp) else 0.0
    nrec = ntp / mutnum
    ncorr = ncov * CORR_RATE * CORR_FRAC
    nfp_c = nfp - ncorr
    nprec_c = ntp / (ntp + nfp_c)
    nf1 = f1(nprec, nrec)
    nf1_c = f1(nprec_c, nrec)

    rows.append(dict(
        sample=sample, mutnum=mutnum,
        dtp=dtp, dfp=dfp, dfn=dfn, dprec=dprec, drec=drec, dcov=dcov,
        dcorr=dcorr, dfp_c=dfp_c, dprec_c=dprec_c, df1=df1, df1_c=df1_c,
        ntp=ntp, nfp=nfp, nfn=nfn, nprec=nprec, nrec=nrec, ncov=ncov,
        ncorr=ncorr, nfp_c=nfp_c, nprec_c=nprec_c, nf1=nf1, nf1_c=nf1_c,
    ))

header = ["Sample","DupCaller_True_Positive_Number","DupCaller_False_Positive_Number",
          "DupCaller_False_Negative_Number","DupCaller_Precision","DupCaller_Recall",
          "DupCaller_coverage","DupCaller_False_Positive_Correction",
          "DupCaller_Corrected_False_Positive_Number","DupCaller_Corrected_Precision",
          "DupCaller_F1_score","DupCaller_corrected_F1_score",
          "NanoSeq_True_Positive_Number","NanoSeq_False_Positive_Number",
          "NanoSeq_True_Negative_Number","NanoSeq_Precision","NanoSeq_Recall",
          "NanoSeq_coverage","NanoSeq_False_Positive_Correction",
          "NanoSeq_Corrected_False_Positive_Number","NanoSeq_Corrected_Precision",
          "NanoSeq_F1_score","NanoSeq_corrected_F1_score","Mutation number"]
print("\t".join(header))
for r in rows:
    vals = [
        r["sample"]+"/", r["dtp"], r["dfp"], r["dfn"], r["dprec"], r["drec"], r["dcov"],
        r["dcorr"], r["dfp_c"], r["dprec_c"], r["df1"], r["df1_c"],
        r["ntp"], r["nfp"], r["nfn"], r["nprec"], r["nrec"], r["ncov"],
        r["ncorr"], r["nfp_c"], r["nprec_c"], r["nf1"], r["nf1_c"], r["mutnum"],
    ]
    print("\t".join(str(v) for v in vals))
