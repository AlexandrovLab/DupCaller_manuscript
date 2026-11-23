#!/usr/bin/env python3
"""Tabulate NanoSeq-specific (role = ns_only) sites by cohort and DupCaller --rescue FILTER.

Reads site_classification.csv (written by classify.py) and writes
  NS_only_rescue_flag_SBS.tsv    (category snv)
  NS_only_rescue_flag_indel.tsv  (category indel)
with columns cohort, rescue_flag, number; rows ordered by cohort, then count
(descending), then flag name.

This last counting step was not saved as a script in the original run
(ns_only_rescue_investigation_20261117_dmgfix); this script was written afterwards
and reproduces both tables in DupCaller_figure_reproduction_1.2.13/data/NS_only_rescue/
exactly (same cohorts, flags and counts).

Usage: summarize_rescue_flags.py [site_classification.csv] [out_dir]
"""
import csv
import os
import sys
from collections import Counter

src = sys.argv[1] if len(sys.argv) > 1 else "site_classification.csv"
out_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(src))

rows = [r for r in csv.DictReader(open(src)) if r["role"] == "ns_only"]
for category, tag in [("snv", "SBS"), ("indel", "indel")]:
    counts = Counter((r["cohort"], r["rescue_filter"]) for r in rows if r["category"] == category)
    path = os.path.join(out_dir, f"NS_only_rescue_flag_{tag}.tsv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cohort", "rescue_flag", "number"])
        for (cohort, flag), n in sorted(counts.items(), key=lambda kv: (kv[0][0], -kv[1], kv[0][1])):
            w.writerow([cohort, flag, n])
    print(f"wrote {path} ({sum(counts.values())} {category} sites)")
