#!/usr/bin/env python3
"""Why did NanoSeq (NAS) not call the neuron DC_only 2x2+ SBS of compare_20261010?

For each site in snv_intersect_2x2+/<s>/0000.vcf, dsa_rows/<s>/<chr>_<pos>.bed
holds NAS's own dsa pileup rows (run_dsa_site.sh, exact args of the real run).
The row of the DC read family (dsa barcode == DC TAG1|TAG2) is pushed through a
line-by-line port of variantcaller.cc (CollectMetrics/CallDuplex/ApplyFilters/
VafFilter/IsMasked) with the run's thresholds
(-a 50 -b 0 -c 0.02 -d 2 -f 0.9 -i 1.0 -m 8 -n 3 -p 0 -q 60 -r 144 -v 0.01 -x 8 -z 10).

Consensus failures asked about:
  frac_lt_0.9 : on F1R2 or F2R1 no base reaches 90% of the canonical reads
                (>=10% of the strand's bases disagree) -> CallDuplex gives 'N' and
                variantcaller drops the row before any filter.
  bq_lt_60    : consensus reached on both strands, but the summed base quality
                of the called base is < 60 on F1R2 or F2R1
                (ConsensusBaseQualityFilter).
All other variantcaller filters are evaluated too (even for frac-failed rows) so
we can say whether the consensus failure is the sole reason.
"""
import csv, glob, os
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ISEC = os.path.join(os.path.dirname(HERE), "snv_intersect_2x2+")
A, BULK, Z, CLIP, D, FRAC, INDEL, MINC, NMMS, PPAIR, Q, RLEN, VAF, MAXC = \
    50, 0, 10, 0.02, 2, 0.9, 1.0, 8, 3, 0, 60, 144, 0.01, 8
ALPH = "ACGT"
FILTERS = ["dplx_clip", "alignment_score", "mismatch", "matched_normal", "duplex",
           "consensus_bq_lt60", "indel", "five_prime_trim", "three_prime_trim", "proper_pair", "vaf"]


def parse(f):
    r = {"beg": int(f[1]), "ctx": f[3], "snp": int(f[4]), "sw": int(f[5]),
         "basxs": float(f[6]), "bnm": float(f[7]),
         "bf": [int(x) for x in f[8:13]], "br": [int(x) for x in f[13:18]],
         "bp_beg": int(f[18]), "bp_end": int(f[19]), "bc": f[20], "type": int(f[21]),
         "dasxs": float(f[22]), "dclip": float(f[23]), "dnm": float(f[24]),
         "f1": [int(x) for x in f[25:30]], "f2": [int(x) for x in f[30:35]],
         "q1": [int(x) for x in f[35:39]], "q2": [int(x) for x in f[39:43]],
         "bpp": float(f[43]), "dpp": float(f[44])}
    r["left"] = r["beg"] - r["bp_beg"] + 1
    r["right"] = r["bp_end"] - r["beg"] - 1
    return r


def strand_call(counts):
    can = sum(counts[:4])
    call = "N"
    for i in range(4):
        if can and counts[i] / can >= FRAC:
            call = ALPH[i]
    return call


def evaluate(r, alt):
    """Port of variantcaller.cc for one dsa row. Returns (stage, failed filter list).
    For rows without consensus, filters are computed with the DC alt as the call."""
    c1, c2 = strand_call(r["f1"]), strand_call(r["f2"])
    if sum(r["f1"][:4]) == 0 or sum(r["f2"][:4]) == 0:
        stage = "strand_no_canonical_base"
    elif c1 == "N" or c2 == "N":
        stage = "frac_lt_0.9"
    elif c1 != c2:
        stage = "strands_disagree"
    elif c1 == r["ctx"][1]:
        stage = "consensus_is_ref"
    elif c1 != alt:
        stage = "consensus_other_alt"
    else:
        stage = "consensus_alt"
    call = alt if stage != "consensus_alt" else c1
    i = ALPH.index(call)
    isvar = r["ctx"][1] != call
    bft, brt = sum(r["bf"]), sum(r["br"])
    f1t, f2t = sum(r["f1"]), sum(r["f2"])
    bfc, brc = sum(r["bf"][:4]), sum(r["br"][:4])
    mx = RLEN - MAXC
    sides = {1: [r["left"]], 2: [r["right"]], 3: [r["left"], r["right"]]}[r["type"]]
    ok = {
        "dplx_clip": r["dclip"] <= CLIP,
        "alignment_score": r["basxs"] >= A and r["dasxs"] >= A,
        "mismatch": r["bnm"] <= NMMS and r["dnm"] <= NMMS + (1 if isvar else 0),
        "matched_normal": bft >= BULK and brt >= BULK and bft + brt >= Z,
        "duplex": f1t >= D and f2t >= D,
        "consensus_bq_lt60": r["q1"][i] >= Q and r["q2"][i] >= Q,
        "indel": not any(t > 0 and ind / t > INDEL for ind, t in
                         [(r["bf"][4], bft), (r["br"][4], brt), (r["f1"][4], f1t), (r["f2"][4], f2t)]),
        "five_prime_trim": all(x >= MINC for x in sides),
        "three_prime_trim": all(x <= mx for x in sides),
        "proper_pair": r["bpp"] >= PPAIR and r["dpp"] >= PPAIR,
        "vaf": not ((bfc > 0 and r["bf"][i] / bfc > VAF) or (brc > 0 and r["br"][i] / brc > VAF)),
    }
    failed = [k for k in FILTERS if not ok[k]]
    masked = r["snp"] + r["sw"] > 0
    return stage, failed, masked, c1, c2


def info(s):
    return dict(kv.split("=", 1) for kv in s.split(";") if "=" in kv)


rows = []
for vcf in sorted(glob.glob(os.path.join(ISEC, "PD432*", "0000.vcf"))):
    s = vcf.split("/")[-2]
    for l in open(vcf):
        if l[0] == "#":
            continue
        x = l.rstrip("\n").split("\t")
        c, p, ref, alt = x[0], int(x[1]), x[3], x[4]
        inf = info(x[7])
        dcb = {inf["TAG1"], inf["TAG2"]}
        dsa = [parse(f.split("\t")) for f in open(os.path.join(HERE, "dsa_rows", s, f"{c}_{p}.bed")).read().splitlines() if f]
        dsa = [r for r in dsa if r["beg"] == p - 1]
        fam = [r for r in dsa if set(r["bc"].split("|")) == dcb]
        # does any family at the site pass NAS completely with this alt? (sanity: should be ~0)
        any_pass = any(e[0] == "consensus_alt" and not e[1] and not e[2]
                       for e in (evaluate(r, alt) for r in dsa))
        out = {"sample": s, "chrom": c, "pos": p, "ref": ref, "alt": alt,
               "dc_F1R2": inf["F1R2"], "dc_F2R1": inf["F2R1"], "dc_TC": inf.get("TC", ""), "dc_BC": inf.get("BC", ""),
               "dc_barcodes": "|".join(sorted(dcb)), "n_dsa_families": len(dsa), "other_family_passes": any_pass}
        if not fam:
            out.update(stage="dc_family_absent_from_dsa", failed_filters="", masked="", detail="")
        else:
            r = fam[0]
            stage, failed, masked, c1, c2 = evaluate(r, alt)
            out.update(stage=stage, failed_filters=",".join(failed), masked=int(masked),
                       detail=(f"f1r2 ACGTI={r['f1']} Q={r['q1']} call={c1}; f2r1 ACGTI={r['f2']} Q={r['q2']} call={c2}; "
                               f"normal fwd={r['bf']} rev={r['br']}; ASXS={r['dasxs']}/{r['basxs']} NM={r['dnm']}/{r['bnm']} "
                               f"clip={r['dclip']} left={r['left']} right={r['right']} type={r['type']}"))
        rows.append(out)

with open(os.path.join(HERE, "DC_only_SBS_2x2plus_neuron_NAS_reason.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

n = len(rows)
print(f"neuron DC_only 2x2+ SBS: n={n}")
print("stage of the DC family's dsa row:")
for k, v in Counter(r["stage"] for r in rows).most_common():
    print(f"  {k:28s} {v:5d}")
frac = [r for r in rows if r["stage"] == "frac_lt_0.9"]
bq = [r for r in rows if r["stage"] == "consensus_alt" and "consensus_bq_lt60" in r["failed_filters"].split(",")]
print(f"\n>=10% of a strand's bases disagree (frac<0.9): {len(frac)}")
print(f"   of which would also fail >=1 other filter: {sum(1 for r in frac if r['failed_filters'] or r['masked'])}")
print(f"   of which ALSO have summed BQ<60 on a strand: {sum(1 for r in frac if 'consensus_bq_lt60' in r['failed_filters'])}")
print(f"consensus reached, summed BQ<60 on a strand:    {len(bq)}")
print(f"   BQ<60 is the ONLY failing filter:           {sum(1 for r in bq if r['failed_filters'] == 'consensus_bq_lt60' and not r['masked'])}")
print(f"union (frac<0.9 OR BQ<60):                      {len(frac) + len(bq)}")
ca = [r for r in rows if r["stage"] == "consensus_alt"]
print("\nconsensus_alt rows, failing filters (counted per filter):")
for k, v in Counter(f for r in ca for f in r["failed_filters"].split(",") if f).most_common():
    print(f"  {k:28s} {v:5d}")
print(f"  masked (snp/noise)           {sum(1 for r in ca if r['masked'])}")
print(f"  passes everything            {sum(1 for r in ca if not r['failed_filters'] and not r['masked'])}")
print(f"\nsites where SOME family passes NAS fully: {sum(1 for r in rows if r['other_family_passes'])}")
