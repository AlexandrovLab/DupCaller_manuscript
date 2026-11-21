#!/usr/bin/env python3
"""Why did NanoSeq (NAS) not call the DC_only 2x2+ indels of compare_20261010?

NAS indel calling (tmpNanoSeq_oldmask_20260813_112213/indel, args.json:
--rb 2 --t5 8 --t3 136 -z 10 -v 0.01 -a 50 -c 0.02) is three steps per chunk:
  step1 indelCaller_step1.pl : dsa rows -> rows of read bundles carrying an indel
         (bulk cov >= 10, dplx clip <= 0.02, dplx NM <= 20, AS-XS >= 50 dplx+bulk,
          >= 2 reads per strand, indel in >= 90% of the bundle's reads,
          8 <= qpos <= 136); tags BULK_SEEN if bulk indel reads > 1% of bulk.
  step2 indelCaller_step2.pl : per read bundle, drop BULK_SEEN positions, skip the
         bundle if every remaining position is SNP/noise masked, then
         bcftools mpileup|call on the bundle's reads; keeps ONE record per bundle
         (best-QUAL PASS); FILTER MASKED (>25% of REF span masked) / BULK_SEEN.
  step3 indelCaller_step3.R  : matched normal +-5 bp: NEI_IND if any position has
         indel fraction > 0.01, MISSINGBULK if no bases.
Final NAS set = step3 PASS (post), then norm -D / noChrY / bcftools isec vs DC.

For each DC_only site: dsa_rows/<s>/<chr>_<pos>.bed = NAS's own dsa rows over the
indel span (run_dsa_site.sh, exact args of the real run); nas_indel_stages/<s>.*
= the real run's step1/2/3 outputs. The DC family (dsa barcode == TAG1|TAG2) is
followed through every step and the first step at which it is lost is reported.
"""
import csv, glob, os, re
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ISEC = os.path.join(ROOT, "indel_intersect_2x2+")
STG = os.path.join(HERE, "nas_indel_stages")
RB, T5, T3, MC, VAF, ASXS, CLIP, NM = 2, 8, 136, 10, 0.01, 50, 0.02, 20
STEP1 = ["bulk_cov_lt10", "dplx_clip", "dplx_nm_gt20", "asxs_lt50", "strand_lt2_reads",
         "indel_frac_lt0.9", "qpos_lt8", "qpos_gt136"]
NEURON = {"PD43272", "PD43273", "PD43274", "PD43275", "PD43276", "PD43277", "PD43279", "PD43280",
          "PD43281", "PD43282", "PD43283", "PD43284", "PD43285", "PD43287", "PD43288", "PD43289", "PD43290"}


def parse(f):
    r = {"beg": int(f[1]), "end": int(f[2]), "ctx": f[3], "snp": int(f[4]), "sw": int(f[5]),
         "basxs": float(f[6]), "bnm": float(f[7]),
         "bf": [int(x) for x in f[8:13]], "br": [int(x) for x in f[13:18]],
         "bp_beg": int(f[18]), "bp_end": int(f[19]), "bc": f[20],
         "dasxs": float(f[22]), "dclip": float(f[23]), "dnm": float(f[24]),
         "f1": [int(x) for x in f[25:30]], "f2": [int(x) for x in f[30:35]]}
    return r


def step1(r):
    """Port of indelCaller_step1.pl for one dsa row -> (failed filters, bulk_seen, qpos)."""
    r1, r2 = sum(r["f1"]), sum(r["f2"])
    bt = sum(r["bf"]) + sum(r["br"])
    qf, qr = r["beg"] - r["bp_beg"] + 1, r["bp_end"] - r["end"]
    q = min(qf, qr)
    ok = {
        "bulk_cov_lt10": bt >= MC,
        "dplx_clip": r["dclip"] <= CLIP,
        "dplx_nm_gt20": r["dnm"] <= NM,
        "asxs_lt50": r["dasxs"] >= ASXS and r["basxs"] >= ASXS,
        "strand_lt2_reads": r1 >= RB and r2 >= RB,
        "indel_frac_lt0.9": r["f1"][4] + r["f2"][4] >= 0.9 * (r1 + r2),
        "qpos_lt8": q >= T5,
        "qpos_gt136": q <= T3,
    }
    bulk_seen = r["bf"][4] + r["br"][4] > VAF * bt
    return [k for k in STEP1 if not ok[k]], bulk_seen, q


def info(s):
    return dict(kv.split("=", 1) for kv in s.split(";") if "=" in kv)


def load_stage(s):
    s1 = []   # (chr, pos1, rbkey, line)
    for l in open(os.path.join(STG, f"{s}.step1.bed")):
        x = l.rstrip("\n").split("\t")
        if len(x) < 4:
            continue
        rbid = x[3].split(";")[0]                    # chr:beg-end:BC1|BC2
        c, be, bc = rbid.split(":")
        b, e = be.split("-")
        s1.append((x[0], int(x[2]), (c, b, e, frozenset(bc.split("|"))), x[3]))
    vcf = {}
    for st in ("step2", "step3"):
        recs = []
        for l in open(os.path.join(STG, f"{s}.{st}.vcf")):
            x = l.rstrip("\n").split("\t")
            m = re.search(r"RB=([^;]+)", x[7])
            c, b, e, *bc = m.group(1).split(",")
            recs.append({"chr": x[0], "pos": int(x[1]), "ref": x[3], "alt": x[4], "qual": x[5],
                         "filter": x[6], "rb": (c, b, e, frozenset(bc)),
                         "nn": (re.search(r"NN=(\[[^\]]*\])", x[7]) or [None, ""])[1]})
        vcf[st] = recs
    fin = []
    p = os.path.join(ROOT, "NAS_indel", "norm", "noChrY", f"{s}.vcf.gz")
    for l in os.popen(f"zcat {p}"):
        if l[0] != "#":
            x = l.split("\t")
            fin.append((x[0], int(x[1]), x[3], x[4]))
    return s1, vcf["step2"], vcf["step3"], fin


rows = []
cache = {}
for vcf in sorted(glob.glob(os.path.join(ISEC, "PD4*", "0000.vcf"))):
    s = vcf.split("/")[-2]
    for l in open(vcf):
        if l[0] == "#":
            continue
        if s not in cache:
            cache[s] = load_stage(s)
        s1, s2, s3, fin = cache[s]
        x = l.rstrip("\n").split("\t")
        c, p, ref, alt = x[0], int(x[1]), x[3], x[4]
        inf = info(x[7])
        dcb = frozenset({inf["TAG1"], inf["TAG2"]})
        lo, hi = p - 1, p + len(ref)               # 1-based window around the indel
        out = {"group": "neuron" if s in NEURON else "cordblood", "sample": s, "chrom": c, "pos": p,
               "ref": ref, "alt": alt, "dc_F1R2": inf["F1R2"], "dc_F2R1": inf["F2R1"],
               "dc_TC": inf.get("TC", ""), "dc_BC": inf.get("BC", ""), "dc_HP": inf.get("HP", ""),
               "dc_barcodes": "|".join(sorted(dcb))}
        dsa = [parse(f.split("\t")) for f in open(os.path.join(HERE, "dsa_rows", s, f"{c}_{p}.bed")).read().splitlines() if f]
        fam = [r for r in dsa if frozenset(r["bc"].split("|")) == dcb]
        ind = [r for r in fam if r["f1"][4] + r["f2"][4] > 0]
        ev = [(step1(r), r) for r in ind]
        passing = [(e, r) for e, r in ev if not e[0]]
        fam_s1 = [t for t in s1 if t[0] == c and lo <= t[1] <= hi and t[2][0] == c and t[2][3] == dcb]
        rbkeys = {t[2] for t in fam_s1}
        fam_s2 = [v for v in s2 if v["rb"] in rbkeys]
        fam_s3 = [v for v in s3 if v["rb"] in rbkeys]
        near_fin = [f for f in fin if f[0] == c and abs(f[1] - p) <= 20]
        # best dsa row to report = fewest step1 failures
        best = min(ev, key=lambda t: len(t[0][0])) if ev else None
        detail = ""
        if best:
            (fl, bs, q), r = best
            detail = (f"pos={r['end']} f1r2 ACGTI={r['f1']} f2r1 ACGTI={r['f2']} normal fwd={r['bf']} rev={r['br']} "
                      f"ASXS={r['dasxs']}/{r['basxs']} NM={r['dnm']}/{r['bnm']} clip={r['dclip']} qpos={q} "
                      f"snp={r['snp']} sw={r['sw']}")
        out.update(n_dsa_families=len({r["bc"] for r in dsa}),
                   dsa_indel_rows=len(ind),
                   step1_failed=",".join(best[0][0]) if best and not passing else "",
                   bulk_seen=int(any(e[1] for e, _ in passing)) if passing else "",
                   in_step1=len(fam_s1),
                   step2=";".join(f"{v['pos']}:{v['ref']}>{v['alt']}:{v['filter']}:Q{v['qual']}" for v in fam_s2),
                   step3=";".join(f"{v['pos']}:{v['ref']}>{v['alt']}:{v['filter']}:NN{v['nn']}" for v in fam_s3),
                   nas_final_within20=";".join(f"{f[1]}:{f[2]}>{f[3]}" for f in near_fin),
                   detail=detail)
        same = lambda v: v["pos"] == p and v["ref"] == ref and v["alt"] == alt
        if not fam:
            stage = "dc_family_absent_from_dsa"
        elif not ind:
            stage = "no_indel_in_family_dsa"
        elif not passing:
            stage = "step1_filtered"
        elif not fam_s1:
            stage = "step1_pass_but_absent_from_step1_bed"     # sanity: should not happen
        elif all(e[1] for e, _ in passing):
            stage = "bulk_seen_normal_indel_gt1pct"
        elif not fam_s2:
            stage = "step2_no_call"                            # all masked, or bcftools called nothing
        elif not any(same(v) for v in fam_s2):
            stage = "step2_called_different_indel"
        elif not any(same(v) and v["filter"] == "PASS" for v in fam_s2):
            stage = "step2_" + next(v["filter"] for v in fam_s2 if same(v))
        elif not any(same(v) and v["filter"] == "PASS" for v in fam_s3):
            stage = "step3_" + next(v["filter"] for v in fam_s3 if same(v))
        else:
            stage = "nas_pass_but_not_matched_in_isec"
        sub = ""
        if stage == "step3_NEI_IND":
            nn = next(v["nn"] for v in fam_s3 if same(v)).strip("[]").split(":")
            # step3 compares max_per_site to max_vaf passed as a *string* ("0.01"):
            # a zero-depth position in the +-5 bp normal window gives 0/0 = NaN and
            # "NaN" > "0.01" is TRUE, so the site is flagged with 0 normal indel reads.
            sub = ("NaN_string_compare_bug" if nn[2] == "NaN" else
                   "normal_indel_ge20pct_nearby" if float(nn[2]) >= 0.2 else "normal_indel_lt20pct")
        elif stage == "step2_called_different_indel":
            d = min(abs(v["pos"] - p) for v in fam_s2)
            sub = "same_event_other_representation" if d <= 25 else "other_indel_in_bundle_kept"
        elif stage == "step1_filtered":
            sub = out["step1_failed"]
        out["stage"] = stage
        out["sub_reason"] = sub
        rows.append(out)

fields = list(rows[0])
fields.insert(fields.index("dc_barcodes") + 1, fields.pop(fields.index("stage")))
fields.insert(fields.index("stage") + 1, fields.pop(fields.index("sub_reason")))
with open(os.path.join(HERE, "DC_only_indel_2x2plus_NAS_reason.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)

for g in ("neuron", "cordblood"):
    rs = [r for r in rows if r["group"] == g]
    print(f"\n{g} DC_only 2x2+ indels: n={len(rs)}")
    for k, v in Counter(r["stage"] for r in rs).most_common():
        print(f"  {k:40s} {v:4d}")
    for k, v in Counter((r["stage"], r["sub_reason"]) for r in rs
                        if r["stage"] in ("step3_NEI_IND", "step2_called_different_indel")).most_common():
        print(f"  {k[0]:28s} {k[1]:34s} {v:4d}")
    f1 = [r for r in rs if r["stage"] == "step1_filtered"]
    if f1:
        print("  step1 filters failed (per filter, best dsa row):")
        for k, v in Counter(f for r in f1 for f in r["step1_failed"].split(",") if f).most_common():
            print(f"    {k:38s} {v:4d}")
        print("  step1 combos:", Counter(r["step1_failed"] for r in f1).most_common())
