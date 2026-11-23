#!/usr/bin/env python3
"""Label every site in sites.tsv (NS_only + shared controls) with DupCaller's
FILTER from (a) the genome-wide <sample>_dmgfix run (PASS + _fail vcfs, no
--rescue) and (b) this directory's region-restricted --rescue run.

Matching: SNVs look in SBS + atomized DBS records, indels in INDEL records
(exact allele, else same pos/same length, else same length within +-10bp).
When a site has several DupCaller read-family records, the one whose
TAG1/TAG2 barcode pair equals NanoSeq's read-bundle barcode (SNV BTAG,
indel RB) is used (match=family); otherwise PASS if any, else the most common
label (match=site, i.e. DupCaller only saw the alt in a different family)."""
import csv, os, re
from collections import Counter, defaultdict

ROOT = "/tscc/lustre/restricted/alexandrov-ddn/users/yuc211/NanoSeq_Sanger_Analysis/runtime_bechmarking"
OUT = os.path.join(ROOT, "ns_only_rescue_investigation_20261117_dmgfix")
AA = {"40mM-72h", "70mM-48h", "80mM-48h"}
WINDOW = 10


def cohort(s):
    return "Aristolochic Acid" if s in AA else ("Cord blood" if s.startswith(("PD47269", "PD48442")) else "Neuron")


def info_dict(info):
    return dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in info.split(";"))


def load_vcf(path, kind, out, wanted):
    """Append records near a wanted (chrom,pos) to out[(chrom,pos,ref,alt)];
    DBS atomized to two SNVs."""
    if not os.path.isfile(path):
        return False
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 2)
            if (f[0], int(f[1])) not in wanted:
                continue
            f = line.rstrip("\n").split("\t")
            c, p, r, a, flt, info = f[0], int(f[1]), f[3], f[4], f[6], info_dict(f[7])
            rec = {"filter": flt, "kind": kind, "info": info,
                   "tags": frozenset((info.get("TAG1"), info.get("TAG2"))),
                   "tumor": f[9] if len(f) > 9 else "", "normal": f[10] if len(f) > 10 else "",
                   "allele": f"{p}:{r}>{a}"}
            if kind == "dbs":
                for i in range(len(r)):
                    if r[i] != a[i]:
                        out[(c, p + i, r[i], a[i])].append(rec)
            else:
                out[(c, p, r, a)].append(rec)
    return True


def load_run(prefix_dir, tag, wanted):
    snv, indel = defaultdict(list), defaultdict(list)
    ok = True
    for sub, kind, store in [("SBS", "sbs", snv), ("DBS", "dbs", snv), ("INDEL", "indel", indel)]:
        for suf in ["", "_fail"]:
            ok &= load_vcf(os.path.join(prefix_dir, sub, f"{tag}_{kind}{suf}.vcf"), kind, store, wanted)
    return snv, indel, ok


def nas_tags(cat, info):
    d = info_dict(info)
    if cat == "snv" and "BTAG" in d:
        return frozenset(d["BTAG"].split("|"))
    if "RB" in d:
        rb = d["RB"].split(",")
        return frozenset(rb[3:5])
    return None


def find(cat, c, p, r, a, snv, indel):
    if cat == "snv":
        return snv.get((c, p, r, a), []), "exact"
    if (c, p, r, a) in indel:
        return indel[(c, p, r, a)], "exact"
    dlen = len(a) - len(r)
    near = sorted((abs(k[1] - p), k) for k in indel
                  if k[0] == c and abs(k[1] - p) <= WINDOW and len(k[3]) - len(k[2]) == dlen)
    if near:
        return indel[near[0][1]], f"shifted_{near[0][0]}bp"
    return [], ""


def pick(recs, tags):
    fam = [x for x in recs if tags and x["tags"] == tags]
    pool, how = (fam, "family") if fam else (recs, "site")
    best = next((x for x in pool if x["filter"] == "PASS"), None)
    if best is None:
        top = Counter(x["filter"] for x in pool).most_common(1)[0][0]
        best = next(x for x in pool if x["filter"] == top)
    return best, how


def label(cat, c, p, r, a, tags, run):
    snv, indel, ok = run
    if not ok:
        return {"filter": "run_missing"}
    recs, allele_how = find(cat, c, p, r, a, snv, indel)
    if not recs:
        return {"filter": "not_found"}
    best, fam_how = pick(recs, tags)
    i = best["info"]
    return {"filter": best["filter"], "match": f"{fam_how}/{allele_how}", "kind": best["kind"],
            "allele": best["allele"], "n_families": len(recs),
            "all_filters": ";".join(f"{k}x{n}" for k, n in Counter(x["filter"] for x in recs).most_common()),
            "F1R2": i.get("F1R2", ""), "F2R1": i.get("F2R1", ""), "LR": i.get("LR", ""),
            "FDR": i.get("FDR", ""), "TC": i.get("TC", ""), "BC": i.get("BC", ""),
            "MSP": i.get("MSP", ""), "SNPM": i.get("SNPM", ""), "NOISEM": i.get("NOISEM", ""),
            "HP": i.get("HP", ""), "STR": i.get("STR", ""),
            "tumor_AC_RC_DP": best["tumor"], "normal_AC_RC_DP": best["normal"]}


def main():
    sites = list(csv.DictReader(open(os.path.join(OUT, "sites.tsv")), delimiter="\t"))
    # (chrom,pos) keys any site can match: +-WINDOW for indel shifts, +-1 so an
    # SNV still sees a DBS whose first base is the preceding position
    wanted = defaultdict(set)
    for st in sites:
        c, p = st["chrom"], int(st["pos"])
        for d in range(-WINDOW, WINDOW + 1):
            wanted[st["sample"]].add((c, p + d))
    runs = {}
    rows = []
    for st in sites:
        s = st["sample"]
        if s not in runs:
            runs.clear()
            runs[s] = (load_run(os.path.join(ROOT, f"{s}_dmgfix"), f"{s}_dmgfix", wanted[s]),
                       load_run(os.path.join(OUT, f"rescue_{s}", f"rescue_{s}"), f"rescue_{s}", wanted[s]))
        latest, rescue = runs[s]
        c, p, r, a, cat = st["chrom"], int(st["pos"]), st["ref"], st["alt"], st["category"]
        tags = nas_tags(cat, st["nas_info"])
        lt = label(cat, c, p, r, a, tags, latest)
        rs = label(cat, c, p, r, a, tags, rescue)
        nd = info_dict(st["nas_info"])
        row = {"sample": s, "cohort": cohort(s), "role": st["role"], "category": cat,
               "chrom": c, "pos": p, "ref": r, "alt": a,
               "nas_tags": "|".join(sorted(tags)) if tags else "",
               "nas_depth_fwd": nd.get("DEPTH_FWD", ""), "nas_depth_rev": nd.get("DEPTH_REV", ""),
               "nas_norm_fwd": nd.get("DEPTH_NORM_FWD", ""), "nas_norm_rev": nd.get("DEPTH_NORM_REV", ""),
               "nas_qpos": nd.get("QPOS", ""), "nas_dplx_nm": nd.get("DPLX_NM", ""),
               "nas_dplx_asxs": nd.get("DPLX_ASXS", ""), "nas_bulk_asxs": nd.get("BULK_ASXS", ""),
               "latest_filter": lt["filter"], "latest_match": lt.get("match", ""),
               "latest_all_filters": lt.get("all_filters", "")}
        row.update({f"rescue_{k}": v for k, v in rs.items()})
        rows.append(row)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    path = os.path.join(OUT, "site_classification.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")

    for role in ["control", "ns_only"]:
        for cat in ["snv", "indel"]:
            sub = [r for r in rows if r["role"] == role and r["category"] == cat]
            print(f"\n=== {role} {cat} (n={len(sub)}) ===")
            print("  latest :", dict(Counter(r["latest_filter"] for r in sub).most_common()))
            print("  rescue :", dict(Counter(r["rescue_filter"] for r in sub).most_common()))
            print("  match  :", dict(Counter(r.get("rescue_match", "") for r in sub).most_common()))


if __name__ == "__main__":
    main()
