"""Per-channel comparison of dmg.tn profiles: new (_dmgfix) vs earlier (_latest, results_20250110).
usage: compare_dmg_tn.py OUTDIR SAMPLE new.txt old_label old.txt [old_label old.txt ...]
Writes OUTDIR/SAMPLE.dmg_compare.tsv (192 channels) and prints per-class pooled summaries."""
import sys, numpy as np, pandas as pd
COMP = str.maketrans("ACGT", "TGCA"); rc = lambda s: s.translate(COMP)[::-1]
SUBS = [("C","A"),("C","G"),("C","T"),("T","A"),("T","C"),("T","G")]; B = "ACGT"
def table(path):
    df = pd.read_csv(path, sep="\t", index_col=0); rows = []
    for strand, pur in (("pyri", False), ("puri", True)):
        for ref, alt in SUBS:
            for l in B:
                for r in B:
                    ctx = l+ref+r; row, a, rb = (rc(ctx), rc(alt), rc(ref)) if pur else (ctx, alt, ref)
                    rows.append((strand, f"{ref}>{alt}", ctx, df.loc[row, a], df.loc[row, rb]))
    return pd.DataFrame(rows, columns=["strand","sub","ctx","alt","opp"])
out, sample, newp, rest = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
new = table(newp); res = new.rename(columns={"alt":"alt_new","opp":"opp_new"}); res["rate_new"] = res.alt_new/res.opp_new
for lab, p in zip(rest[::2], rest[1::2]):
    o = table(p); res[f"alt_{lab}"] = o.alt; res[f"opp_{lab}"] = o.opp; res[f"rate_{lab}"] = o.alt/o.opp
    res[f"ratio_vs_{lab}"] = res.rate_new/res[f"rate_{lab}"]
    print(f"\n== {sample}: new / {lab} ==")
    print(f"opportunity ratio  C-centre(pyri rows) {res.opp_new[(res.strand=='pyri')&res['sub'].str[0].eq('C')].sum()/res[f'opp_{lab}'][(res.strand=='pyri')&res['sub'].str[0].eq('C')].sum():.3f}"
          f"  T-centre {res.opp_new[(res.strand=='pyri')&res['sub'].str[0].eq('T')].sum()/res[f'opp_{lab}'][(res.strand=='pyri')&res['sub'].str[0].eq('T')].sum():.3f}"
          f"  G-centre(puri rows) {res.opp_new[(res.strand=='puri')&res['sub'].str[0].eq('C')].sum()/res[f'opp_{lab}'][(res.strand=='puri')&res['sub'].str[0].eq('C')].sum():.3f}"
          f"  A-centre {res.opp_new[(res.strand=='puri')&res['sub'].str[0].eq('T')].sum()/res[f'opp_{lab}'][(res.strand=='puri')&res['sub'].str[0].eq('T')].sum():.3f}")
    g = res.groupby(["strand","sub"], sort=False)
    s = pd.DataFrame({"alt_new": g.alt_new.sum(), f"alt_{lab}": g[f"alt_{lab}"].sum(),
                      "rate_new": g.alt_new.sum()/g.opp_new.sum(), f"rate_{lab}": g[f"alt_{lab}"].sum()/g[f"opp_{lab}"].sum()})
    s["pooled_ratio"] = s.rate_new/s[f"rate_{lab}"]
    ok = res[(res.alt_new >= 20) & (res[f"alt_{lab}"] >= 20)]
    s["median_ch_ratio(alt>=20)"] = ok.groupby(["strand","sub"], sort=False)[f"ratio_vs_{lab}"].median()
    s["n_ch(alt>=20)"] = ok.groupby(["strand","sub"], sort=False).size()
    s["rate_corr_ch"] = g.apply(lambda d: np.corrcoef(d.rate_new, d[f"rate_{lab}"])[0,1] if d.alt_new.sum()>50 else np.nan)
    pd.set_option("display.width", 250); print(s.to_string(float_format=lambda x: f"{x:.3g}"))
res.to_csv(f"{out}/{sample}.dmg_compare.tsv", sep="\t", index=False, float_format="%.4g")
