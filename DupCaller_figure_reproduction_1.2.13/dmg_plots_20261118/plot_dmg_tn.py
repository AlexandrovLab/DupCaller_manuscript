"""Plot DupCaller *.dmg.tn.txt as SBS96-style rate bars (pyrimidine- and purine-centre panels).
rate = count(alt) / count(centre base) in the same row; purine rows are reverse-complemented
(e.g. G>T on the read strand -> C>A panel), matching results_20250110/*/dmg.{pyri,puri}.png.
Y axis is auto-scaled per panel (max * 1.15).
usage: plot_dmg_tn.py OUTDIR LABEL dmg.tn.txt [LABEL dmg.tn.txt ...]
"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]
SUBS = [("C", "A"), ("C", "G"), ("C", "T"), ("T", "A"), ("T", "C"), ("T", "G")]
COLS = ["#03BCEE", "#010101", "#E32926", "#CAC9C9", "#A1CE63", "#EBC6C4"]
B = "ACGT"

def rates(df, purine):
    out, lab = [], []
    for ref, alt in SUBS:
        for l in B:
            for r in B:
                ctx = l + ref + r
                if purine:
                    row, a, rb = rc(ctx), rc(alt), rc(ref)
                else:
                    row, a, rb = ctx, alt, ref
                out.append(df.loc[row, a] / df.loc[row, rb])
                lab.append(ctx)
    return np.array(out), lab

def panel(ax, y, lab, ymax, title):
    x = np.arange(96)
    ax.bar(x, y, color=np.repeat(COLS, 16), width=0.8)
    ax.set_xlim(-4.5, 99.5); ax.set_ylim(0, ymax)
    ax.set_xticks(x); ax.set_xticklabels(lab, rotation=90, fontsize=11, fontweight="bold", family="monospace")
    ax.tick_params(axis="y", labelsize=22)
    ax.yaxis.get_offset_text().set_fontsize(18)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    for i, ((r, a), c) in enumerate(zip(SUBS, COLS)):
        ax.add_patch(Rectangle((i * 16 - 0.5, 1.0), 16, 0.13, transform=ax.get_xaxis_transform(),
                               color=c, clip_on=False))
        ax.text(i * 16 + 7.5, 1.065, f"{r}>{a}", transform=ax.get_xaxis_transform(), ha="center",
                va="center", fontsize=28, color="white" if c != "#CAC9C9" else "white")
    ax.set_title(title, fontsize=20, loc="left", pad=48)

def main():
    outdir, args = sys.argv[1], sys.argv[2:]
    for lab_s, path in zip(args[::2], args[1::2]):
        df = pd.read_csv(path, sep="\t", index_col=0)
        both = {}
        for name, pur in (("pyri", False), ("puri", True)):
            y, lab = rates(df, pur)
            both[name] = (y, lab)
            ymax = y.max() * 1.15
            fig, ax = plt.subplots(figsize=(20, 4.2))
            panel(ax, y, lab, ymax, f"{lab_s}  dmg.{name}  (rate per centre base; ymax {ymax:.2e})")
            fig.tight_layout(); fig.savefig(f"{outdir}/{lab_s}.dmg.{name}.png", dpi=200); plt.close(fig)
            pd.DataFrame({"context": lab, "sub": np.repeat([f"{r}>{a}" for r, a in SUBS], 16), "rate": y}
                         ).to_csv(f"{outdir}/{lab_s}.dmg.{name}.tsv", sep="\t", index=False)
        fig, axs = plt.subplots(2, 1, figsize=(20, 8.6))
        for ax, name in zip(axs, ("pyri", "puri")):
            y, lab = both[name]
            panel(ax, y, lab, y.max() * 1.15, f"{lab_s}  dmg.{name}")
        fig.tight_layout(); fig.savefig(f"{outdir}/{lab_s}.dmg.png", dpi=200); plt.close(fig)
        print(lab_s, {k: f"{v[0].max():.2e} @ {v[1][v[0].argmax()]}" for k, v in both.items()})

main()
