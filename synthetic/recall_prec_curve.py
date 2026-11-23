from pysam import VariantFile as VCF
import sys
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

DEFAULT_FDR = 0.05   # DupCaller default operating point: calls with lFDR <= 0.05

def fdr_values(path):
    return np.array([rec.info["FDR"] for rec in VCF(path, 'r')], dtype=float)

tp_fdr = np.sort(fdr_values(sys.argv[1]))
fp_fdr = np.sort(fdr_values(sys.argv[2]))
fn_num = sum(1 for _ in VCF(sys.argv[3], 'r'))

# One threshold per distinct FDR value. Row i counts the calls with FDR <= cuts[i]
# (side="right" keeps calls whose FDR equals the threshold).
cuts = np.unique(np.concatenate([tp_fdr, fp_fdr]))
tp_mat = np.searchsorted(tp_fdr, cuts, side="right").astype(float)
fp_mat = np.searchsorted(fp_fdr, cuts, side="right").astype(float)
truth_num = fn_num + tp_fdr.size          # FN outside the call set + every TP call
fn_mat = truth_num - tp_mat
prec = tp_mat/(tp_mat+fp_mat)
rec = tp_mat/(tp_mat+fn_mat)
f1 = 2* tp_mat / (2*tp_mat + fp_mat + fn_mat)

# Default operating point: the last threshold <= DEFAULT_FDR (inclusive).
thr_ind = np.searchsorted(cuts, DEFAULT_FDR, side="right") - 1
plt.plot(rec,prec)
plt.ylabel("precision")
plt.xlabel("recall")
plt.ylim(0, 1.05)
plt.xlim(0, 0.8)
# stdout (saved as <sample>.recall.log): [cut at/below DEFAULT_FDR, next cut] precision recall tp fp fn
if thr_ind < 0:
    # no call has FDR <= DEFAULT_FDR: nothing is called, precision is undefined
    print(f"WARNING: no calls with FDR <= {DEFAULT_FDR}", file=sys.stderr)
    print(cuts[:1], np.nan, 0.0, 0.0, 0.0, float(truth_num))
else:
    print(cuts[thr_ind:thr_ind+2],prec[thr_ind],rec[thr_ind],tp_mat[thr_ind],fp_mat[thr_ind],fn_mat[thr_ind])
output_pd = pd.DataFrame({"precision":prec.copy(),"recall":rec.copy(),"tp":tp_mat.copy(),"fp":fp_mat.copy(),"fn":fn_mat.copy(),"thr":cuts.copy()})
output_pd.to_csv(sys.argv[6],sep='\t',index=False)
plt.figure()
plt.plot(cuts,f1)
plt.xlabel("FDR threshold")
plt.ylabel("F1 score")
plt.savefig(sys.argv[5],dpi=300)
