from pysam import VariantFile as VCF
import sys
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
tp_vcf = VCF(sys.argv[1],'r')
fp_vcf = VCF(sys.argv[2],'r')
fn_vcf = VCF(sys.argv[3],'r')
cuts = set()
fn_num = 0
for rec in fn_vcf:
    fn_num += 1
for rec in tp_vcf:
    GTmin = rec.info["FDR"]
    cuts.update((GTmin,))
for rec in fp_vcf:
    GTmin = rec.info["FDR"]
    cuts.update((GTmin,)) 
cuts = np.array(sorted(list(cuts)))
tp_mat = np.zeros(cuts.size)
fp_mat = np.zeros(cuts.size)
#fn_mat = np.ones(cuts.size)*fn_num
tp_vcf = VCF(sys.argv[1],'r')
fp_vcf = VCF(sys.argv[2],'r')
for rec in tp_vcf:
    GTmin = rec.info["FDR"]
    tp_mat[np.searchsorted(cuts,GTmin)] += 1
for rec in fp_vcf:
    GTmin = rec.info["FDR"]
    fp_mat[np.searchsorted(cuts,GTmin)] += 1
tp_mat = np.cumsum(tp_mat)
fp_mat = np.cumsum(fp_mat)
fn_num = fn_num + max(tp_mat)
fn_mat = fn_num - tp_mat
prec = tp_mat/(tp_mat+fp_mat)
rec = tp_mat/(tp_mat+fn_mat)
f1 = 2* tp_mat / (2*tp_mat + fp_mat + fn_mat)
#thr_ind = np.searchsorted(cuts,np.arange(0,10,0.1))
thr_ind = np.searchsorted(cuts,0.05)
#plt.plot(rec[thr_ind[0]:thr_ind[-1]],prec[thr_ind[0]:thr_ind[-1]])
plt.plot(rec,prec)
plt.ylabel("preicision")
plt.xlabel("recall")
plt.ylim(0, 1.05)
plt.xlim(0, 0.8)
print(cuts[thr_ind-1:thr_ind+1],prec[thr_ind-1],rec[thr_ind-1],tp_mat[thr_ind-1],fp_mat[thr_ind-1],fn_mat[thr_ind-1])
output_pd = pd.DataFrame({"precision":prec.copy(),"recall":rec.copy(),"tp":tp_mat.copy(),"fp":fp_mat.copy(),"fn":fn_mat.copy(),"thr":cuts.copy()})
#print(output_pd)
output_pd.to_csv(sys.argv[6],sep='\t',index=False)
#plt.scatter(rec[thr_ind[0:35]],prec[thr_ind[0:35]],s=50,marker="o")
#print(prec,rec,tp_mat,fp_mat,cuts)
#plt.savefig(sys.argv[4],dpi=300)
plt.figure()
plt.plot(cuts,f1)
plt.xlabel("FDR threshold")
plt.ylabel("F1 score")
plt.savefig(sys.argv[5],dpi=300)
