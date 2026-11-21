"""SBS96 cosine similarities using ../Figure2_response.ipynb's multinomial method.
Reference comparisons resample the observed counts only; strata comparisons
independently resample both count vectors and compare replicate pairs.
Seed fixed for reproducibility; reference notebook did not fix a seed.
"""
from pathlib import Path
import csv
import numpy as np
import openpyxl

ROOT = Path(__file__).resolve().parent
SHEET = openpyxl.load_workbook(ROOT / 'DupCaller_SBS96_ID83_spectra.xlsx', data_only=True)['SBS96']
rows = list(SHEET.values)
headers = rows[1]
channels = [r[0] for r in rows[2:]]
assert len(channels) == len(set(channels)) == 96
vectors = {name: np.array([r[i] for r in rows[2:]], dtype=float) for i, name in enumerate(headers) if i and name}
rng = np.random.RandomState(20261119)

def compare(a_name, b_name, repeats, both=False):
    a, b = vectors[a_name], vectors[b_name]
    assert np.isfinite(a).all() and np.isfinite(b).all()
    assert (a >= 0).all() and (b >= 0).all() and a.sum() > 0 and b.sum() > 0
    assert np.array_equal(a, a.astype(np.int64))
    if both:
        assert np.array_equal(b, b.astype(np.int64))
    observed = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    # Batches bound memory while retaining the notebook's resampling distribution.
    sims = []
    for start in range(0, repeats, 10000):
        n = min(10000, repeats-start)
        ra = rng.multinomial(int(a.sum()), a/a.sum(), size=n)
        rb = rng.multinomial(int(b.sum()), b/b.sum(), size=n) if both else b
        dot = (ra*rb).sum(axis=1) if both else ra @ rb
        denom = np.linalg.norm(ra,axis=1)*np.linalg.norm(rb,axis=1) if both else np.linalg.norm(ra,axis=1)*np.linalg.norm(rb)
        sims.extend(dot/denom)
    low, median, high = np.percentile(sims, [2.5,50,97.5])
    return dict(spectrum=a_name, reference=b_name, observed_cosine=observed,
                bootstrap_median=median, lower_2_5=low, upper_97_5=high,
                n_spectrum=int(a.sum()), n_reference=float(b.sum()),
                replicates=repeats, resample_both=both, seed=20261119)

results = []
for name in ['AA_DupCaller_only','AA_NanoSeq_only','AA_both']:
    results.append(compare(name,'COSMIC_v3.5_SBS22a (proportion)',10000))
results.append(compare('Sperm_DupCaller','Parental_de_novo_reference (decode trios)',300000))
for group in ['AA','Neuron','Sperm']:
    results.append(compare(group+'_2x2+',group+'_2x2-',100000,both=True))
out = ROOT / 'SBS96_cosine_similarity.csv'
with out.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(results[0]));writer.writeheader();writer.writerows(results)
for r in results:
    print(f"{r['spectrum']} vs {r['reference']}: observed={r['observed_cosine']:.6f}; median={r['bootstrap_median']:.6f}; 95% interval=[{r['lower_2_5']:.6f}, {r['upper_97_5']:.6f}]")
print('Saved:',out)
