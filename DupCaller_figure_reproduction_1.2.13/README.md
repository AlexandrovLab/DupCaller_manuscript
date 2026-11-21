# DupCaller manuscript figures — reproduction package

This folder is self-contained: it holds every input needed to regenerate the Python-made
panels of the DupCaller manuscript (main_v8 Figures 2–6 and supplementary_figures_v8 S1–S6),
the code, and the generated plots.

```
DupCaller_figure_reproduction/
├── README.md
├── requirements.txt / environment.yml   Python dependencies
├── DupCaller_manuscript_plots.py         script: regenerates every plot
├── DupCaller_manuscript_plots.ipynb      same code as a notebook, one section per figure
├── data/                                 all inputs (read-only)
└── DupCaller_manuscript_plots/           outputs (PNG, 600 dpi, 8-pt text)
    └── Supplementary/
```

## How to run

```bash
cd DupCaller_figure_reproduction
pip install -r requirements.txt          # or: conda env create -f environment.yml && conda activate dupcaller-figures
python DupCaller_manuscript_plots.py     # ~1 min; rewrites DupCaller_manuscript_plots/
```

Or open `DupCaller_manuscript_plots.ipynb` in Jupyter and choose *Run All*.
Paths are relative to this folder, so it can be moved or zipped as a whole.
Tested with Python 3.13 using both pandas 2.2 / matplotlib 3.9 and pandas 3.0 / matplotlib 3.11 (seaborn 0.13.2).

## Inputs (`data/`)

| File | Contents | Used for |
|---|---|---|
| `DupCaller_figure_data.xlsx` | Per-sample data, one tab per analysis (stored as values): `synthetic_SBS_summary`, `synthetic_indel_summary` (precision/recall/F1 per synthetic dataset), `AA_neuron_cordblood` (shared/caller-specific counts, burdens, duplication rate per sample), `sperm_blood_saliva` (per-sample burdens and duplication rate), `runtime` (hours per sample, step and pipeline), `downsampling` (one row per downsampled dataset), `cordblood_CI_calculation` (Sample, Mutation Type, number, correction ratio, Caller — SBS96 calls and correction ratios per cord-blood sample for DupCaller and NanoSeq Software; Fig 4d) | all main figures, S2, S4, S6 |
| `synthetic_SBS_benchmark.xlsx` | `README` tab; `NanoSeq_Software_v3.5.5` (TP/FP/FN/precision/recall per dataset); one DupCaller precision–recall curve per dataset, tabs `SBS_n<N>_rep<R>` | Fig S1 |
| `synthetic_indel_benchmark.xlsx` | same layout for indels (tabs `indel_n<N>_rep<R>`) | Fig S3 |
| `NS_only_rescue/NS_only_rescue_SBS.tsv`, `NS_only_rescue_indel.tsv` | NanoSeq-specific calls per cohort, by the DupCaller step that excluded them (`cohort`, `rescue_flag`, `number`); not plotted by the script | Results/Methods "Analysis of NanoSeq-specific mutations", Response R1.3 |
| `source_data_spectra/DupCaller_SBS96_ID83_spectra.xlsx` | SBS96 / ID83 counts behind the spectrum panels (spectra not plotted by the script; the SBS96 tab feeds the cosine-similarity table) | Figs 3b,d,e,g · 4c · 5c · 6c,d · S5 · Response Fig 1.3 · `SBS96_cosine_similarity.csv` |

Tab ↔ previous supplementary-table numbering: `synthetic_SBS_summary` = Supp Table 1 (v5), `synthetic_indel_summary` = 2, `AA_neuron_cordblood` = 3, `sperm_blood_saliva` = 4, `runtime` = 5, `downsampling` = 6.
Synthetic dataset names map to the original pipeline IDs as `SBS_n100_rep1` = `mock_100snp_1`, `indel_n10_rep1` = `mock_10indel_1`.

## Outputs (`DupCaller_manuscript_plots/`)

| File | Panel | Data |
|---|---|---|
| `Figure2b_precision.png`, `Figure2b_recall.png` | Fig 2b | `downsampling` |
| `Figure2c.png` | Fig 2c | `downsampling` |
| `Figure2d.png` | Fig 2d | `downsampling` |
| `Figure3c_SBS_count.png`, `Figure3c_SBS_burden.png` | Fig 3c | `AA_neuron_cordblood` |
| `Figure3f_indel_count.png`, `Figure3f_indel_burden.png` | Fig 3f | `AA_neuron_cordblood` |
| `Figure4a.png`, `Figure4b.png`, `Figure4d.png` | Fig 4a, b, d | `AA_neuron_cordblood` |
| `Figure5a.png`, `Figure5b.png`, `Figure5d.png` | Fig 5a, b, d | `AA_neuron_cordblood` |
| `Figure6a.png`, `Figure6b.png` | Fig 6a, b | `sperm_blood_saliva` |
| `Supplementary/FigureS1.png` (+ `FigureS1_S3_legend.svg`) | Fig S1 | `synthetic_SBS_benchmark.xlsx` |
| `Supplementary/FigureS2b_{precision,recall,F1}.png` | Fig S2b (SBS) | `synthetic_SBS_summary` |
| `Supplementary/FigureS2c_{precision,recall,F1}.png` | Fig S2c (indel) | `synthetic_indel_summary` |
| `Supplementary/FigureS3.png` | Fig S3 | `synthetic_indel_benchmark.xlsx` |
| `Supplementary/FigureS4.png` | Fig S4 | `AA_neuron_cordblood` |
| `Supplementary/FigureS6.png` | Fig S6 | `runtime` |
| `SBS96_cosine_similarity.csv` | SBS96 cosine similarities with 95% multinomial-resampling intervals: AA calls vs SBS22a, sperm vs parental de novo reference, 2×2+ vs 2×2− (AA, neuron, sperm) | `source_data_spectra/DupCaller_SBS96_ID83_spectra.xlsx` (SBS96) |

Not made by this code: schematics and illustrations (Figs 1, 2a, 3a, S2a, S7) and the mutational-spectrum
panels (Figs 3b/d/e/g, 4c, 5c, 6c/d, S5); the counts for the spectra are in `data/source_data_spectra/`.

## Notes

- Figure 4d: for both DupCaller and NanoSeq Software, the corrected SBS burden and its 95% Fay–Feuer (gamma) interval are computed from `cordblood_CI_calculation` (pooled over the six PD48442 samples, divided by duplex coverage, ×1/0.76 sensitivity factor); only the standard-WGS reference values (Abascal et al.) are entered in the script: SBS 1.18 × 10⁻⁸ (1.10–1.27 × 10⁻⁸) for Fig 4d and indel 6.5 × 10⁻¹⁰ (4.8–8.7 × 10⁻¹⁰) for Fig 5d. Figure 5d uses Poisson intervals computed from `AA_neuron_cordblood`.
- Fig S4's call-out text is generated from the data (the highlighted neuron sample, PD43289, has a duplication rate of 75.4%).
- Error bars in Figs 2b–d and S6 are mean ± SD (deterministic; matches the figure legends).
- `SBS96_cosine_similarity.csv` uses the same method, comparisons and seed (20261119) as `data/source_data_spectra/calculate_cosine_similarity.py` and reproduces `data/source_data_spectra/SBS96_cosine_similarity.csv` (differences ≤ 1e-15). The resampling adds ~20 s to a run.
- Fig S6 is drawn from the `runtime` tab; the NanoSeq "total" bar there (~21 h) differs from the v8 PDF panel (~14.7 h).
