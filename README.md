# DupCaller manuscript: analysis scripts

Scripts behind the DupCaller manuscript analyses. Python 3 with pysam, pandas, numpy and matplotlib
(see `DupCaller_figure_reproduction_1.2.13/requirements.txt` for the plotting package).

```
.
├── downsampling/                        read-family downsampling and family-level truth filters
├── synthetic/                           synthetic (spike-in) NanoSeq benchmark
├── nanoseq_comparison/                  why NanoSeq missed DupCaller-only calls
└── DupCaller_figure_reproduction_1.2.13/  self-contained figure reproduction package (own README)
```

## downsampling/

| Script | What it does | Input | Output |
| --- | --- | --- | --- |
| `randomly_downsample_read_family.py <in.bam> <out_prefix> <seed>` | Groups reads into duplex families (barcode, start, template length). Each family with ≥ 3 reads on both strands is downsampled to five F1R2 × F2R1 compositions; the reads left out are flagged as duplicates. Smaller families are copied unchanged to all five BAMs. | DupCaller-ready BAM (barcodes in the read name) | `<out_prefix>_{1x1,1x2,2x1,1xn,nx1}.bam` |
| `family_consensus_filter.py -i <in.vcf> -b <in.bam> -o <out.vcf> [-f ref.fa] [--min-cq 60] [--all]` | Rebuilds the read family behind each DupCaller call the way DupCaller groups reads. A call passes only when every read in its family supports the variant and the consensus base quality on each strand is > 60. | DupCaller SBS/indel VCF, the BAM it was called from, reference FASTA (needed for indels) | VCF of passing calls (`--all`: every call, annotated) |
| `family_bq_filter.py -i <in.vcf> -b <in.bam> -o <out.vcf> [-f ref.fa] [--min-bq 18] [--min-support 3] [--max-nonsupport 1] [--min-cq 60]` | BQ-aware version of the family filter, used to build the high-confidence truth set. A call passes with ≥ 3 supporting reads at BQ ≥ 18, at most one non-supporting read (and only if its BQ < 18), and per-strand consensus BQ > 60. Imports `family_consensus_filter.py` from the same directory. | DupCaller SBS/indel VCF, BAM, reference FASTA (indels) | VCF of passing calls |

## synthetic/

| Script | What it does | Input | Output |
| --- | --- | --- | --- |
| `filter_variants_by_normal_depth.py <in.vcf> <normal.bam> <sample.bam> <out.vcf>` | Picks candidate spike-in sites. Keeps biallelic SNVs and indels with depth ≥ 10 and no alt reads in donor B, and ≥ 1 alt read in donor A (MAPQ > 30, BQ ≥ 30). | Candidate VCF, donor B BAM, donor A BAM | VCF of sites present in A and absent from B |
| `generate_synthetic_NanoSeq.py <A_R1> <A_R2> <B_R1> <B_R2> <A.bam> <B.bam> <candidates.vcf> <unused> <truth.vcf> <out_R1> <out_R2> <n_mutations> <B_undiluted.bam>` | Builds a synthetic NanoSeq sample. For each of `n_mutations` random candidate sites, it takes one duplex family carrying the variant from donor A and removes the same number of reads covering the site from donor B. Indel sites with any indel read in `B_undiluted.bam` are skipped. Records the family size in INFO `RF`, `RT`, `RB`. | Donor A/B FASTQ pairs and BAMs, candidate VCF (from `filter_variants_by_normal_depth.py`), undiluted donor B BAM | Mixed FASTQ pair (donor B minus removed reads, plus donor A families), truth VCF of the spiked-in variants |
| `recall_prec_curve.py <TP.vcf> <FP.vcf> <FN.vcf> <unused> <F1.png> <prec_rec.tsv>` | Scores DupCaller calls against the truth across INFO/FDR thresholds. | TP/FP/FN VCFs (e.g. from `bcftools isec`); TP and FP need INFO/FDR | F1-vs-FDR-threshold plot; table of precision, recall, TP, FP, FN per threshold; stdout: metrics at FDR 0.05 (saved as `<sample>.recall.log`) |
| `fn_trim_ncov.py -b <mock.bam> -v <FN.vcf> -n <normal.bam> [...] -o <out.tsv> [-s summary.tsv]` | Explains why each false negative was missed, using DupCaller's own rules. Reasons: `n_cov` (normal depth < 10), `trim` (every alt duplex family falls in the read-end trim), `no_duplex` (no family has the alt on both strands), or `other`. | Synthetic sample BAM, false-negative VCF, matched-normal BAM(s) | Per-site TSV with a `reason` column; optional summary TSV of reason counts |
| `nanoseq_indel_build_chart.py <benchmark_dir> > chart.tsv` | Builds the DupCaller vs NanoSeq indel benchmark table for the 25 `mock_{10,50,100,200,500}indel_{1..5}` samples. Reports TP/FP/FN, precision, recall and F1, raw and with FPs corrected for expected true somatic mutations (duplex coverage × 1.18 × 10⁻⁸ × 0.76). | Under `<benchmark_dir>`: `compare_20260910/logs/<s>.recall.log`, `DupCaller_new_indels/.../<s>_latest_sbs_burden.txt` (duplex coverage), `NAS_indel_z10_rerun_20260911/compare/truth/<s>/counts.tsv`, `NAS/sams/<s>/tmpNanoSeq/post/results.mut_burden.tsv` | TSV on stdout, one row per sample |

## nanoseq_comparison/

Both scripts read their inputs relative to their own location (sibling directories one level up, and subdirectories next to the script) and write the CSV next to the script.

| Script | What it does | Input | Output |
| --- | --- | --- | --- |
| `classify_sbs.py` | For each DupCaller-only 2×2+ SBS in neurons, passes the DupCaller read family's NanoSeq `dsa` row through a port of NanoSeq's `variantcaller.cc` filters, using the run's thresholds. Reports the step at which NanoSeq loses the call (e.g. strand consensus < 90%, consensus BQ < 60) and every other filter it fails. | `../snv_intersect_2x2+/PD432*/0000.vcf`, `dsa_rows/<sample>/<chr>_<pos>.bed` | `DC_only_SBS_2x2plus_neuron_NAS_reason.csv`; stdout: summary counts by stage and filter |
| `classify_indel.py` | Same analysis for DupCaller-only 2×2+ indels (neuron and cord blood): follows the DupCaller family through NanoSeq indel calling steps 1–3 and the final normalised set, and reports the first step where it is lost. | `../indel_intersect_2x2+/PD4*/0000.vcf`, `dsa_rows/`, `nas_indel_stages/<sample>.{step1.bed,step2.vcf,step3.vcf}`, `../NAS_indel/norm/noChrY/<sample>.vcf.gz` | `DC_only_indel_2x2plus_NAS_reason.csv`; stdout: summary counts per group by stage and sub-reason |

## DupCaller_figure_reproduction_1.2.13/

Self-contained: inputs, code and generated plots. See the folder's `README.md` for the panel-by-panel list.

| Script | What it does | Input | Output |
| --- | --- | --- | --- |
| `DupCaller_manuscript_plots.py` (or `.ipynb`) | Regenerates every Python-made panel of main Figures 2–6 and supplementary Figures S1–S4 and S6, plus the SBS96 cosine-similarity table. | `data/DupCaller_figure_data.xlsx`, `data/synthetic_{SBS,indel}_benchmark.xlsx`, `data/source_data_spectra/DupCaller_SBS96_ID83_spectra.xlsx` | `DupCaller_manuscript_plots/*.png`, `DupCaller_manuscript_plots/Supplementary/*.png`, `DupCaller_manuscript_plots/SBS96_cosine_similarity.csv` |
| `data/source_data_spectra/calculate_cosine_similarity.py` | Standalone SBS96 cosine similarities with 95% multinomial-resampling intervals (seed 20261119). Compares AA calls with SBS22a, sperm with the parental de novo reference, and 2×2+ with 2×2− for AA, neuron and sperm. | `DupCaller_SBS96_ID83_spectra.xlsx` (SBS96 sheet), next to the script | `SBS96_cosine_similarity.csv` next to the script; stdout: one line per comparison |
| `dmg_plots_20261118/plot_dmg_tn.py <outdir> <label> <dmg.tn.txt> [<label> <dmg.tn.txt> ...]` | Plots DupCaller's damage profile as SBS96-style rate bars (alt count / centre-base count), with pyrimidine-centre and purine-centre (reverse-complemented) panels. | One or more DupCaller `*.dmg.tn.txt` | `<label>.dmg.png` (both panels), `<label>.dmg.{pyri,puri}.png`, `<label>.dmg.{pyri,puri}.tsv` |
| `dmg_plots_20261118/compare_dmg_tn.py <outdir> <sample> <new.txt> <label> <old.txt> [...]` | Compares a new damage profile with one or more earlier runs, channel by channel (192 channels). | Two or more `*.dmg.tn.txt` | `<sample>.dmg_compare.tsv` (per-channel counts, rates and ratios); stdout: opportunity ratios and pooled per-class summary |
