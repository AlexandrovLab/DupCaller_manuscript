# DupCaller manuscript: analysis scripts

Scripts behind the DupCaller manuscript analyses. Python 3 with pysam, pandas, numpy and matplotlib
(see `DupCaller_figure_reproduction_1.2.13/requirements.txt` for the plotting package).

```
.
├── downsampling/                        read-family downsampling and family-level truth filters
├── synthetic/                           synthetic (spike-in) NanoSeq benchmark
├── nanoseq_comparison/                  why NanoSeq missed DupCaller-only calls, and why DupCaller missed NanoSeq-only calls
└── DupCaller_figure_reproduction_1.2.13/  self-contained figure reproduction package (own README)
```

## downsampling/

| Script | What it does | Input | Output |
| --- | --- | --- | --- |
| `randomly_downsample_read_family.py <in.bam> <out_prefix> <seed>` | Groups reads into duplex families (barcode, start, template length). Each family with ≥ 3 reads on both strands and a non-duplicate (primary) read is downsampled to five F1R2 × F2R1 compositions; only the sampled reads are written (reads left out are dropped, not flagged), with their original duplicate flags, and if no sampled read is primary one is picked at random and made primary. Families below 3×3 that have a primary read contribute only that one non-duplicate read (plus its mate) to all five BAMs; their other reads are dropped. Families with no primary read are written unchanged (all reads still duplicate-flagged) to all five BAMs. Secondary, supplementary, QC-fail, improperly paired and `DT`-tagged reads, and reads whose 5′ end is not an aligned match (e.g. soft-clipped), are skipped. | DupCaller-ready BAM (barcodes in the read name) | `<out_prefix>_{1x1,1x2,2x1,1xn,nx1}.bam` |
| `family_consensus_filter.py -i <in.vcf> -b <in.bam> -o <out.vcf> [-f ref.fa] [--min-cq 60] [--all]` | Rebuilds the read family behind each DupCaller call the way DupCaller groups reads. A call passes only when every read in its family supports the variant and the consensus base quality on each strand is > 60. | DupCaller SBS/indel VCF, the BAM it was called from, reference FASTA (needed for indels) | VCF of passing calls (`--all`: every call, annotated) |
| `family_bq_filter.py -i <in.vcf> -b <in.bam> -o <out.vcf> [-f ref.fa] [--min-bq 18] [--max-nonsupport 1] [--min-cq 60]` | BQ-aware version of the family filter, used to build the high-confidence truth set. A call passes with at most one non-supporting read (and only if its BQ < 18) and per-strand consensus BQ > 60. Imports `family_consensus_filter.py` from the same directory. | DupCaller SBS/indel VCF, BAM, reference FASTA (indels) | VCF of passing calls |

The downsampled BAMs are not written in coordinate order (each read is written when its family is processed, and its
mate later at the mate's position), so sort and index them before running DupCaller:

```bash
python downsampling/randomly_downsample_read_family.py <in.bam> <out_prefix> <seed>
for c in 1x1 1x2 2x1 1xn nx1; do
    samtools sort -@ 8 -o <out_prefix>_${c}.sorted.bam <out_prefix>_${c}.bam
    mv <out_prefix>_${c}.sorted.bam <out_prefix>_${c}.bam
    samtools index <out_prefix>_${c}.bam
done
```

## synthetic/

| Script | What it does | Input | Output |
| --- | --- | --- | --- |
| `filter_variants_by_normal_depth.py <in.vcf> <normal.bam> <sample.bam> <out.vcf>` | Picks candidate spike-in sites. Keeps biallelic SNVs and indels with depth ≥ 10 and no alt reads in donor B, and ≥ 1 alt read in donor A (MAPQ > 30, BQ ≥ 30). | Candidate VCF, donor B BAM, donor A BAM | VCF of sites present in A and absent from B |
| `generate_synthetic_NanoSeq.py <A_R1> <A_R2> <B_R1> <B_R2> <A.bam> <B.bam> <candidates.vcf> <unused> <truth.vcf> <out_R1> <out_R2> <n_mutations> <B_undiluted.bam>` | Builds a synthetic NanoSeq sample. For each of `n_mutations` random candidate sites, it takes one duplex family carrying the variant from donor A and removes the same number of reads covering the site from donor B. Indel sites with any indel read in `B_undiluted.bam` are skipped. Records the family size in INFO `RF`, `RT`, `RB`. | Donor A/B FASTQ pairs and BAMs, candidate VCF (from `filter_variants_by_normal_depth.py`), undiluted donor B BAM | Mixed FASTQ pair (donor B minus removed reads, plus donor A families), truth VCF of the spiked-in variants |
| `recall_prec_curve.py <TP.vcf> <FP.vcf> <FN.vcf> <unused> <F1.png> <prec_rec.tsv>` | Scores DupCaller calls against the truth across INFO/FDR thresholds. | TP/FP/FN VCFs (e.g. from `bcftools isec`); TP and FP need INFO/FDR | F1-vs-FDR-threshold plot; table of precision, recall, TP, FP, FN per threshold; stdout: `[thr, next thr] precision recall TP FP FN` at the last threshold with FDR ≤ 0.05, inclusive (saved as `<sample>.recall.log`; if no call passes, precision is `nan` and TP = FP = 0) |
| `fn_trim_ncov.py -b <mock.bam> -v <FN.vcf> -n <normal.bam> [...] -o <out.tsv> [-s summary.tsv]` | Explains why each false negative was missed, using DupCaller's own rules. Reasons: `n_cov` (normal depth < 10), `trim` (every alt duplex family falls in the read-end trim), `no_duplex` (no family has the alt on both strands), or `other`. | Synthetic sample BAM, false-negative VCF, matched-normal BAM(s) | Per-site TSV with a `reason` column; optional summary TSV of reason counts |

## nanoseq_comparison/

### DupCaller-only calls (why NanoSeq missed them)

Both scripts read their inputs relative to their own location (sibling directories one level up, and subdirectories next to the script) and write the CSV next to the script.

| Script | What it does | Input | Output |
| --- | --- | --- | --- |
| `classify_sbs.py` | For each DupCaller-only 2×2+ SBS in neurons, passes the DupCaller read family's NanoSeq `dsa` row through a port of NanoSeq's `variantcaller.cc` filters, using the run's thresholds. Reports the step at which NanoSeq loses the call (e.g. strand consensus < 90%, consensus BQ < 60) and every other filter it fails. | `../snv_intersect_2x2+/PD432*/0000.vcf`, `dsa_rows/<sample>/<chr>_<pos>.bed` | `DC_only_SBS_2x2plus_neuron_NAS_reason.csv`; stdout: summary counts by stage and filter |
| `classify_indel.py` | Same analysis for DupCaller-only 2×2+ indels (neuron and cord blood): follows the DupCaller family through NanoSeq indel calling steps 1–3 and the final normalised set, and reports the first step where it is lost. | `../indel_intersect_2x2+/PD4*/0000.vcf`, `dsa_rows/`, `nas_indel_stages/<sample>.{step1.bed,step2.vcf,step3.vcf}`, `../NAS_indel/norm/noChrY/<sample>.vcf.gz` | `DC_only_indel_2x2plus_NAS_reason.csv`; stdout: summary counts per group by stage and sub-reason |

### `ns_only_rescue/`: NanoSeq-specific calls (why DupCaller missed them)

Produces `DupCaller_figure_reproduction_1.2.13/data/NS_only_rescue/NS_only_rescue_{SBS,indel}.tsv`. Each NanoSeq-only site
is labelled with the FILTER DupCaller gives it when re-run with `--rescue` (which reports every candidate with its
failing filter instead of dropping it) over ±2 kb around the site. `build_regions.py`, `run_rescue.py`, `classify.py`
and `samples.txt` are copied unchanged from the original run
(`$RES_HOME/NanoSeq_Sanger_Analysis/runtime_bechmarking/ns_only_rescue_investigation_20261117_dmgfix`), so their paths are
absolute paths on TSCC.
`summarize_rescue_flags.py` was written afterwards for the final counting step (not saved as a script in the original
run); it reproduces both TSVs byte for byte from that run's `site_classification.csv`.

| Step | Script | What it does | Output |
| --- | --- | --- | --- |
| 1 | `build_regions.py` | Takes the NanoSeq-only SNVs and indels (`compare_20261117_dmgfix/{snv,indel}_intersect/<sample>/0001.vcf`) plus up to 30 shared sites per type (`0002.vcf`, seed 20261010) as positive controls, and merges ±2 kb around them. | `<sample>_regions.bed.gz`, `sites.tsv` |
| 2 | `run_rescue.py [sample ...]` | Runs every sample that has a `<sample>_regions.bed.gz` (or only the samples given; `samples.txt` lists the 26 used), up to 8 at a time with 4 threads each. Runs the same DupCaller code as the genome-wide `<sample>_dmgfix` calls (1.2.12, main 64dac19 + damage uint8 fix) with `--rescue -R <regions>`, the same SNP/noise masks, and that run's error profiles (`-E`) and mutation-rate tables (`-mr`). Unlike the genome-wide runs, `-nms` is left off so masked candidates are still reported with their FILTER. | `rescue_<sample>/` (SBS/DBS/INDEL VCFs incl. `_fail`) |
| 3 | `classify.py` | Matches each site to DupCaller records (SNVs: SBS and atomised DBS records; indels: exact allele, else same length within ±10 bp). If several read families carry it, uses the one whose TAG1/TAG2 barcodes equal NanoSeq's read bundle (`BTAG`/`RB`), else PASS if any, else the most common FILTER; `not_found` if there is no record. All shared controls should come back PASS. | `site_classification.csv` |
| 4 | `summarize_rescue_flags.py [site_classification.csv] [out_dir]` | Counts NanoSeq-only sites by cohort and rescue FILTER, separately for SBS and indels. | `NS_only_rescue_flag_{SBS,indel}.tsv` (saved in the figure package as `data/NS_only_rescue/NS_only_rescue_{SBS,indel}.tsv`) |

Rescue flags: `cov_mask` matched-normal depth < 10; `normal_vaf` matched normal has the alt above `--naf` (1%);
`high_nm` read family blacklisted because a read (usually the mate) has NM ≥ 5; `nm_mask` position-level NM or AS−XS mask;
`low_ASXS` / `low_mapq` family-level alignment filters; `trim_mask` read-end trim zone; `masked` noise mask;
`underpowered` below the FDR/LR threshold; `strand_independence`; `no_good_alt_read` no usable alt read in the pileup;
`not_found` no DupCaller record at the site even under `--rescue`.

## DupCaller_figure_reproduction_1.2.13/

Self-contained: inputs, code and generated plots. See the folder's `README.md` for the panel-by-panel list.

| Script | What it does | Input | Output |
| --- | --- | --- | --- |
| `DupCaller_manuscript_plots.py` (or `.ipynb`) | Regenerates every Python-made panel of main Figures 2–6 and supplementary Figures S1–S4 and S6, plus the SBS96 cosine-similarity table. | `data/DupCaller_figure_data.xlsx`, `data/synthetic_{SBS,indel}_benchmark.xlsx`, `data/source_data_spectra/DupCaller_SBS96_ID83_spectra.xlsx` | `DupCaller_manuscript_plots/*.png`, `DupCaller_manuscript_plots/Supplementary/*.png`, `DupCaller_manuscript_plots/SBS96_cosine_similarity.csv` |
| `data/source_data_spectra/calculate_cosine_similarity.py` | Standalone SBS96 cosine similarities with 95% multinomial-resampling intervals (seed 20261119). Compares AA calls with SBS22a, sperm with the parental de novo reference, and 2×2+ with 2×2− for AA, neuron and sperm. | `DupCaller_SBS96_ID83_spectra.xlsx` (SBS96 sheet), next to the script | `SBS96_cosine_similarity.csv` next to the script; stdout: one line per comparison |
