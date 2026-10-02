# Reproducible analysis: ferroptosis defense-axis and m6A readers in lung adenocarcinoma

Scripts, intermediate result tables and figures for the manuscript
**"A reproducible, targetable ferroptosis defense-axis signature in lung adenocarcinoma"**
(working title; manuscript v3.5, bilingual Chinese/English edition).

The repository accompanies the paper's reproducibility commitment: every numeric claim in the
revision modules (subgroup survival, time-dependent AUC, drug sensitivity, CRISPR dependency,
immune infiltration, mutation landscape, GSEA, multivariable Cox, and spatial analyses) is
traceable to a script in this repository plus a table in `data/`.

## Repository layout

| Path | Content |
|---|---|
| `scripts/` | One script per analysis module (see mapping table below) |
| `scripts/singlecell_knk/` | Single-cell covariation and in-silico-knockout pipeline (pathB module B2) |
| `data/` | Intermediate result tables (TSV/CSV/JSON), small enough for version control |
| `figures/` | English-labelled, paper-grade figures produced by the scripts |
| `LICENSE` | MIT |

## Script-to-section mapping

| Script | Manuscript section | Question addressed |
|---|---|---|
| `tcga_univariate_cox.R` | §3.6 | Univariate OS association of the 9-gene defense-axis score (TCGA-LUAD) |
| `wp2_subgroup_km.R` | §3.6 (new) | Score-high vs score-low KM within age/sex/stage/KEAP1-NFE2L2/EGFR subgroups |
| `wp7_5_multivariable_cox.R` | §3.6 (new) | Multivariable Cox adjusting for age, stage, TMB, KEAP1, TP53; score-by-stage interaction |
| `wp3_timeroc_curves.R` | §3.6/§2.7 | 1/3/5-year time-dependent AUC of the fixed 9-gene score (TCGA + GSE68465) |
| `wp4_drug_spearman.R` | §3.10 (new) | Spearman correlation, 9-gene score vs GDSC/PRISM IC50 in LUAD cell lines (5 PI3K/AKT/mTOR compounds) |
| `wp5_dependency_groups.R` | §3.10 (new) | Chronos dependency, high- vs low-score LUAD cell-line groups (DepMap 23Q4) |
| `wp6_cibersort.R` / `wp6_mcp_checkpoint.R` / `wp6_compare.R` | §3.6 (new) | Immune infiltration (CIBERSORT LM22 + MCP-counter), checkpoint genes, CYT proxy |
| `wp7a_mutation_landscape.R` | Discussion / supplement | Mutation frequency (KEAP1/NFE2L2/TP53/KRAS) and TMB by score group |
| `wp7b_gsea.R` | Discussion / supplement | fgsea Hallmark + targeted sets (NRF2, KEAP1-NFE2L2, ferroptosis) |
| `ftl_external_hr.R` | §3.6 | External validation of individual axis genes (FTL focus) |
| `spatial_defense_axis.py` | §3.15 | Spot-level defense/reader scoring in 4 NSCLC Visium sections (E-MTAB-13530) |
| `spatial_section_level_figs.py` | §3.15 | Section-level statistics for spatial figures (n = 4 sections) |
| `spatial_stratified_wp1.py` | §3.15 (update) | Per-section stratified re-analysis + lepidic-pattern section signature scoring |
| `singlecell_knk/*` | §3.14 | IGF2BP3-defense-axis covariation and in-silico knockout in single cells (GSE131907) |

## The 9-gene defense-axis score (pre-specified, unweighted)

`SLC7A11, SLC3A2, GCLC, GCLM, GOT2, TFRC, VDAC2, IGF2BP2, IGF2BP3`

Definition: within-cohort z-score of log2(x+1) expression per gene, arithmetic mean across the
9 genes, frozen before external evaluation. The same definition is used in every script;
several scripts contain hard assertions that re-derive reference numbers (e.g. cohort HR = 1.81)
to guard against silent口径 drift. TCGA expression: UCSC Xena HiSeqV2 (log2 scale).

## Data sources (accessions; raw data not included due to size)

| Source | Accession / release | Used in |
|---|---|---|
| TCGA-LUAD expression | UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` | all TCGA modules |
| TCGA-LUAD survival | UCSC Xena survival phenotype | WP2, WP3, WP7.5 |
| TCGA MC3 mutations / clinical | GDC API (`ssm_occurrences`, `cases`) | WP2, WP7a, WP7.5 |
| External cohort | GSE68465 | WP3 |
| GDSC drug sensitivity | Sanger GDSC release-8.4 | WP4 |
| PRISM drug sensitivity | Broad figshare 9393293 v4 | WP4 |
| DepMap CRISPR Chronos | DepMap Public 23Q4 | WP5 |
| Spatial transcriptomics | E-MTAB-13530 (De Zuani et al., Nat Commun 2024;15:4613) | spatial modules |
| Lepidic-pattern section | supplementary Visium matrix (Space Ranger 2.0.1 output; no coordinates) | WP1 |
| Histologic-pattern signatures | Xie et al., Clin Transl Med 2024;14:e1573, Table S1 | WP1 |
| Single-cell RNA-seq | GSE131907 | singlecell_knk |
| Immune signatures | LM22 (CIBERSORT), MCP-counter official signatures | WP6 |

## Environment

- R 4.6.1 (macOS): `survival`, `timeROC`, `fgsea`, `msigdbr`, `data.table`, `R.utils`,
  plus `CIBERSORT.R` v1.04 (place `CIBERSORT.R` and `LM22.txt` next to `wp6_cibersort.R`)
- Python 3.13: `scanpy 1.12.4`, `pandas`, `numpy`, `scipy`, `matplotlib`, `openpyxl`

## Reproduction notes

1. Scripts use absolute paths at the top (`DIR`, `BASE`, `PB`); edit these three variables to
   point to your local copies of the raw data before running.
2. Raw matrices (Visium `.h5`/`.mtx`, DepMap ~450 MB expression CSV, GDSC/PRISM matrices) are
   downloaded from the accessions above and are **not** vendored.
3. Negative results are part of the analysis: WP4 (drug sensitivity) and WP5 (dependency) are
   reported as-is, with Benjamini–Hochberg-adjusted P-values, in the corresponding tables under
   `data/`.

## License

MIT (see `LICENSE`). Data files remain subject to their original sources' terms.
