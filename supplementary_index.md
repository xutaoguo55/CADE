

## Table of Contents

1. **Main Figures (1-7)**: workflow, benchmark validation, FHL CCI ranking, cell-composition profile, pathway activity, WGCNA, target-gene signatures
2. **Supplementary Figures (S1-S7)**: cross-disease heatmap, sepsis validation, MAS pseudobulk validation, sensitivity analyses, extended benchmark, parameter sensitivity, CADE-ILR rank stability
3. **Main Tables (1A-5)**: dose-response benchmark, scRNA-seq pseudobulk benchmark, MAS pseudobulk validation, modular contribution, FHL ferroptosis panel, immune cell-type scores, ferroptosis gene expression, targeted enrichment
4. **Supplementary Tables (S1-S8 groups)**: GSE26050 DE and reproducibility checks; marker sets and pathway scores; WGCNA and cross-disease context; primary FHL CADE outputs and CCI calibration; benchmark outputs; external validation and MAS pseudobulk validation outputs; parameter sensitivity; CADE-ILR robustness
5. **Software**: `CADE_public_software_v1.1.0_2026-06-01.zip` (MIT license)

---
# CADE Submission Package — Supplementary Index

**Manuscript:** CADE: a reference-free coefficient-sensitivity ranking workflow for composition-aware differential expression analysis in small-sample bulk transcriptomics

**Target Journal:** NAR Genomics and Bioinformatics

**Last update:** June 4, 2026

---

## Main Figures (single multi-panel TIF per figure)

| # | File | Description |
|---|------|-------------|
| Figure 1 | `Figure1_Workflow_Summary.tif` | (A) Cell-composition confounding; (B) CADE three-stage framework; (C) CCI interpretation; (D) Outputs |
| Figure 2 | `Figure2_CADE_Benchmark.tif` | (A) Dose-response synthetic benchmark; (B) Sparse-marker scRNA-seq pseudo-bulk; (C) CCI calibration; (D) Per-cell-type weight accuracy (PBMC + MAS); (E) Five-method AUROC on GSE26050 probe-level data |
| Figure 2E | `Figure2E_5method_AUROC.tif` | Standalone high-resolution version of Figure 2E for upload systems that request single-panel image files |
| Figure 3 | `Figure3_CADE_FHL_Multi_Panel.tif` | (A) CCI waterfall (19 genes); (B) Adjusted vs raw logFC scatter; (C) Bootstrap & approx intervals; (D) Permutation null distributions |
| Figure 4 | `Figure4_Composition_Correlation.tif` | (A) Cell-type marker score Δ; (B) Direction-only summary; (C) Ferroptosis-gene × cell-type correlation heatmap |
| Figure 5 | `Figure5_Pathway_Analysis.tif` | (A) Ferroptosis/iron GSVA pathways; (B) Multi-cell-death comparison; (C) Pathway-level Δ vs significance |
| Figure 6 | `Figure6_WGCNA_Modules.tif` | (A) All 9 module–FHL correlations; (B) Module size × |r| × FDR bubble plot; (C) Module 2 hub network; (D) Module 6 hub network |
| Figure 7 | `Figure7_TF_TargetGene_Signature.tif` | (A) Per-target Δ; (B) Aggregate target Δ; (C) Target-gene scoring scope |

## Supplementary Figures

| # | File | Description |
|---|------|-------------|
| Figure S1 | `SuppFigure_S1_CrossDisease_Heatmap.tif` | Cross-disease ferroptosis expression heatmap (platform-qualified) |
| Figure S2 | `SuppFigure_S2_Sepsis_Validation.tif` | (A) Genome-wide CCI distribution; (B) Adjusted vs raw logFC; (C) Sepsis panel waterfall; (D) Cohort and CCI summary panel |
| Figure S3 | `SuppFigure_S3_MAS_Validation.tif` | 2 × 4 per-cell-type scatter — CADE weights vs scRNA-seq ground truth, with annotated r/RMSE/MAE/p |
| Figure S4 | `SuppFigure_S4_Sensitivity_Analyses.tif` | (A) CCI threshold sensitivity; (B) Parametric bootstrap CIs for 19 panel genes; (C) Null-control vs panel CCI distribution; (D) WGCNA module × CCI tier composition; (E) PCA on CADE weights; (F) CCI vs cell-type correlation strength |
| Figure S5 | `SuppFigure_S5_Extended_Benchmark.tif` | (A) Synthetic 5-method benchmark; (B) Panel CCI by cell-type tracking category; (C) Genome-wide CCI vs \|logFC\| landscape; (D) Per-module CCI Mann-Whitney with BH-FDR; (E) CCI density across composition-sensitivity module classes; (F) Synthesis text panel |
| Figure S6 | `SuppFigure_S6_Parameter_Sensitivity.tif` | Parameter sensitivity analyses: (A) marker number and noise fraction effects on AUROC, proportion correlation, and CCI separation; (B) QP tolerance and max-iteration effects on convergence, AUROC, and runtime; (C) number of cell-type covariates (`top_cts`) vs AUROC and cumulative variance explained |
| Figure S7 | `SuppFigure_S7_CADE_ILR_RankStability.tif` | CADE-ILR marker-dropout bootstrap CCI and probability of top-5 low-CCI rank stability for the 20-gene FHL panel |

## Main Text Tables

| Table | File | Section | Description |
|-------|------|---------|-------------|
| Table 1A | `Table_01A_Benchmark_Gradient.csv` | Section 3.1 | Dose-response benchmark performance |
| Table 1B | `Table_01B_scRNA_Pseudobulk_Benchmark.csv` | Section 3.1 | Sparse-marker scRNA-seq benchmark |
| Table 1C | `Table_01C_MAS_Pseudobulk_Validation.csv` | Section 3.1 | MAS pseudobulk validation |
| Table 1D | `Table_01D_CADE_Modular_Contribution.csv` | Section 3.1 | Modular contribution of CADE beyond a conventional marker-covariate adjusted DE table |
| Table 2  | `Table_02_CADE_Ferroptosis_Panel.csv` | Section 3.2 | CADE composition-adjusted DE for 19-gene ferroptosis panel |
| Table 3  | `Table_03_Immune_CellType_Scores.csv` | Section 3.3 | Immune cell-type signature scores |
| Table 4  | `Table_04_Ferroptosis_Gene_Expression.csv` | Section 3.4 | Ferroptosis gene expression |
| Table 5  | `Table_05_Targeted_Enrichment.csv` | Section 3.5 | Targeted over-representation results |

## Supplementary Tables (S1-S8; provided in `CADE_Supplementary_Tables_S1-S8.xlsx`)

The formal supplementary tables have been consolidated into eight themed worksheet groups to keep the initial review package focused. Each group contains the indicated source-table worksheets; the original CSV components are retained under `supplementary/raw_csv_components/` for auditability and exact reuse.

| # | Worksheet group | Included source tables | Description |
|---|-----------------|------------------------|-------------|
| S1 | GSE26050 differential expression and reproducibility checks | S01, S02, S08 | Full GSE26050 DE results, cross-subtype consistency summary, and file-level reproducibility cross-check |
| S2 | Marker sets, pathway scores, cell-type summaries, and TF target signatures | S03-S07, S12 | Immune marker scores, ferroptosis-immune correlations, targeted enrichment/CAMERA context, GSVA and cell-death pathway scores, and TF target-gene scoring |
| S3 | WGCNA and cross-disease context | S09-S11, S13 | WGCNA module associations, gene-module assignments, hub genes, and platform-qualified cross-disease panel |
| S4 | Primary FHL CADE outputs and CCI calibration | S14, S15, S18, S24 | Full CADE DE results, ferroptosis/iron/immune panel, permutation calibration, and FHL marker-derived weights |
| S5 | Benchmark outputs | S16, S17, S19, S20 | Synthetic benchmark summaries and replicate-level outputs, PBMC benchmark context, and sparse-marker scRNA-seq pseudobulk benchmark |
| S6 | External validation and MAS pseudobulk validation | S21-S23, S25-S29 | Sepsis validation, GSE66099 cross-context CCI, MAS pseudobulk full results, weight validation, summaries, and estimated weights |
| S7 | Parameter sensitivity analyses | `S7_Sensitivity_MarkerQuality`, `S7_Sensitivity_QPParams`, `S7_Sensitivity_TopCTs`, `S7_Sensitivity_Summary` | Marker quality/noise, QP convergence, number of covariates, and narrative sensitivity summary |
| S8 | CADE-ILR robustness and multi-seed stability | S34-S40 | CADE-ILR full results, focused panel, bootstrap rank stability, raw-vs-ILR comparison, ILR covariates, and multi-seed stability outputs |
| S8 (extension) | Full-permutation sensitivity for CCI null | Table_S_FullPermutation_vs_Shared, Table_S_FullPermutation_RawPermCCI | Comparison of shared-weight vs full-permutation null distributions; supports the qualitative CCI ranking under both designs |
| S9 (extension) | ILR basis-stability analysis | Table_S42_ILR_BasisStability | CADE-ILR CCI ranking across 4 alternative Helmert bases; shows moderate-to-high rank stability for the most common contrast designs |
| S10 (extension) | GSE26050 pathway-enrichment comparison | Table_S43_GSE26050_PathwayEnrichment, Table_S44_CADE_Top20_vs_GSE26050 | Top-100 CADE DE genes are enriched for inflammatory (19x) and inflammasome (30x) pathways but NOT for cytotoxic, IFN, or JAK-STAT pathways; supports FHL as secondary inflammation |
| S11 (extension) | Null-control seed sensitivity | Table_S45_NullControl_10Seeds, Table_S45_NullControl_10Seeds_Summary | 10-seed null-control CCI distribution is stable (mean 0.405-0.530, SD of seed means 0.034) |
| S12 (extension) | Treatment sensitivity analysis | Table_S12_C2_TreatmentConfound_Summary, Table_S12_C2_TreatmentConfound_3of11, Table_S12_C2_TreatmentConfound_6of11, Table_S12_C2_TreatmentConfound_11of11 | 100-replicate sensitivity simulation of hypothetical dexamethasone/etoposide exposure at 25%, 50%, and 100% treatment coverage; quantifies potential logFC shifts for 32 drug-responsive genes if treatment-naïve status were misclassified (Supplementary Table S12) |

## Core Analysis Scripts

| Script | Description |
|--------|-------------|
| `run_all.R` | End-to-end reproducibility script (`--step N`, `--skip-benchmarks`) |
| `cade_method.R` | CADE core: weight estimation, adjusted DE, CCI |
| `cade_ilr_uncertainty.R` | CADE-ILR compositional covariates, approximate CCI intervals, and bootstrap rank-stability outputs |
| `validate_cade_ilr_multiseed.R` | Multi-seed CADE-ILR bootstrap reproducibility validation |
| `test_cade_ilr_uncertainty.R` | Quick regression checks for ILR transform and uncertainty/rank-stability columns |
| `geo_de_analysis.R` | GEO data download (GSE26050) + standard DE |
| `fix_cade_expression.R` | Expression matrix correction for CADE |
| `cci_permutation_test.R` | Permutation test for CCI calibration |
| `deeper_analysis_v3.R` | GSVA, immune scoring, targeted enrichment |
| `wgcna_standalone.R` | WGCNA co-expression network analysis |
| `additional_analysis.R` | TF target-gene scoring, cross-disease comparison |
| `cade_real_scRNA_benchmark_v2.R` | Real scRNA-seq pseudo-bulk benchmark |
| `cade_pbmc_benchmark.R` | PBMC ground-truth benchmark |
| `cade_external_validation_v4.R` | Sepsis external validation (GSE28750) |
| `generate_figure2_method_comparison.py` | Figure 2E five-method AUROC bar chart |
| `generate_supp_figure_S2_sepsis.py` | Sepsis Supplementary Figure S2 generation |

---

*All scripts use `PROJECT_ROOT`-relative paths. Cross-disease comparison (S13) is platform-qualified (microarray vs 10x).*
