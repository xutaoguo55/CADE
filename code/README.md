# CADE: Composition-Aware Differential Expression

Analysis scripts and supplementary data tables for the manuscript:

> **CADE: reference-free coefficient-sensitivity ranking for composition-aware differential expression in bulk transcriptomics**
>
> Haiqing Zheng, Qi Wei, Hongbing Jiang, Junwei Huang, Xiaolei Wei, Yongqiang Wei, Ru Feng, Xutao Guo

## Contents

- **Core CADE method**: `cade_method.R`
- **Benchmarks**: synthetic dose-response, PBMC ground-truth, scRNA-seq pseudobulk, real-data AUROC
- **FHL application**: GEO download, DE, GSVA, WGCNA, TF target-gene scoring
- **External validation**: sepsis (GSE28750, GSE66099), MAS pseudobulk
- **Supplementary tables S1-S8**: consolidated workbook of full DE results, CCI outputs, benchmarks, validation outputs, sensitivity analyses, CADE-ILR rank-stability outputs, and multi-seed bootstrap checks; original CSV components are retained for reproducibility
- **Reproducibility**: `run_all.R` end-to-end script

## Core workflow interface

CADE is implemented as a staged workflow rather than a single black-box test. The main callable components are:

- `estimate_proportions_cade(expr, markers)`: derives reference-free marker-based relative composition weights.
- `cade_de_analysis(expr, group, weights, composition_transform="raw")`: fits paired unadjusted and composition-adjusted limma models and reports coefficient-sensitivity outputs.
- `cade_de_analysis(expr, group, weights, composition_transform="ilr")`: repeats the paired model using ILR balance covariates for closed-composition sensitivity analysis.
- `cade_bootstrap(...)`: estimates marker-dropout uncertainty and rank-stability summaries for selected genes.
- `cade_ilr_uncertainty.R` and `validate_cade_ilr_multiseed.R`: generate the ILR sensitivity, rank-stability, and multi-seed robustness tables used in the manuscript.

The standard output is a per-gene audit table with paired raw and adjusted coefficients (`logFC_unadj`, `logFC_adj`), coefficient-change fields, CCI mode, signed direction labels, selected covariates, approximate intervals, bootstrap rank summaries, and ILR raw-vs-coordinate comparisons. This is the software layer that distinguishes CADE from simply adding marker covariates to a differential-expression model.

## License

MIT License. See the package-level `LICENSE` file.

## Code fixes (v1.0.1)

- `cade_de_analysis`: added robust `group` type handling (`logical`, `character`, `factor`) with auto-detection of case/control levels
- `cade_bootstrap`: `key_genes` parameterised instead of hard-coded
- `benchmark_cade`: fixed `roc(direction)` so AUROC scores are directionally correct
- `run_cade_fhl`: removed duplicate `Delta_logFC` definition that had the opposite sign to `cade_de_analysis`
- `generate_synthetic_data`: `group` length now scales with `n_samples` instead of hard-coded 44

## Method upgrade (v1.1.0)

- `cade_de_analysis`: added `composition_transform="ilr"` for isometric log-ratio covariates from closed marker-derived weights
- `cade_de_analysis`: reports limma-derived adjusted SE, approximate CCI intervals, active composition transform, and covariates used
- `cade_bootstrap`: reports CCI rank intervals and probabilities of remaining in the top-5 lowest-CCI panel genes
- `cade_ilr_uncertainty.R`: generates CADE-ILR full DE, panel, bootstrap rank-stability, raw-vs-ILR comparison, covariate, and figure outputs
- `validate_cade_ilr_multiseed.R`: reruns CADE-ILR bootstrap under multiple random seeds and summarizes rank-stability reproducibility
- `test_cade_ilr_uncertainty.R`: quick regression checks for ILR transform dimensions, closed compositions, CCI interval bounds, and rank-stability probabilities

## Availability Note

The submission package includes the complete analysis scripts, run-order notes, dependency summary, supplementary output tables, `CITATION.cff`, Zenodo metadata, and a release-ready supplementary software archive. Cite the supplementary software archive unless a public repository or DOI is created before submission.

## Contact

Corresponding author: Xutao Guo <gxt827@126.com>
