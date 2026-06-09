# CADE Run Order

Run commands from the package root:

```bash
Rscript code/run_all.R
```

For a shorter core-analysis rerun:

```bash
Rscript code/run_all.R --skip-benchmarks
```

Use `--skip-benchmarks` for routine reproducibility checks of the FHL worked example, CADE-ILR rank-stability layer, CCI calibration, WGCNA integration, and manuscript-output generation. Use the full command when the benchmark tables and validation figures need to be regenerated.

## Expected Inputs

- Public GEO datasets listed in the manuscript and README.
- Included processed supplementary tables for manuscript consistency checks.
- R package versions recorded in `renv.lock` and summarized in `code/DEPENDENCIES.md`.

## Interpretation Boundaries

- CADE estimates marker-derived relative weights for within-study adjustment; it does not estimate externally calibrated absolute cell-type proportions.
- CCI is a coefficient-sensitivity ranking metric for bulk DE coefficients; causal mediation and cell-type-specific differential-expression tests are separate analysis targets.
- Synthetic and pseudobulk modules provide ground-truth performance checks. Real-data rank-recovery modules are internal consistency checks because their labels are defined from CADE-derived CCI rankings.
- The matched-reference nnls benchmark defines an ideal matched-reference comparator; CIBERSORTx-style workflows with independent references are separate literature comparators.

## Pipeline Steps

1. Download or load GEO expression matrices.
2. Reconstruct the FHL expression matrix and standard differential-expression analysis.
3. Estimate marker-derived CADE composition weights, run composition-adjusted limma models, and calculate raw-weight CCI.
4. Run CADE-ILR compositional covariates, approximate CCI intervals, and marker-dropout rank stability.
5. Run permutation calibration for key genes.
6. Generate WGCNA co-expression outputs before downstream integration.
7. Generate GSVA, immune scoring, targeted enrichment, and WGCNA-integrated summaries.
8. Generate additional pathway and cross-context analyses.
9. Generate real scRNA-seq benchmark outputs unless `--skip-benchmarks` is used.
10. Generate PBMC benchmark outputs unless `--skip-benchmarks` is used.
11. Generate sepsis external-validation outputs unless `--skip-benchmarks` is used.
12. Generate empirical comparator and runtime/scalability benchmark outputs unless `--skip-benchmarks` is used.

## Runtime

Expected runtime is approximately 20-21 minutes on a recent laptop for the full workflow. The `--skip-benchmarks` mode runs the core analysis in approximately 12 minutes; the CADE-ILR step uses 200 marker-dropout bootstrap iterations and takes about 2 minutes on the current machine. The empirical comparator/runtime benchmark adds approximately 1-2 minutes on the current machine.

All stochastic steps use fixed seeds where applicable. The CADE-ILR step accepts `--seed`, `--n-bootstrap`, `--top-cts`, and `--out-dir` arguments for reproducibility checks without overwriting the main outputs. For example:

```bash
Rscript code/cade_ilr_uncertainty.R --n-bootstrap 200 --seed 42
Rscript code/cade_ilr_uncertainty.R --n-bootstrap 200 --seed 2026 --out-dir analysis_output/CADE/seed_2026
```

To check whether the ILR bootstrap ranking is stable across random marker-dropout streams, run:

```bash
Rscript code/validate_cade_ilr_multiseed.R --n-bootstrap 200 --seeds 42,2026,31415
```

Outputs are written under `analysis_output/CADE/` and can be copied to `figures/`, `tables/`, and `supplementary/` for submission packaging.
