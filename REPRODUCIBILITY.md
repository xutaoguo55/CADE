# CADE Reproducibility Guide

**Version:** 1.1.0  
**Last updated:** 2026-06-14

---

## Quick Start (3 options)

### Option 1: Docker (recommended)

```bash
docker build -t cade:1.1.0 .
docker run --rm -v $(pwd)/analysis_output:/cade/analysis_output cade:1.1.0
```

To run with benchmarks:

```bash
docker run --rm -v $(pwd)/analysis_output:/cade/analysis_output cade:1.1.0 \
  Rscript code/run_all.R
```

Expected runtime: ~20-21 min (full), ~12 min (--skip-benchmarks).

### Option 2: Conda

```bash
conda env create -f environment.yml
conda activate cade
Rscript code/run_all.R --skip-benchmarks
```

### Option 3: Manual + renv

```r
# In R:
install.packages("renv", repos = "https://cloud.r-project.org")
renv::restore(lockfile = "renv.lock")
source("code/run_all.R")
```

---

## Environment Verification

```r
# Check that all required packages load
required <- c("limma", "GSVA", "WGCNA", "GEOquery", "sva", "RUVSeq",
              "pROC", "quadprog", "EPIC", "fgsea", "dplyr", "ggplot2")
for (p in required) library(p, character.only = TRUE)
cat("All packages OK\n")
sessionInfo()
```

## Expected Outputs

After a successful run, the following directories are populated locally:

| Path | Contents |
|------|----------|
| `analysis_output/CADE/` | CADE DE results, CCI tables, ILR outputs |
| `figures/` | Regenerated figures (TIF format; generated locally, not tracked in this public repository) |

## Key Scripts

| Script | Purpose | Runtime |
|--------|---------|---------|
| `code/run_all.R` | End-to-end driver | 12-21 min |
| `code/cade_method.R` | Core CADE/CCI computation | < 1 min |
| `code/cade_ilr_uncertainty.R` | ILR sensitivity | 3-5 min |
| `code/geo_de_analysis.R` | GEO download + standard DE | 3-5 min |
| `code/cade_pbmc_benchmark.R` | PBMC benchmark | 2-3 min |
| `code/benchmark_method_comparison.R` | Multi-method comparison | 2-3 min |
| `code/empirical_comparator_runtime_benchmark.R` | Empirical comparator + runtime/scalability benchmark | 1-2 min |
| `code/generate_model_metadata_table.R` | Model metadata, parameter and collinearity diagnostics | < 1 sec |

## Data Sources

All input data are publicly available from GEO:
- GSE26050: FHL PBMC (n=11 FHL vs 30 controls)
- GSE28750: Sepsis whole blood
- GSE66099: Paediatric SIRS/sepsis
- GSE207633: MAS scRNA-seq

## Random Seeds

All stochastic steps use `set.seed(42)` or equivalent.

## License

MIT License. See `LICENSE` and `LICENSE_CODE.md`.
