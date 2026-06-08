# CADE Computational Environment

## R

- **Version:** >= 4.4
- **Submitted lockfile environment:** R 4.5.3 with Bioconductor 3.22, as recorded in the package-level `renv.lock`
- **Local smoke-check environment:** R 4.5.3

## Required R Packages

| Package | Version | Source | Used in |
|---------|---------|--------|---------|
| limma | 3.66.0 | Bioconductor | Core DE analysis |
| GSVA | 2.4.9 | Bioconductor | Pathway scoring |
| WGCNA | 1.74 | CRAN | Co-expression network |
| GEOquery | 2.78.0 | Bioconductor | GEO data download |
| quadprog | 1.5-8 | CRAN | QP refinement (optional) |
| pROC | 1.19.0.1 | CRAN | AUROC computation |
| celldex | >= 1.12 | Bioconductor | Monaco immune reference benchmark |
| sva | >= 3.50 | Bioconductor | SVA-like benchmark comparison |
| RUVSeq | >= 1.36 | Bioconductor | RUVg-like benchmark comparison |
| fgsea | >= 1.28 | Bioconductor | Targeted enrichment |
| msigdbr | >= 7.5 | CRAN | Gene set database |
| dplyr | >= 1.1 | CRAN | Data manipulation |
| tibble | >= 3.2 | CRAN | Data frames |
| ggplot2 | >= 3.4 | CRAN | Plotting |
| patchwork | >= 1.1 | CRAN | Multi-panel figures |
| pheatmap | >= 1.0 | CRAN | Heatmaps |
| gplots | >= 3.1 | CRAN | Plotting utilities |
| affy | >= 1.78 | Bioconductor | Microarray preprocessing |
| annotate | >= 1.78 | Bioconductor | Annotation utilities |
| hgu133plus2.db | >= 3.13 | Bioconductor | Platform annotation (GSE26050) |
| GSEABase | >= 1.62 | Bioconductor | GSEA infrastructure |
| MASS | >= 7.3 | CRAN | Statistical functions |

## Python

- **Version:** 3.9+
- **Packages:** numpy, pandas, matplotlib

## Installation

```r
# CRAN packages
install.packages(c("WGCNA", "quadprog", "pROC", "dplyr", "tibble", "ggplot2",
                   "patchwork", "pheatmap", "gplots", "msigdbr", "MASS"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "GSVA", "GEOquery", "celldex", "sva",
                       "RUVSeq", "fgsea", "affy", "annotate",
                       "hgu133plus2.db", "GSEABase"))
```

## Reproducibility Notes

- All random processes use `set.seed(42)`.
- Full pipeline runtime: ~19 minutes on a recent laptop; the `--skip-benchmarks` core-analysis rerun takes approximately 12 minutes.
- GEO downloads require an active internet connection.
- See `README_RUN_ORDER.md` for step-by-step execution order and intermediate outputs.
- Use the package-level `renv.lock` as the authoritative version record for reproducing the submitted package.
