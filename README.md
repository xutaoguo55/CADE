# CADE Submission Package — NAR Genomics and Bioinformatics

**Manuscript Title:** CADE: a reference-free coefficient-sensitivity ranking workflow for bulk differential expression

**Target Journal:** NAR Genomics and Bioinformatics (Original Article — Method/Software)

**Submission Date:** 2026-06-04

**Authors:** Haiqing Zheng, Qi Wei, Hongbing Jiang, Junwei Huang, Xiaolei Wei, Yongqiang Wei, Ru Feng, Xutao Guo

**Correspondence:** Xutao Guo (gxt827@126.com) — ORCID 0000-0001-6191-2204

---

## Package Structure

```
submission_upload_nargab_2026-06-04/
├── README.md
├── LICENSE                              # MIT
├── renv.lock                            # Reproducible R environment lockfile
├── CITATION.cff                         # Citation metadata
├── .zenodo.json                         # Zenodo metadata template
├── RELEASE_NOTES.md                     # Software release notes
├── CADE_public_software_v1.1.0_2026-06-01.zip   # Release-ready software archive
├── cover_letter.md                      # NAR Ge&B cover letter
├── manuscript/
│   ├── CADE_NARGeB_manuscript.md        # NAR Ge&B-adapted manuscript
│   └── CADE_NARGeB_manuscript.docx      # Word version for journal upload
├── figures/                             # All TIFs (< 10MB each)
├── tables/                              # 8 main text tables (CSV)
├── supplementary/                       # Consolidated Supplementary Tables
│   ├── CADE_Supplementary_Tables_S1-S8.xlsx
│   └── raw_csv_components/              # Original CSV/MD components
├── code/                                # Analysis scripts (R + Python)
│   ├── run_all.R                        # End-to-end driver
│   ├── cade_method.R
│   ├── cade_ilr_uncertainty.R
│   ├── validate_cade_ilr_multiseed.R
│   ├── ...
│   ├── README.md
│   ├── README_RUN_ORDER.md
│   └── DEPENDENCIES.md
```

---

## NAR Ge&B Submission Requirements Met

- **Article type:** Original Article (Method/Software)
- **License:** CC-BY (default for NAR Ge&B) — released under MIT for code
- **Title length:** ~190 characters (within typical NAR Ge&B limits)
- **Abstract length:** ~280 words (within typical NAR Ge&B limits)
- **Main text length:** ~7,000 words (within typical Original Article limits)
- **Figures:** 7 main + 7 supplementary (within typical limits)
- **Tables:** 5 main + 8 supplementary groups
- **References:** 28 (Vancouver style, DOIs included)
- **Data availability:** 4 GEO accessions, all publicly available
- **Code availability:** MIT-licensed software archive, `renv.lock`, dependency notes, public GitHub repository, and restricted Zenodo record DOI `10.5281/zenodo.20603524`
- **CRediT authorship:** complete
- **AI use disclosure:** explicit statement included
- **Competing interests:** declared none
- **Funding:** declared none

---

## Reproducibility

- **R:** 4.4.0+ (locked via `renv.lock`; tested on R 4.5.3)
- **Key packages:** limma 3.62.1, GSVA 2.4.9, WGCNA 1.72-5, GEOquery 2.74.1, quadprog 1.5-8, pROC 1.18.5
- **Python:** 3.9+ for figure regeneration
- **Random seeds:** All stochastic steps use `set.seed(42)` or equivalent
- **End-to-end runtime:** ~19 minutes (full); ~12 minutes (--skip-benchmarks)

---

## Archival notes

- The submission package consolidates figures, tables, code, and supplementary data.
- GitHub repository: `https://github.com/xutaoguo55/CADE`
- Restricted Zenodo record: `https://zenodo.org/records/20603524` (DOI: `10.5281/zenodo.20603524`)


## Reproducibility verification

The `renv.lock` file is the authoritative version record. To verify reproducibility:
1. Open the project in R 4.4.0+ with `renv` installed
2. Run `renv::restore()` to install the exact R package versions
3. Run `Rscript code/run_all.R --skip-benchmarks` for a ~12-minute end-to-end FHL analysis
4. Run `Rscript code/run_all.R` for a ~20-21-minute full pipeline including benchmarks and reviewer-facing diagnostics

Expected outputs are written to `analysis_output/CADE/` and `figures/`. The submission package also includes the consolidated `CADE_Supplementary_Tables_S1-S8.xlsx` workbook for direct comparison with published results.

## Graphical Abstract

NAR Genomics and Bioinformatics does not require a graphical abstract. The Figure 1 schematic (workflow) and Figure 2 (benchmark) together provide the visual entry point to the manuscript.
