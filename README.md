# CADE

CADE is a reference-free coefficient-sensitivity ranking workflow for bulk transcriptomic differential-expression studies where cell-composition shifts are expected but disease-matched single-cell references are unavailable.

The workflow derives marker-based relative composition weights, fits paired unadjusted and composition-adjusted limma models, and reports per-gene coefficient-sensitivity outputs: Composition Confounding Index (CCI), signed coefficient-change labels, denominator handling, bootstrap rank stability, permutation calibration, threshold scans, and ILR coordinate sensitivity.

## Repository Scope

This public GitHub repository contains the code, environment files, metadata, main result tables, and supplementary result tables needed to inspect and rerun the CADE workflow. It intentionally does not include the journal submission portal package, cover letter, manuscript DOCX/PDF files, internal review notes, or generated upload archives.

The archived review package is available as a restricted Zenodo record:

- Zenodo record: https://zenodo.org/records/20603524
- DOI: `10.5281/zenodo.20603524`
- Access right: restricted, for editorial/reviewer access

## Contents

```text
.
├── code/                  # R and Python analysis scripts
├── tables/                # Main result and benchmark CSV tables
├── supplementary/         # Supplementary workbook and raw CSV components
├── supplementary_index.md # Supplementary table index
├── Dockerfile             # Containerized reproduction entry point
├── environment.yml        # Conda environment
├── renv.lock              # R package lockfile
├── REPRODUCIBILITY.md     # Run instructions and expected outputs
├── CITATION.cff           # Citation metadata
├── .zenodo.json           # Zenodo metadata
├── LICENSE                # MIT license
└── LICENSE_CODE.md        # Code license details
```

## Quick Start

With Docker:

```bash
docker build -t cade:1.1.0 .
docker run --rm -v "$(pwd)/analysis_output:/cade/analysis_output" cade:1.1.0
```

With Conda:

```bash
conda env create -f environment.yml
conda activate cade
Rscript code/run_all.R --skip-benchmarks
```

Full benchmark reproduction:

```bash
Rscript code/run_all.R
```

Expected runtime is approximately 12 minutes for `--skip-benchmarks` and 20-21 minutes for the full workflow on the current development machine.

## Key Scripts

- `code/run_all.R`: end-to-end driver
- `code/cade_method.R`: core CADE weight estimation, paired limma models, and CCI outputs
- `code/cade_ilr_uncertainty.R`: ILR coordinate-sensitivity and rank-stability workflow
- `code/validate_cade_ilr_multiseed.R`: multi-seed bootstrap robustness checks
- `code/benchmark_method_comparison.R`: multi-method benchmark comparison
- `code/empirical_comparator_runtime_benchmark.R`: empirical comparator and runtime/scalability benchmark
- `code/generate_model_metadata_table.R`: model formula, parameter, and collinearity diagnostics

## Data Sources

The workflow uses public GEO datasets:

- GSE26050: FHL PBMC worked example
- GSE28750: sepsis whole-blood validation
- GSE66099: paediatric SIRS/sepsis validation
- GSE207633: MAS scRNA-seq pseudobulk validation

## Citation

If you use CADE, cite the manuscript and archived software package. Citation metadata are provided in `CITATION.cff`.

```text
CADE: a reference-free coefficient-sensitivity ranking workflow for bulk differential expression.
Zenodo DOI: 10.5281/zenodo.20603524
Repository: https://github.com/xutaoguo55/CADE
```

## License

The code is released under the MIT License. See `LICENSE` and `LICENSE_CODE.md`.
