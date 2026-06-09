# CADE Release Notes

## v1.1.0 - 2026-06-01

This release accompanies the CADE manuscript as a reference-free coefficient-sensitivity ranking workflow for bulk transcriptomic studies.

### Main Features

- Marker-derived relative composition weight estimation from curated marker sets.
- Paired unadjusted and composition-adjusted limma models.
- Composition Confounding Index (CCI) outputs with legacy and stabilized denominator handling.
- Signed coefficient-change labels for attenuation, amplification, reversal, and stable effects.
- Marker-dropout bootstrap rank-stability summaries.
- Label-permutation calibration and threshold sensitivity analyses.
- ILR coordinate-sensitivity analysis for closed composition weights.
- End-to-end reproducibility driver through `code/run_all.R`.

### Reproducibility Materials

- `renv.lock` records the R package environment.
- `code/README_RUN_ORDER.md` documents the run order and expected runtime.
- `code/test_cade_ilr_uncertainty.R` provides quick regression checks for ILR transform dimensions, CCI interval bounds, and bootstrap probability outputs.
- Supplementary tables contain full CADE, CADE-ILR, benchmark, and validation outputs.

### Scope Boundaries

- CADE estimates marker-derived relative weights, not externally calibrated absolute cell-type proportions.
- CCI is a coefficient-sensitivity ranking metric rather than a causal mediation fraction.
- CADE evaluates sensitivity of bulk DE coefficients to marker-derived composition adjustment rather than cell-type-specific DE.
- Matched nnls benchmark results define an ideal matched-reference comparator and are reported separately from CIBERSORTx-style workflows.

### Implementation Notes

- The current implementation is script-based rather than a formal R package.
- Public GitHub repository: `https://github.com/xutaoguo55/CADE`.
- Restricted Zenodo record: `https://zenodo.org/records/20603524` (DOI: `10.5281/zenodo.20603524`).
- The FHL biological application is used as a reproducible disease-focused demonstration and is interpreted through the CADE sensitivity framework.
