# Cover Letter — NAR Genomics and Bioinformatics

**Date:** June 9, 2026

**To:** The Editors, *NAR Genomics and Bioinformatics*

**Manuscript Title:** CADE: a reference-free coefficient-sensitivity ranking workflow for bulk differential expression

**Article Type:** Research Article (Method / Software)

Dear Editors,

We submit the manuscript "CADE: a reference-free coefficient-sensitivity ranking workflow for bulk differential expression" for consideration as a Research Article in *NAR Genomics and Bioinformatics*. All authors have approved the submission; the manuscript has not been published previously and is not under consideration elsewhere.

CADE addresses a specific gap in rare-disease bulk transcriptomics: matched single-cell references are unavailable, yet cell-composition shifts can dominate apparent DE. CADE reports a per-gene Composition Confounding Index (CCI) that ranks coefficients from adjustment-stable to adjustment-sensitive, with signed attenuation/amplification/reversal labels, marker-dropout bootstrap rank stability, label-permutation calibration, threshold scans, stabilised-denominator diagnostics and an isometric log-ratio (ILR) coordinate-sensitivity layer. The workflow is benchmarked on six datasets (synthetic composition-bias gradient: AUROC 0.996 vs 0.902 for unadjusted limma at 100% bias; PBMC: r = 0.993; sparse-marker scRNA-seq stress test; two independent sepsis cohorts; MAS pseudobulk validation: r = 0.779, 7/8 cell types p < 0.001) and applied to FHL peripheral blood mononuclear cells (GSE26050; n = 11 untreated FHL vs 30 healthy controls), identifying *SLC7A11* as the most composition-stable ferroptosis candidate (rank-stability probability 0.835).

All data are public (GEO: GSE26050, GSE28750, GSE66099, GSE207633). The code archive (MIT license), `renv.lock`, public GitHub repository (`https://github.com/xutaoguo55/CADE`), restricted Zenodo archive (DOI: `10.5281/zenodo.20603524`), and a `run_all.R` driver are provided. The study uses de-identified data, declares no competing interests, no specific funding, and discloses AI-assisted language editing (Claude) with author responsibility for content.

Sincerely,

Xutao Guo
Department of Hematology, Nanfang Hospital, Southern Medical University
Clinical Medical Research Center of Hematological Diseases of Guangdong Province
Email: gxt827@126.com
ORCID: https://orcid.org/0000-0001-6191-2204
