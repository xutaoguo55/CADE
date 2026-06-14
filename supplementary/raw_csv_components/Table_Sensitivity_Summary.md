# Table S7. Parameter sensitivity summary

# This supplementary summary links the three parameter-sensitivity grids reported in
# Supplementary Table S7 and the parameter-sensitivity analyses shown in Figure S8.

| Dimension | Source table | Main diagnostic | Summary |
|---|---|---|---|
| Marker number and marker noise | Table_Sensitivity_MarkerQuality.csv | AUROC, mean proportion correlation, CCI separation, runtime | CADE retained high DE-recovery AUROC across marker-number and noise settings, while estimated-weight correlation was more sensitive to marker quality. |
| QP convergence tolerance and maximum iterations | Table_Sensitivity_QPParams.csv | Convergence rate, AUROC, runtime | The linear softmax initialization captured most of the DE-relevant composition signal; QP refinement improved convergence behavior and provided incremental gains when additional computation was acceptable. |
| Number of composition covariates | Table_Sensitivity_TopCTs.csv | AUROC and cumulative variance explained | Including the top three highest-variance composition weights captured most cross-sample composition variance; adding more covariates produced little additional DE-recovery gain. |

These summaries support the default analysis choices used in the manuscript: curated marker sets with at least 5-7 markers per cell type where possible, CADE-lite for routine screening, optional QP refinement when weight ranking is prioritized, and a small number of high-variance composition covariates in the adjusted DE model.
