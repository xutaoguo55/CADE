
suppressPackageStartupMessages({library(limma); library(quadprog)})
PR <- "/Users/guoxutao/.openclaw/workspace/HLH_Research/03_Papers_Manuscripts/CADE_Submission_Package"
cm <- readLines(file.path(PR, "submission_upload_nargab_2026-06-04/code/cade_method.R"))
fn <- cm[1:728]
cp <- tempfile(fileext=".R")
writeLines(fn, cp)
source(cp, echo=FALSE); unlink(cp)

exprs_mat <- as.matrix(read.csv(file.path(PR, "geo_analysis_output/GSE26050_expression_corrected.csv"), row.names=1, check.names=FALSE))
true_group <- c(rep(1, 11), rep(0, 33))
marker_list <- list(
  CD8_Tcells=c("CD8A","CD8B","CD3D","CD3E","TRAC","CD2","GZMK","CCL5","PRF1"),
  CD4_Tcells=c("CD4","IL7R","CCR7","LEF1","MAL","TCF7","LDHB"),
  NK_cells=c("NKG7","GNLY","KLRD1","KLRB1","GZMB","CTSW"),
  B_cells=c("CD19","MS4A1","CD79A","CD79B","BANK1","CD22","PAX5"),
  Monocytes=c("LYZ","CD14","FCGR3A","MS4A7","ITGAM","CCR2","CD163","CSF1R","S100A8"),
  Macrophages=c("CD68","CD163","MRC1","MSR1","MARCO","CSF1R"),
  Neutrophils=c("FCGR3B","CSF3R","S100A8","S100A9","CXCR2","ITGAM","MMP9"),
  Erythrocytes=c("HBB","HBA1","HBA2","HBD","AHSP","ALAS2","SLC25A37")
)
marker_list <- lapply(marker_list, function(g) intersect(g, rownames(exprs_mat)))
marker_list <- marker_list[sapply(marker_list, length) >= 4]
ferroptosis_genes <- c("SLC7A11","IFNG","TFRC","FTH1","TNF","IL1B","NFKB1","FTL",
                        "STAT3","JAK2","NFE2L2","SLC25A37","CXCL8","SLC40A1","IL6",
                        "GPX4","GCLM","HMOX1","STAT1","NCOA4")
ferroptosis_genes <- intersect(ferroptosis_genes, rownames(exprs_mat))

# Test n_bootstrap values
n_values <- c(50, 100, 200, 400, 800)
results <- data.frame(n_bootstrap=integer(), Gene=character(), Mean_CCI=numeric(), SD_CCI=numeric(), Top5_Prob=numeric(), stringsAsFactors=FALSE)

for(n_boot in n_values) {
  cat(sprintf("Running n_bootstrap=%d ...\n", n_boot))
  set.seed(42)
  boot <- cade_bootstrap(exprs_mat, marker_list, true_group, n_bootstrap=n_boot, top_cts=4, cci_variant="stabilized", verbose=FALSE)
  rownames(boot$summary) <- boot$summary$Gene
  for(g in ferroptosis_genes) {
    if(g %in% rownames(boot$summary)) {
      boot_summary <- boot$summary[boot$summary$Gene == g, ]
      results <- rbind(results, data.frame(
        n_bootstrap=n_boot,
        Gene=g,
        Mean_CCI=boot_summary$Boot_CCI_median,
        SD_CCI=(boot_summary$Boot_CCI_upper - boot_summary$Boot_CCI_lower) / 3.29,  # Approx SD from 95% CI
        Top5_Prob=boot_summary$Prob_top5_low_CCI
      ))
    }
  }
}

# Save
write.csv(results, file.path(PR, "submission_upload_nargab_2026-06-04/analysis_output/CADE/Table_S46_BootstrapConvergence.csv"), row.names=FALSE)
cat("Saved bootstrap convergence results\n")

# Summary: how much do the CCI values change with n?
suppressPackageStartupMessages(library(dplyr))
conv <- results %>% group_by(Gene) %>% summarise(
  Range_Mean_CCI = max(Mean_CCI) - min(Mean_CCI),
  Range_Top5 = max(Top5_Prob) - min(Top5_Prob),
  CV_Mean_CCI = sd(Mean_CCI) / mean(Mean_CCI) * 100
)
cat("\n=== Bootstrap convergence summary ===\n")
print(head(conv, 20))
cat(sprintf("\nMean range of Mean_CCI across n_boot values: %.3f\n", mean(conv$Range_Mean_CCI)))
cat(sprintf("Mean range of Top5_Prob across n_boot values: %.3f\n", mean(conv$Range_Top5)))
cat(sprintf("Mean CV of Mean_CCI: %.1f%%\n", mean(conv$CV_Mean_CCI)))
cat("If ranges are small, bootstrap has converged; if large, more iterations needed.\n")
