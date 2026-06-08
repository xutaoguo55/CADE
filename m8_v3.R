
suppressPackageStartupMessages({library(limma); library(quadprog)})
# Source from the same dir the script lives in
src <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "cade_method.R"))
cm <- readLines(src)
fn <- cm[1:728]
cp <- tempfile(fileext=".R")
writeLines(fn, cp)
source(cp, echo=FALSE)
unlink(cp)

# Data
data_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "..", "..", ".."))
exprs_mat <- as.matrix(read.csv(file.path(data_root, "geo_analysis_output", "GSE26050_expression_corrected.csv"), row.names=1, check.names=FALSE))
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

all_results <- list()
for(seed in c(42, 123, 456, 789, 2026, 31415, 271828, 141421, 618034, 999999)) {
  set.seed(seed)
  candidate_genes <- setdiff(rownames(exprs_mat), ferroptosis_genes)
  null_genes <- sample(candidate_genes, 100)
  prop <- estimate_proportions_cade(exprs_mat, marker_list, max_iter=30, tol=1e-5, verbose=FALSE)
  de <- cade_de_analysis(exprs_mat, true_group, prop$proportions, top_cts=4, cci_variant="stabilized", verbose=FALSE)
  rownames(de) <- de$Gene
  cci_col <- "CCI_stabilized"
  null_cci <- sapply(null_genes, function(g) {
    if(g %in% rownames(de)) de[g, cci_col] else NA
  })
  all_results[[as.character(seed)]] <- data.frame(
    Seed = seed,
    Gene = null_genes,
    CCI = as.numeric(null_cci),
    stringsAsFactors = FALSE
  )
}

combined <- do.call(rbind, all_results)
# Write to current dir
out_file <- "Table_S_NullControl_10Seeds.csv"
write.csv(combined, out_file, row.names=FALSE)
cat("Saved:", out_file, "\n")

suppressPackageStartupMessages(library(dplyr))
summary_stats <- combined %>% group_by(Seed) %>% summarise(
  N_genes = n(),
  Mean_CCI = mean(CCI, na.rm=TRUE),
  Median_CCI = median(CCI, na.rm=TRUE),
  SD_CCI = sd(CCI, na.rm=TRUE),
  Pct_CCI_gt_0_5 = mean(CCI > 0.5, na.rm=TRUE) * 100,
  Pct_CCI_lt_0_2 = mean(CCI < 0.2, na.rm=TRUE) * 100
)
print(summary_stats)
write.csv(summary_stats, "Table_S_NullControl_10Seeds_Summary.csv", row.names=FALSE)
