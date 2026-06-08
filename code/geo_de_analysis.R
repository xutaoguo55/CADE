#!/usr/bin/env Rscript
# ============================================================
# GSE26050 Full Differential Expression Analysis
# Downloads real data from GEO, runs limma, performs GSEA
# ============================================================
library(GEOquery)
library(limma)
library(affy)
library(hgu133plus2.db)
library(annotate)
library(dplyr)
library(tibble)
library(ggplot2)

set.seed(42)

# ── Package-local path configuration ──
if (exists("CODE_DIR", envir = .GlobalEnv)) {
  SCRIPT_DIR <- normalizePath(get("CODE_DIR", envir = .GlobalEnv), mustWork = TRUE)
} else {
  SCRIPT_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
}
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
OUT_DIR <- file.path(PROJECT_ROOT, "geo_analysis_output")
DATA_FILES <- normalizePath(file.path(PROJECT_ROOT, "..", "..", "05_Data_Files"), mustWork = FALSE)
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# ============================================================
# Step 1: Download GSE26050
# ============================================================
cat("=== Step 1: Downloading GSE26050 from GEO ===\n")
gse <- getGEO("GSE26050", GSEMatrix=TRUE, getGPL=TRUE, destdir=OUT_DIR)
eset <- gse[[1]]

cat(sprintf("  Dimensions: %d probes x %d samples\n", nrow(eset), ncol(eset)))

# Extract phenotype data
pdata <- pData(eset)
cat(sprintf("  Phenotype columns: %s\n", paste(colnames(pdata), collapse=", ")))

# Identify FHL vs HC groups
# GSE26050 title format typically: "FHL patient N" or "Healthy control N"
sample_titles <- pdata$title
cat("\n  Sample titles:\n")
for(i in seq_along(sample_titles)) {
  cat(sprintf("    %d: %s\n", i, sample_titles[i]))
}

# Try to find group column
group_col <- NULL
for(col in colnames(pdata)) {
  vals <- as.character(pdata[[col]])
  if(any(grepl("FHL|control|healthy|normal", vals, ignore.case=TRUE))) {
    group_col <- col
    break
  }
}

if(is.null(group_col)) {
  # Use title field
  group <- ifelse(grepl("FHL|HLH|patient", sample_titles, ignore.case=TRUE), "FHL", "HC")
} else {
  vals <- as.character(pdata[[group_col]])
  group <- ifelse(grepl("FHL|HLH|patient", vals, ignore.case=TRUE), "FHL", "HC")
}

cat(sprintf("\n  Group assignment: FHL=%d, HC=%d\n", sum(group=="FHL"), sum(group=="HC")))

# ============================================================
# Step 2: Process expression data
# ============================================================
cat("\n=== Step 2: Expression data processing ===\n")

# Get expression matrix
exprs_mat <- exprs(eset)
cat(sprintf("  Expression matrix: %d probes x %d samples\n", nrow(exprs_mat), ncol(exprs_mat)))

# Check if already log2 transformed (typical for RMA/DWD microarray data on GEO).
# DWD or median-centered log2 matrices can contain negative values, so requiring
# min >= 0 incorrectly double-transforms valid centered log2 data.
expr_range <- range(exprs_mat, na.rm=TRUE)
is_log2 <- expr_range[2] <= 30 && expr_range[1] >= -30 &&
  (expr_range[2] - expr_range[1]) <= 60
cat(sprintf("  Expression range: [%.2f, %.2f]\n", expr_range[1], expr_range[2]))
cat(sprintf("  Data appears %s log2-scale\n",
            if(is_log2) "already in" else "not to be in"))

# If not log2-scale, transform raw positive intensities.
if(!is_log2) {
  exprs_mat[exprs_mat <= 0] <- min(exprs_mat[exprs_mat > 0], na.rm=TRUE) / 2
  exprs_mat <- log2(exprs_mat)
}

# ============================================================
# Step 3: Probe-to-gene mapping
# ============================================================
cat("\n=== Step 3: Probe-to-gene mapping ===\n")

probe_ids <- rownames(exprs_mat)
cat(sprintf("  Total probes: %d\n", length(probe_ids)))

# Map probes to gene symbols
gene_symbols <- getSYMBOL(probe_ids, "hgu133plus2.db")
cat(sprintf("  Probes mapped to genes: %d\n", sum(!is.na(gene_symbols))))

# For genes with multiple probes, keep the one with highest IQR
exprs_with_genes <- data.frame(
  Probe = probe_ids,
  Gene = gene_symbols,
  expr = apply(exprs_mat, 1, IQR, na.rm=TRUE),
  stringsAsFactors = FALSE
)
exprs_with_genes <- exprs_with_genes[!is.na(exprs_with_genes$Gene), ]

# Keep highest IQR probe per gene
exprs_with_genes <- exprs_with_genes %>%
  group_by(Gene) %>%
  slice_max(order_by = expr, n = 1) %>%
  ungroup()

cat(sprintf("  Unique genes after filtering: %d\n", nrow(exprs_with_genes)))

# Subset expression matrix
exprs_final <- exprs_mat[exprs_with_genes$Probe, ]
rownames(exprs_final) <- exprs_with_genes$Gene

# Save processed matrix
write.csv(exprs_final, file.path(OUT_DIR, "GSE26050_processed_matrix.csv"))
cat("  Processed matrix saved\n")

# ============================================================
# Step 4: limma differential expression
# ============================================================
cat("\n=== Step 4: limma differential expression ===\n")

# Design matrix
design <- model.matrix(~ 0 + factor(group))
colnames(design) <- c("FHL", "HC")

# Fit linear model
fit <- lmFit(exprs_final, design)

# Contrast: FHL vs HC
contrast_mat <- makeContrasts(FHL_vs_HC = FHL - HC, levels = design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2, trend=TRUE)

# Get all results
de_results <- topTable(fit2, coef=1, number=Inf, adjust.method="BH")
de_results$Gene <- rownames(de_results)

cat(sprintf("  Total genes tested: %d\n", nrow(de_results)))
cat(sprintf("  Significant (FDR<0.05): %d\n", sum(de_results$adj.P.Val < 0.05, na.rm=TRUE)))

# Create ranked list for GSEA (by t-statistic)
ranked_list <- de_results$t
names(ranked_list) <- de_results$Gene
ranked_list <- ranked_list[!is.na(ranked_list)]
ranked_list <- sort(ranked_list, decreasing=TRUE)

# Save full results
write.csv(de_results, file.path(OUT_DIR, "GSE26050_full_DE_results.csv"), row.names=FALSE)
saveRDS(ranked_list, file.path(OUT_DIR, "GSE26050_ranked_list.rds"))
cat("  Full DE results saved\n")

# ============================================================
# Step 5: Verify against known data
# ============================================================
cat("\n=== Step 5: Verification against four_disease_complete_panel ===\n")

panel_candidates <- c(
  file.path(PROJECT_ROOT, "data", "four_disease_complete_panel.csv"),
  file.path(PROJECT_ROOT, "supplementary", "four_disease_complete_panel.csv"),
  file.path(DATA_FILES, "four_disease_complete_panel.csv")
)
panel_file <- panel_candidates[file.exists(panel_candidates)][1]

if (!is.na(panel_file)) {
  panel <- read.csv(panel_file, stringsAsFactors=FALSE)

  # Extract our values for panel genes
  verification <- de_results[de_results$Gene %in% panel$Gene, ]
  verification <- verification[, c("Gene", "logFC", "adj.P.Val", "P.Value", "t", "B")]

  cat(sprintf("  Panel file: %s\n", panel_file))
  cat("\n  Comparison of key genes:\n")
  cat(sprintf("  %-12s %12s %12s %12s %s\n", "Gene", "Our_logFC", "Panel_logFC", "Our_FDR", "Match?"))
  cat("  ------------------------------------------------------------\n")
  for(g in c("SLC7A11", "SLC25A37", "GCLM", "FTH1", "GPX4", "FTL", "SLC40A1",
             "STAT3", "JAK2", "IFNG", "CXCL8", "IL1B", "TNF")) {
    our_row <- de_results[de_results$Gene == g, ]
    panel_row <- panel[panel$Gene == g, ]
    if(nrow(our_row) > 0 && nrow(panel_row) > 0) {
      match_str <- if(abs(our_row$logFC - panel_row$FHL_logFC) < 0.1) "MATCH" else "DIFF"
      cat(sprintf("  %-12s %+12.4f %+12.4f %12.2e %s\n",
                  g, our_row$logFC, panel_row$FHL_logFC, our_row$adj.P.Val, match_str))
    } else if(nrow(our_row) > 0) {
      cat(sprintf("  %-12s %+12.4f %12s %12.2e (not in panel)\n",
                  g, our_row$logFC, "N/A", our_row$adj.P.Val))
    } else {
      cat(sprintf("  %-12s %12s %+12.4f %12s (not in our results)\n",
                  g, "N/A", panel_row$FHL_logFC, "N/A"))
    }
  }

  write.csv(verification, file.path(OUT_DIR, "verification_against_panel.csv"), row.names=FALSE)
} else {
  cat("  Optional panel file not found; skipping legacy panel verification.\n")
}

cat("\n=== Analysis complete ===\n")
cat(sprintf("Output directory: %s\n", OUT_DIR))
