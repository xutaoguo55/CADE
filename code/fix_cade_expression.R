#!/usr/bin/env Rscript
# Fix: Re-extract GSE26050 expression matrix correctly for CADE
# The processed_matrix.csv was double-log2-transformed (range [-24, 2.8]).
# Correct data from GEO has range [-4.6, 7.0] (RMA+DWD+median-centered in log2).
library(GEOquery)
library(hgu133plus2.db)
library(dplyr)

# ── Package-local path configuration ──
if (exists("CODE_DIR", envir = .GlobalEnv)) {
  SCRIPT_DIR <- normalizePath(get("CODE_DIR", envir = .GlobalEnv), mustWork = TRUE)
} else {
  SCRIPT_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
}
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)

OUT_DIR <- file.path(PROJECT_ROOT, "geo_analysis_output")
DATA_DIR <- file.path(PROJECT_ROOT, "analysis_output", "CADE")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== Fixing expression matrix for CADE ===\n")
corrected_path <- file.path(OUT_DIR, "GSE26050_expression_corrected.csv")

if (file.exists(corrected_path)) {
  cat(sprintf("Using existing corrected expression matrix: %s\n", corrected_path))
  exprs_final <- as.matrix(read.csv(corrected_path, row.names = 1, check.names = FALSE))
  cat(sprintf("  Cached matrix: %d genes x %d samples, range [%.2f, %.2f]\n",
      nrow(exprs_final), ncol(exprs_final), min(exprs_final), max(exprs_final)))
} else {
  # ── Step 1: Download GSE26050 from GEO ──
  cat("Loading GSE26050 from GEO cache or remote source...\n")
  gse <- getGEO("GSE26050", GSEMatrix=TRUE, getGPL=FALSE, destdir=OUT_DIR)
  eset <- gse[[1]]
  exprs_mat <- exprs(eset)
  cat(sprintf("  Raw matrix: %d probes x %d samples, range [%.2f, %.2f]\n",
      nrow(exprs_mat), ncol(exprs_mat), min(exprs_mat), max(exprs_mat)))

  # ── Step 2: Probe-to-gene mapping ──
  cat("Mapping probes to genes...\n")
  probe_ids <- rownames(exprs_mat)
  gene_map <- AnnotationDbi::select(hgu133plus2.db, keys=probe_ids,
                     columns=c("SYMBOL"), keytype="PROBEID")
  gene_map <- gene_map[!is.na(gene_map$SYMBOL) & gene_map$SYMBOL != "", ]

  # Keep probe with max IQR per gene
  exprs_with_genes <- data.frame(
    Probe = probe_ids,
    Gene = gene_map$SYMBOL[match(probe_ids, gene_map$PROBEID)],
    IQR = apply(exprs_mat, 1, IQR, na.rm=TRUE),
    stringsAsFactors = FALSE
  )
  exprs_with_genes <- exprs_with_genes[!is.na(exprs_with_genes$Gene), ]
  exprs_with_genes <- exprs_with_genes %>%
    group_by(Gene) %>%
    slice_max(IQR, n=1, with_ties=FALSE) %>%
    ungroup()

  cat(sprintf("  Unique genes: %d\n", nrow(exprs_with_genes)))

  # Subset to keep only mapped probes
  exprs_final <- exprs_mat[exprs_with_genes$Probe, ]
  rownames(exprs_final) <- exprs_with_genes$Gene
  cat(sprintf("  Final matrix: %d genes x %d samples, range [%.2f, %.2f]\n",
      nrow(exprs_final), ncol(exprs_final), min(exprs_final), max(exprs_final)))

  # ── Step 3: Save corrected matrix ──
  write.csv(exprs_final, corrected_path)
  cat(sprintf("  Saved corrected matrix: %s\n", corrected_path))
}

# ── Step 4: Verify against standard DE results ──
cat("\nVerifying against standard DE results...\n")
std_de <- read.csv(file.path(OUT_DIR, "GSE26050_full_DE_results.csv"))
cat(sprintf("  Standard DE file: %d genes\n", nrow(std_de)))

# Run quick limma on corrected matrix to confirm
library(limma)
group <- colnames(exprs_final) %in% c(
  "GSM639703", "GSM639704", "GSM639705", "GSM639706", "GSM639707",
  "GSM639708", "GSM639709", "GSM639710", "GSM639711", "GSM639712", "GSM639713"
)
design <- model.matrix(~ group)
fit <- lmFit(exprs_final, design)
fit <- eBayes(fit, trend=TRUE)
de_check <- topTable(fit, coef=2, number=Inf, adjust.method="BH")
de_check$Gene <- rownames(de_check)

# Compare key genes
key_genes <- c("SLC7A11", "GPX4", "FTH1", "FTL", "STAT3", "IL6", "TNF", "IL1B", "IFNG", "HMOX1", "CXCL8")
cat(sprintf("\n  %-12s %10s %10s %10s\n", "Gene", "logFC_corrected", "logFC_std_DE", "Delta"))
for(g in key_genes) {
  new_logfc <- de_check$logFC[de_check$Gene == g]
  old_logfc <- std_de$logFC[std_de$Gene == g]
  if(length(new_logfc) && length(old_logfc)) {
    cat(sprintf("  %-12s %+10.3f %+10.3f %+10.3f\n", g, new_logfc, old_logfc, new_logfc - old_logfc))
  }
}

cat("\n=== Correction complete ===\n")
cat("Next: re-run CADE using corrected expression matrix\n")
