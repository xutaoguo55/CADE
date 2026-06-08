#!/usr/bin/env Rscript
# CCI Permutation Test v2: two-sided statistical significance framework
# Tests whether CCI values deviate from null in EITHER direction:
#   P_lower: P(CCI_perm <= CCI_obs) -> observed CCI is lower than the permuted-label null
#   P_upper: P(CCI_perm >= CCI_obs) -> observed CCI is higher than the permuted-label null
library(limma)
library(quadprog)

# ── Path configuration (adjust PROJECT_ROOT for your setup) ──
get_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args_full, value = TRUE)
  if(length(hit) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", hit[1]))))
  }
  try_frame <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if(!is.null(try_frame) && nzchar(try_frame)) {
    return(dirname(normalizePath(try_frame)))
  }
  getwd()
}
resolve_project_root <- function(script_dir) {
  is_pkg_root <- function(p) {
    dir.exists(p) &&
      file.exists(file.path(p, "README.md")) &&
      dir.exists(file.path(p, "code"))
  }
  candidates <- unique(c(
    normalizePath(file.path(script_dir, ".."), mustWork = FALSE),
    normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
  ))
  for(p in candidates) {
    if(is_pkg_root(p)) {
      return(normalizePath(p, mustWork = TRUE))
    }
  }
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}
SCRIPT_DIR <- get_script_dir()
PROJECT_ROOT <- resolve_project_root(SCRIPT_DIR)

# ── Runtime options ──
# Usage:
#   Rscript cci_permutation_test.R
#   Rscript cci_permutation_test.R --mode legacy
#   Rscript cci_permutation_test.R --mode stabilized
#   Rscript cci_permutation_test.R --mode stabilized --n-perm 1000
args <- commandArgs(trailingOnly = TRUE)
mode <- "stabilized"
n_perm <- 1000
if("--mode" %in% args) {
  idx <- which(args == "--mode")
  if(length(idx) > 0 && idx[1] < length(args)) {
    mode_in <- tolower(args[idx[1] + 1])
    if(mode_in %in% c("legacy", "stabilized")) {
      mode <- mode_in
    } else {
      stop("Invalid --mode. Use 'legacy' or 'stabilized'.")
    }
  }
}
if("--n-perm" %in% args) {
  idx <- which(args == "--n-perm")
  if(length(idx) > 0 && idx[1] < length(args)) {
    n_perm <- as.integer(args[idx[1] + 1])
  }
}
if(is.na(n_perm) || n_perm < 10) {
  stop("--n-perm must be an integer >= 10.")
}

old_skip_main <- Sys.getenv("CADE_SKIP_MAIN", unset = NA)
Sys.setenv(CADE_SKIP_MAIN = "1")
source(file.path(PROJECT_ROOT, "code", "cade_method.R"))
if (is.na(old_skip_main)) {
  Sys.unsetenv("CADE_SKIP_MAIN")
} else {
  Sys.setenv(CADE_SKIP_MAIN = old_skip_main)
}

OUT_DIR <- file.path(PROJECT_ROOT, "analysis_output", "CADE")

cat(sprintf("=== CCI Permutation Test v2 (Two-sided, mode=%s) ===\n", mode))

# Load corrected expression matrix
exprs_mat <- as.matrix(read.csv(
  file.path(PROJECT_ROOT, "geo_analysis_output", "GSE26050_expression_corrected.csv"),
  row.names=1, check.names=FALSE
))
cat(sprintf("Loaded: %d genes x %d samples\n", nrow(exprs_mat), ncol(exprs_mat)))

# True group
true_group <- colnames(exprs_mat) %in% c(
  "GSM639703", "GSM639704", "GSM639705", "GSM639706", "GSM639707",
  "GSM639708", "GSM639709", "GSM639710", "GSM639711", "GSM639712", "GSM639713"
)

# Marker gene sets
marker_list <- list(
  CD8_Tcells = c("CD8A", "CD8B", "PRF1", "GZMB", "GZMA", "GNLY", "NKG7", "CD3E", "CD3D"),
  CD4_Tcells = c("CD4", "CD3E", "CD3D", "IL7R", "CCR7", "LEF1", "TCF7", "SELL"),
  NK_cells = c("NKG7", "GNLY", "PRF1", "GZMB", "KLRB1", "KLRD1", "KLRF1", "NCR1", "CD160"),
  B_cells = c("CD19", "CD79A", "CD79B", "MS4A1", "PAX5", "BLK", "BANK1", "CD22"),
  Monocytes = c("CD14", "FCGR3A", "CSF1R", "ITGAM", "LYZ", "S100A8", "S100A9", "VCAN"),
  Macrophages = c("CD68", "CD163", "MRC1", "MSR1", "CSF1R", "ITGAM", "TLR2", "TLR4"),
  Neutrophils = c("FCGR3B", "CXCR2", "CXCL8", "CSF3R", "MMP8", "MMP9", "ELANE", "MPO"),
  Erythrocytes = c("HBB", "HBA1", "HBA2", "GYPA", "ALAS2", "CA1", "SLC4A1", "AHSP")
)

ferroptosis_genes <- c("SLC7A11", "SLC25A37", "FTH1", "FTL", "GPX4",
                        "GCLM", "TFRC", "HMOX1", "SLC40A1", "NCOA4",
                        "STAT3", "JAK2", "IFNG", "STAT1", "NFE2L2",
                        "IL1B", "TNF", "IL6", "CXCL8", "NFKB1")
ferroptosis_genes <- intersect(ferroptosis_genes, rownames(exprs_mat))

# ── Run true CADE (for observed CCI) ──
cat("\nRunning true CADE...\n")
true_prop <- estimate_proportions_cade(exprs_mat, marker_list, max_iter=50, tol=1e-6)
true_de <- cade_de_analysis(exprs_mat, true_group, true_prop$proportions,
                            top_cts=4, cci_variant=mode)
rownames(true_de) <- true_de$Gene

obs_cci <- true_de[ferroptosis_genes, "CCI"]
names(obs_cci) <- ferroptosis_genes
cat(sprintf("  True CCI for SLC7A11: %.4f, GPX4: %.4f, FTL: %.4f\n",
    obs_cci["SLC7A11"], obs_cci["GPX4"], obs_cci["FTL"]))

# ── Optimization: estimate proportions ONCE (group-independent) ──
# CADE marker-derived weight estimation depends only on expression matrix + marker genes,
# not on group labels. Re-estimating per permutation is redundant.
cat("\nEstimating marker-derived weights (shared across all permutations)...\n")
shared_prop <- estimate_proportions_cade(exprs_mat, marker_list, max_iter=50, tol=1e-6)
cat(sprintf("  Converged at iteration %d (delta=%.2e)\n",
    shared_prop$n_iter, tail(shared_prop$convergence, 1)))

# ── Permutation test: only shuffle labels in DE step ──
set.seed(42)
cat(sprintf("\nRunning %d permutations (DE step only)...\n", n_perm))

perm_cci <- matrix(NA, nrow=length(ferroptosis_genes), ncol=n_perm,
                   dimnames=list(ferroptosis_genes, NULL))

pb <- txtProgressBar(min=0, max=n_perm, style=3)
for(p in 1:n_perm) {
  perm_group <- sample(true_group)

  perm_de <- tryCatch({
    cade_de_analysis(exprs_mat, perm_group, shared_prop$proportions,
                     top_cts=4, cci_variant=mode, verbose=FALSE)
  }, error=function(e) NULL)

  if(!is.null(perm_de)) {
    rownames(perm_de) <- perm_de$Gene
    perm_cci[, p] <- perm_de[ferroptosis_genes, "CCI"]
  }
  setTxtProgressBar(pb, p)
}
close(pb)

# ── Calculate two-sided empirical P-values ──
cat("\n\nCalculating two-sided empirical P-values...\n")
cci_results <- data.frame(
  Gene = ferroptosis_genes,
  CCI_observed = obs_cci,
  CCI_null_mean = rowMeans(perm_cci, na.rm=TRUE),
  CCI_null_sd = apply(perm_cci, 1, sd, na.rm=TRUE),
  CCI_null_median = apply(perm_cci, 1, median, na.rm=TRUE),
  stringsAsFactors = FALSE
)

# Compute both directional P-values
cci_results$P_lower <- sapply(ferroptosis_genes, function(g) {
  obs <- obs_cci[g]
  null_dist <- perm_cci[g, !is.na(perm_cci[g, ])]
  if(length(null_dist) < 10) return(NA)
  mean(null_dist <= obs, na.rm=TRUE)  # P(CCI_perm <= CCI_obs)
})

cci_results$P_upper <- sapply(ferroptosis_genes, function(g) {
  obs <- obs_cci[g]
  null_dist <- perm_cci[g, !is.na(perm_cci[g, ])]
  if(length(null_dist) < 10) return(NA)
  mean(null_dist >= obs, na.rm=TRUE)  # P(CCI_perm >= CCI_obs)
})

# Two-sided P-value
cci_results$P_two_sided <- with(cci_results, {
  pmin(2 * pmin(P_lower, P_upper), 1.0)
})

# FDR correction on two-sided P-values
cci_results$FDR_two_sided <- p.adjust(cci_results$P_two_sided, method="BH")

# -- Statistical classification (direction-aware sensitivity labels) --
cci_results$Direction <- with(cci_results, {
  dir <- rep("null", nrow(cci_results))
  dir[CCI_observed < CCI_null_mean] <- "lower"
  dir[CCI_observed > CCI_null_mean] <- "higher"
  dir[is.na(CCI_observed)] <- "undefined"
  dir
})

cci_results$Classification <- with(cci_results, {
  cls <- rep("Not significant", nrow(cci_results))
  # Significant AND CCI lower than null: unusually adjustment-stable coefficient
  cls[FDR_two_sided < 0.05 & Direction == "lower"] <-
    "Adjustment-stable outlier (FDR<0.05, CCI < null)"
  # Significant AND CCI higher than null: unusually adjustment-sensitive coefficient
  cls[FDR_two_sided < 0.05 & Direction == "higher"] <-
    "Adjustment-sensitive outlier (FDR<0.05, CCI > null)"
  cls[is.na(FDR_two_sided)] <- "Insufficient data"
  cls
})

# Also report heuristic CCI rank for continuity
cci_results$CCI_Rank <- with(cci_results, {
  r <- rep("Low-moderate", nrow(cci_results))
  r[CCI_observed < 0.2] <- "Lowest"
  r[CCI_observed > 0.5] <- "High"
  r[is.na(CCI_observed)] <- "Undefined"
  r
})

# ── Results ──
cat("\n=== CCI Permutation Test Results (Two-sided) ===\n")
cat(sprintf("%-12s %8s %8s %8s %8s %8s %8s %s\n",
    "Gene", "CCI_obs", "NullMean", "P_lower", "P_upper", "P_2s", "FDR_2s", "Classification"))
for(i in 1:nrow(cci_results)) {
  r <- cci_results[i,]
  cat(sprintf("%-12s %8.3f %8.3f %8.3f %8.3f %8.3f %8.4f %s\n",
      r$Gene, r$CCI_observed, r$CCI_null_mean,
      r$P_lower, r$P_upper, r$P_two_sided, r$FDR_two_sided, r$Classification))
}

# ── Save ──
file_suffix <- paste0("_", mode)
if(n_perm != 1000) {
  file_suffix <- paste0(file_suffix, "_test_n", n_perm)
}
out_file <- file.path(OUT_DIR, paste0("Table_CCI_Permutation_Test", file_suffix, ".csv"))
write.csv(cci_results, out_file, row.names=FALSE)
cat(sprintf("\nSaved: %s\n", out_file))

# ── Summary ──
n_stable <- sum(grepl("Adjustment-stable", cci_results$Classification))
n_sensitive <- sum(grepl("Adjustment-sensitive", cci_results$Classification))
n_null <- sum(grepl("Not significant", cci_results$Classification))
n_insuff <- sum(grepl("Insufficient", cci_results$Classification))
cat(sprintf("\nSummary: %d adjustment-stable outliers, %d adjustment-sensitive outliers, %d not significant, %d insufficient data\n",
    n_stable, n_sensitive, n_null, n_insuff))

# ── Directional summary ──
cat(sprintf("\nDirection: %d genes with CCI below null mean, %d above null mean\n",
    sum(cci_results$Direction == "lower", na.rm=TRUE),
    sum(cci_results$Direction == "higher", na.rm=TRUE)))

cat("=== Permutation test complete ===\n")
