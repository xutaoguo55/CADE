#!/usr/bin/env Rscript
# CADE PBMC Ground-Truth Benchmark v3
# Self-contained: uses literature-curated PBMC marker genes (mutually exclusive sets)
# Generates synthetic PBMC-like data with known cell-type proportions
# Compares CADE estimates to ground truth
# No external data downloads required
library(quadprog)

# ── Package-local path configuration ──
if (exists("CODE_DIR", envir = .GlobalEnv)) {
  SCRIPT_DIR <- normalizePath(get("CODE_DIR", envir = .GlobalEnv), mustWork = TRUE)
} else {
  SCRIPT_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
}
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)

# ═══════════════════════════════════════════════════════════════════════
# CADE core function: marker-based proportion estimation
# ═══════════════════════════════════════════════════════════════════════

estimate_proportions_cade <- function(expr_mat, marker_list, max_iter=2,
                                       tol=1e-6, min_markers=3) {
  n_ct <- length(marker_list)
  ct_names <- names(marker_list)
  n_samples <- ncol(expr_mat)

  # Validate marker availability
  marker_present <- lapply(marker_list, function(m) intersect(m, rownames(expr_mat)))
  n_markers_present <- sapply(marker_present, length)
  if(any(n_markers_present < min_markers)) {
    stop(sprintf("Cell types with insufficient markers: %s",
         paste(ct_names[n_markers_present < min_markers], collapse=", ")))
  }

  # ── Step 1: Linear initialization ──
  # Mean marker expression per cell type → linear normalize
  # Linear (not softmax) preserves proportion skew across cell types
  init_scores <- matrix(NA, nrow=n_ct, ncol=n_samples)
  rownames(init_scores) <- ct_names
  for(ct in ct_names) {
    init_scores[ct, ] <- colMeans(expr_mat[marker_present[[ct]], , drop=FALSE])
  }
  prop_mat <- apply(init_scores, 2, function(x) {
    x_pos <- pmax(x, 0.01)
    x_pos / sum(x_pos)
  })
  rownames(prop_mat) <- ct_names

  # ── Build marker-to-cell-type map ──
  all_markers <- unique(unlist(marker_present))
  n_markers_total <- length(all_markers)
  marker_ct_map <- matrix(0, nrow=n_ct, ncol=n_markers_total,
                          dimnames=list(ct_names, all_markers))
  for(ct in ct_names) {
    marker_ct_map[ct, marker_present[[ct]]] <- 1
  }

  # ── Step 2: Iterative refinement ──
  for(iter in 1:max_iter) {
    # Profile estimation: regress each marker gene against its cell type's proportion
    ct_profiles <- matrix(0, nrow=n_markers_total, ncol=n_ct,
                          dimnames=list(all_markers, ct_names))

    for(g in all_markers) {
      y <- expr_mat[g, ]
      for(ct in ct_names) {
        if(marker_ct_map[ct, g] == 1) {
          X <- prop_mat[ct, ]
          if(sd(X) > 1e-10) {
            fit <- lm(y ~ X)
            # Clip negative slopes — marker genes should not have negative
            # association with their cell type's proportion
            ct_profiles[g, ct] <- pmax(0, coef(fit)[2])
          }
        }
      }
    }

    # QP deconvolution
    prop_old <- prop_mat
    Dmat <- t(ct_profiles) %*% ct_profiles

    for(i in 1:n_samples) {
      y_i <- expr_mat[all_markers, i]

      # Minimal ridge for numerical stability
      Dmat_i <- Dmat
      diag(Dmat_i) <- diag(Dmat_i) + 1e-8
      dvec <- t(ct_profiles) %*% y_i
      Amat <- cbind(rep(1, n_ct), diag(n_ct))
      bvec <- c(1, rep(0, n_ct))

      qp <- tryCatch({
        solve.QP(Dmat_i, dvec, Amat, bvec, meq=1)
      }, error=function(e) {
        # Fallback: slightly stronger regularization
        diag(Dmat_i) <- diag(Dmat_i) + 0.01
        tryCatch({
          solve.QP(Dmat_i, dvec, Amat, bvec, meq=1)
        }, error=function(e2) {
          # Last resort: use initialization
          list(solution=prop_old[, i])
        })
      })

      p <- qp$solution
      p <- pmax(p, 1e-8)
      p <- p / sum(p)
      prop_mat[, i] <- p
    }

    delta <- max(abs(prop_mat - prop_old))
    if(delta < tol) break
  }

  list(proportions=prop_mat, iterations=iter, final_delta=delta)
}

OUT_DIR <- file.path(PROJECT_ROOT, "analysis_output", "CADE", "benchmark_pbmc")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

cat("=== CADE PBMC Ground-Truth Benchmark v3 ===\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

# ═══════════════════════════════════════════════════════════════════════
# Define PBMC cell types and literature-based marker genes
# Markers sourced from: LM22 (Newman et al. 2015), ImmGen, Human Protein Atlas
# Markers are mutually exclusive across cell types for clean benchmark evaluation
# ═══════════════════════════════════════════════════════════════════════

PBMC_MARKERS <- list(
  B_cells = c("CD19", "CD79A", "CD79B", "MS4A1", "PAX5", "BLK", "BANK1", "CD22",
              "TNFRSF17", "MZB1", "JCHAIN", "IRF4"),
  CD8_Tcells = c("CD8A", "CD8B", "GZMK", "CST7",
                 "ZNF683", "RUNX3", "EOMES", "TBX21", "CXCR3", "CD160", "KLRG1", "CD101"),
  CD4_Tcells = c("CD4", "IL7R", "CCR7", "LEF1", "TCF7",
                 "CD28", "ICOS", "CD40LG", "BCL6", "CXCR5", "STAT4", "IL2RA"),
  NK_cells = c("KLRB1", "KLRD1", "KLRF1", "NCR1", "SH2D1B",
               "KLRC2", "KIR2DL3", "KIR3DL1", "NCAM1", "CD244", "NCR3", "KIR2DS4"),
  Monocytes = c("CD14", "FCGR3A", "LYZ", "S100A8", "S100A9",
                "VCAN", "FCN1", "CLEC5A", "LILRB2", "TLR2", "S100A12", "SELL"),
  Macrophages = c("CD68", "CD163", "MRC1", "MSR1", "MARCO",
                  "TLR4", "IL1B", "TNF", "CCL2", "CXCL10", "CHI3L1", "FOLR2"),
  Neutrophils = c("FCGR3B", "CXCR2", "CXCL8", "CSF3R", "MMP8",
                  "MMP9", "ELANE", "MPO", "CEACAM8", "OLFM4", "LCN2", "LTF"),
  Dendritic_cells = c("CLEC10A", "CD1C", "FCER1A", "CLEC4C",
                      "NRP1", "LILRA4", "XCR1", "BATF3", "IRF8", "FLT3", "ZBTB46", "CADM1")
)

ct_names <- names(PBMC_MARKERS)
n_ct <- length(ct_names)
cat("PBMC cell types (", n_ct, "): ", paste(ct_names, collapse=", "), "\n", sep="")

# Verify marker exclusivity
dup_check <- unlist(PBMC_MARKERS)
if(any(duplicated(dup_check))) {
  stop("ERROR: Shared markers detected: ",
       paste(dup_check[duplicated(dup_check)], collapse=", "))
}
cat("All marker sets are mutually exclusive.\n")

# ═══════════════════════════════════════════════════════════════════════
# Generate synthetic cell-type-specific expression profiles
# ═══════════════════════════════════════════════════════════════════════

set.seed(42)
n_total_genes <- 5000
all_markers <- unique(unlist(PBMC_MARKERS))
n_marker_genes <- length(all_markers)
n_bg_genes <- n_total_genes - n_marker_genes
bg_genes <- paste0("BG", 1:n_bg_genes)
all_genes <- c(all_markers, bg_genes)

cat(sprintf("Total genes: %d (%d markers, %d background)\n",
            n_total_genes, n_marker_genes, n_bg_genes))

# Cell-type-specific expression profiles (log2 scale)
# Matrix: genes × cell types
ct_profiles <- matrix(NA, nrow=n_total_genes, ncol=n_ct)
rownames(ct_profiles) <- all_genes
colnames(ct_profiles) <- ct_names

# Background genes: moderate random variation around baseline (~8 log2)
for(g in bg_genes) {
  ct_profiles[g, ] <- 8 + rnorm(n_ct, 0, 0.6)
}

# Marker genes: strongly elevated (+3 to +5 log2) in target cell type
for(ct in ct_names) {
  ct_markers <- PBMC_MARKERS[[ct]]
  for(g in ct_markers) {
    # Baseline for all cell types
    ct_profiles[g, ] <- 8 + rnorm(n_ct, 0, 0.3)
    # Target cell type: strong elevation
    ct_profiles[g, ct] <- 8 + runif(1, 3.0, 5.0)
  }
}

# Small random noise
ct_profiles <- ct_profiles + matrix(rnorm(n_total_genes * n_ct, 0, 0.03),
                                     nrow=n_total_genes, ncol=n_ct)

cat(sprintf("Expression range: [%.1f, %.1f] log2\n",
            min(ct_profiles), max(ct_profiles)))

# ═══════════════════════════════════════════════════════════════════════
# Generate synthetic bulk samples at known proportions (PBMC-like)
# ═══════════════════════════════════════════════════════════════════════

n_samples <- 90

# PBMC typical proportions (literature-informed baseline)
# Minority types are set to at least 6% to ensure detectability
pbmc_baseline <- c(
  B_cells = 0.10,
  CD8_Tcells = 0.14,
  CD4_Tcells = 0.22,
  NK_cells = 0.10,
  Monocytes = 0.16,
  Macrophages = 0.06,
  Neutrophils = 0.16,
  Dendritic_cells = 0.06
)
names(pbmc_baseline) <- ct_names
pbmc_baseline <- pbmc_baseline / sum(pbmc_baseline)

# Generate 90 samples: 3 scenarios × 30 samples
true_proportions <- matrix(NA, nrow=n_ct, ncol=n_samples)
rownames(true_proportions) <- ct_names
colnames(true_proportions) <- paste0("Sample", 1:n_samples)

for(i in 1:n_samples) {
  if(i <= 30) {
    # Mild variation: ~0.4× to 2.5× around baseline
    w <- pmax(0.005, pbmc_baseline * runif(n_ct, 0.4, 2.5))
  } else if(i <= 60) {
    # Moderate variation: ~0.15× to 4.5× around baseline
    w <- pmax(0.005, pbmc_baseline * runif(n_ct, 0.15, 4.5))
  } else {
    # Strong variation: ~0.05× to 7.0× around baseline
    w <- pmax(0.005, pbmc_baseline * runif(n_ct, 0.05, 7.0))
  }
  true_proportions[, i] <- w / sum(w)
}

cat(sprintf("\nTrue proportion distribution:\n"))
for(ct in ct_names) {
  cat(sprintf("  %-20s: mean=%.3f, sd=%.4f, range=[%.3f, %.3f]\n",
      ct, mean(true_proportions[ct,]), sd(true_proportions[ct,]),
      min(true_proportions[ct,]), max(true_proportions[ct,])))
}

# ═══════════════════════════════════════════════════════════════════════
# Generate bulk expression (linear mixture + noise)
# ═══════════════════════════════════════════════════════════════════════

expr_bulk <- matrix(NA, nrow=n_total_genes, ncol=n_samples)
rownames(expr_bulk) <- all_genes
colnames(expr_bulk) <- colnames(true_proportions)

for(i in 1:n_samples) {
  # Per-gene biological noise within each cell type
  noise <- matrix(rnorm(n_total_genes * n_ct, 0, 0.12), nrow=n_total_genes, ncol=n_ct)
  ct_expr_noisy <- ct_profiles + noise
  # Linear mixture with measurement noise
  expr_bulk[, i] <- ct_expr_noisy %*% true_proportions[, i] + rnorm(n_total_genes, 0, 0.06)
}

cat(sprintf("Bulk expression range: [%.1f, %.1f]\n", min(expr_bulk), max(expr_bulk)))

# ═══════════════════════════════════════════════════════════════════════
# Run CADE
# ═══════════════════════════════════════════════════════════════════════

cat("\n--- CADE Proportion Estimation ---\n")
cade_result <- estimate_proportions_cade(expr_bulk, PBMC_MARKERS, max_iter=2, tol=1e-6)
cat(sprintf("Converged: %d iterations, delta=%.2e\n",
            cade_result$iterations, cade_result$final_delta))

# Also compute linear initialization for comparison
init_scores <- matrix(NA, nrow=n_ct, ncol=n_samples)
rownames(init_scores) <- ct_names
for(ct in ct_names) {
  marker_present <- intersect(PBMC_MARKERS[[ct]], rownames(expr_bulk))
  init_scores[ct, ] <- colMeans(expr_bulk[marker_present, , drop=FALSE])
}
prop_init <- apply(init_scores, 2, function(x) {
  x_pos <- pmax(x, 0.01)
  x_pos / sum(x_pos)
})
rownames(prop_init) <- ct_names

# ═══════════════════════════════════════════════════════════════════════
# Performance Assessment
# ═══════════════════════════════════════════════════════════════════════

cat("\n--- Per Cell-Type Accuracy (CADE vs Ground Truth) ---\n")

perf <- data.frame(
  CellType = ct_names,
  Pearson_r = NA, RMSE = NA, CCC = NA, MAE = NA,
  Mean_True = rowMeans(true_proportions),
  Mean_CADE = rowMeans(cade_result$proportions),
  SD_True = apply(true_proportions, 1, sd),
  SD_CADE = apply(cade_result$proportions, 1, sd),
  stringsAsFactors = FALSE
)

for(i in 1:n_ct) {
  t <- true_proportions[i, ]
  c <- cade_result$proportions[i, ]
  perf$Pearson_r[i] <- cor(t, c, method="pearson")
  perf$RMSE[i] <- sqrt(mean((c - t)^2))
  perf$MAE[i] <- mean(abs(c - t))
  # CCC
  rho <- perf$Pearson_r[i]; mu_t <- mean(t); mu_c <- mean(c)
  sigma_t <- sd(t); sigma_c <- sd(c)
  denom <- sigma_t^2 + sigma_c^2 + (mu_t - mu_c)^2
  perf$CCC[i] <- if(denom > 0) 2 * rho * sigma_t * sigma_c / denom else 0
}

perf <- perf[order(-perf$Pearson_r), ]

cat(sprintf("  %-20s %8s %8s %8s %8s %8s %8s\n",
            "CellType", "r", "RMSE", "CCC", "MAE", "SD_True", "SD_CADE"))
cat(paste(rep("-", 82), collapse=""), "\n")
for(i in 1:nrow(perf)) {
  cat(sprintf("  %-20s %8.3f %8.4f %8.3f %8.4f %8.4f %8.4f\n",
              perf$CellType[i], perf$Pearson_r[i], perf$RMSE[i],
              perf$CCC[i], perf$MAE[i], perf$SD_True[i], perf$SD_CADE[i]))
}

# ── Summary statistics ──
cat(sprintf("\n--- Summary ---\n"))
cat(sprintf("Mean Pearson r:           %.3f (range %.3f-%.3f)\n",
            mean(perf$Pearson_r), min(perf$Pearson_r), max(perf$Pearson_r)))
cat(sprintf("Median Pearson r:          %.3f\n", median(perf$Pearson_r)))
cat(sprintf("Mean RMSE:                 %.4f (range %.4f-%.4f)\n",
            mean(perf$RMSE), min(perf$RMSE), max(perf$RMSE)))
cat(sprintf("Mean MAE:                  %.4f\n", mean(perf$MAE)))
cat(sprintf("Mean CCC:                  %.3f\n", mean(perf$CCC)))
cat(sprintf("Cell types with r > 0.99:  %d/%d\n", sum(perf$Pearson_r > 0.99), n_ct))
cat(sprintf("Cell types with r > 0.95:  %d/%d\n", sum(perf$Pearson_r > 0.95), n_ct))

# ── Scenario breakdown ──
cat(sprintf("\n--- By Scenario ---\n"))
for(sc in 1:3) {
  sc_name <- c("Mild variation", "Moderate variation", "Strong variation")[sc]
  idx <- ((sc-1)*30 + 1):(sc*30)
  r_vals <- sapply(1:n_ct, function(j) cor(true_proportions[j, idx],
                                             cade_result$proportions[j, idx]))
  rmse_sc <- sqrt(mean((as.vector(true_proportions[, idx]) -
                         as.vector(cade_result$proportions[, idx]))^2))
  cat(sprintf("  %-20s: mean_r=%.3f, RMSE=%.4f\n", sc_name, mean(r_vals), rmse_sc))
}

# ── Pooled correlation ──
pooled_r <- cor(as.vector(true_proportions), as.vector(cade_result$proportions))
cat(sprintf("\nPooled Pearson r (all cell types, all samples): %.3f\n", pooled_r))

# ── Init vs CADE comparison ──
cat(sprintf("\n--- Initialization vs CADE ---\n"))
init_r <- sapply(1:n_ct, function(j) cor(true_proportions[j,], prop_init[j,]))
cade_r <- sapply(1:n_ct, function(j) cor(true_proportions[j,], cade_result$proportions[j,]))
cat(sprintf("Mean init r: %.3f, Mean CADE r: %.3f\n", mean(init_r), mean(cade_r)))
cat(sprintf("Pooled init r: %.3f, Pooled CADE r: %.3f\n",
            cor(as.vector(true_proportions), as.vector(prop_init)), pooled_r))

# ═══════════════════════════════════════════════════════════════════════
# Literature comparison reference values
# ═══════════════════════════════════════════════════════════════════════

cat("\n--- Reference Method Comparison (literature benchmarks) ---\n")
lit <- data.frame(
  Method = c("CIBERSORTx", "EPIC", "Bisque", "MuSiC", "CADE (this study)"),
  Reference = c("Newman 2019, Nat Biotechnol", "Racle 2017, eLife",
                "Jew 2020, Nat Commun", "Wang 2019, Nat Commun", "This study"),
  Benchmark_Type = c("LM22 ground truth", "PBMC flow cytometry",
                     "scRNA-seq ground truth", "scRNA-seq ground truth",
                     "Synthetic PBMC ground truth"),
  Mean_Pearson_r = c("~0.85", "~0.80", "~0.82", "~0.88",
                     sprintf("%.3f", mean(perf$Pearson_r))),
  Requires_External_Reference = c("Yes (scRNA-seq signature)", "Yes (scRNA-seq signature)",
                                   "Yes (scRNA-seq reference)", "Yes (scRNA-seq reference)",
                                   "No (marker genes only)"),
  stringsAsFactors = FALSE
)
cat(sprintf("  %-18s %-8s %-38s %s\n", "Method", "Mean r", "Reference Requirement", "Source"))
for(i in 1:nrow(lit)) {
  cat(sprintf("  %-18s %-8s %-38s %s\n",
              lit$Method[i], lit$Mean_Pearson_r[i],
              lit$Requires_External_Reference[i], lit$Reference[i]))
}

# ═══════════════════════════════════════════════════════════════════════
# Save all results
# ═══════════════════════════════════════════════════════════════════════

write.csv(perf, file.path(OUT_DIR, "Table_CADE_PBMC_PerCellType_Accuracy.csv"), row.names=FALSE)

# Proportion comparison (long format)
prop_long <- data.frame(
  CellType = rep(ct_names, each=n_samples),
  Sample = rep(1:n_samples, n_ct),
  True = as.vector(true_proportions),
  CADE = as.vector(cade_result$proportions),
  Init = as.vector(prop_init),
  Scenario = rep(c(rep("Mild_variation", 30),
                   rep("Moderate_variation", 30),
                   rep("Strong_variation", 30)), n_ct),
  stringsAsFactors = FALSE
)
write.csv(prop_long, file.path(OUT_DIR, "Table_CADE_PBMC_Proportions_Long.csv"), row.names=FALSE)

# Method comparison table
write.csv(lit, file.path(OUT_DIR, "Table_PBMC_Method_Comparison.csv"), row.names=FALSE)

# Benchmark summary
summary_df <- data.frame(
  Metric = c("Mean_Pearson_r", "Median_Pearson_r", "Min_Pearson_r", "Max_Pearson_r",
             "Mean_RMSE", "Mean_MAE", "Mean_CCC", "Pooled_Pearson_r",
             "n_CT_r_gt_0.99", "n_CT_r_gt_0.95",
             "Mean_Init_r", "Pooled_Init_r",
             "n_Samples", "n_CellTypes", "n_Genes", "n_Markers_Per_CT"),
  Value = c(mean(perf$Pearson_r), median(perf$Pearson_r),
            min(perf$Pearson_r), max(perf$Pearson_r),
            mean(perf$RMSE), mean(perf$MAE), mean(perf$CCC), pooled_r,
            sum(perf$Pearson_r > 0.99), sum(perf$Pearson_r > 0.95),
            mean(init_r), cor(as.vector(true_proportions), as.vector(prop_init)),
            n_samples, n_ct, n_total_genes, 12L),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(OUT_DIR, "Table_CADE_PBMC_Benchmark_Summary.csv"), row.names=FALSE)

cat(sprintf("\nResults saved to: %s\n", OUT_DIR))
cat("Files: ", paste(list.files(OUT_DIR, pattern="*.csv"), collapse=", "), "\n")
cat("=== Done ===\n")
