#!/usr/bin/env Rscript
# CADE-ILR basis-stability analysis
#
# Tests whether the CADE-ILR CCI ranking depends on the choice of ILR basis
# (default Helmert vs. alternative contrasts). A basis-stable ranking
# supports the claim that the ILR sensitivity layer is robust.
#
# Usage:
#   Rscript code/cade_ilr_basis_stability.R [--out-dir <dir>]

suppressPackageStartupMessages({
  library(limma)
  library(quadprog)
})

get_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  hit <- grep(file_arg, args_full, value = TRUE)
  if(length(hit) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", hit[1]))))
  }
  getwd()
}
script_dir <- get_script_dir()
PROJECT_ROOT <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
DATA_ROOT <- normalizePath(file.path(PROJECT_ROOT, ".."), mustWork = FALSE)
if(!dir.exists(file.path(DATA_ROOT, "geo_analysis_output"))) {
  DATA_ROOT <- PROJECT_ROOT
}
out_dir <- file.path(PROJECT_ROOT, "analysis_output", "CADE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Project root: %s\n", PROJECT_ROOT))
cat(sprintf("Data root: %s\n", DATA_ROOT))

# Load CADE functions (lines 1-728)
cat("Loading CADE functions...\n")
cm_content <- readLines(file.path(PROJECT_ROOT, "code", "cade_method.R"))
fn_content <- cm_content[1:728]
clean_path <- tempfile(fileext = ".R")
writeLines(fn_content, clean_path)
source(clean_path, echo=FALSE)
unlink(clean_path)
cat("  Functions loaded\n")

# Load data
cat("Loading GSE26050 expression data...\n")
expr_file <- file.path(DATA_ROOT, "geo_analysis_output", "GSE26050_expression_corrected.csv")
if(!file.exists(expr_file)) {
  stop("Cannot find expression file: ", expr_file)
}
exprs_mat <- as.matrix(read.csv(expr_file, row.names=1, check.names=FALSE))
true_group <- c(rep(1, 11), rep(0, 33))
cat(sprintf("  Loaded: %d genes × %d samples\n", nrow(exprs_mat), ncol(exprs_mat)))

# Marker list
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

# Ferroptosis panel
ferroptosis_genes <- c("SLC7A11","IFNG","TFRC","FTH1","TNF","IL1B","NFKB1","FTL",
                        "STAT3","JAK2","NFE2L2","SLC25A37","CXCL8","SLC40A1","IL6",
                        "GPX4","GCLM","HMOX1","STAT1","NCOA4")
ferroptosis_genes <- intersect(ferroptosis_genes, rownames(exprs_mat))

# ── Estimate weights once ──
cat("\nEstimating weights...\n")
prop <- estimate_proportions_cade(exprs_mat, marker_list, max_iter=30, tol=1e-5, verbose=FALSE)
cat(sprintf("  Converged at iter %d\n", prop$n_iter))

# ── Define different ILR bases to test ──
# Default: contr.helmert (in cade_ilr_transform)
# Alternative 1: contr.helmert with different reference (manually constructed)
# Alternative 2: contr.treatment (treatment contrasts)
# Alternative 3: Log-ratio with first component as reference (geometric mean of rest)
n_comp <- nrow(prop$proportions)
cat(sprintf("Number of components: %d\n", n_comp))

# Default Helmert
basis_default <- contr.helmert(n_comp)
basis_default <- sweep(basis_default, 2, sqrt(colSums(basis_default^2)), "/")

# Alternative 1: Inverse Helmert (start from last component as reference)
basis_inv_helmert <- -basis_default
basis_inv_helmert[1, ] <- -basis_inv_helmert[1, ]  # Flip to maintain orthonormality
# (proper inverse helmert construction)
basis_inv_helmert <- matrix(0, n_comp, n_comp-1)
for(j in 1:(n_comp-1)) {
  basis_inv_helmert[1:(n_comp-j), j] <- 1
  basis_inv_helmert[(n_comp-j+1):n_comp, j] <- -(n_comp-j)
  basis_inv_helmert[, j] <- basis_inv_helmert[, j] / sqrt((n_comp-j)*(n_comp-j+1))
}

# Alternative 2: Forward differences (sequential pairwise)
basis_fwd <- matrix(0, n_comp, n_comp-1)
for(j in 1:(n_comp-1)) {
  basis_fwd[j, j] <- 1
  basis_fwd[j+1, j] <- -1
  basis_fwd[, j] <- basis_fwd[, j] / sqrt(2)
}

# Alternative 3: First-vs-rest (each component vs. the rest)
basis_fvr <- matrix(0, n_comp, n_comp-1)
for(j in 1:(n_comp-1)) {
  basis_fvr[j, j] <- 1
  basis_fvr[(j+1):n_comp, j] <- -1/(n_comp-1)
  norm <- sqrt(1 + (n_comp-1)*(1/(n_comp-1))^2)
  basis_fvr[, j] <- basis_fvr[, j] / norm
}

bases <- list(
  default_helmert = basis_default,
  inverse_helmert = basis_inv_helmert,
  forward_differences = basis_fwd,
  first_vs_rest = basis_fvr
)

# ── Compute CCI for each basis ──
cat("\n=== Computing CCI for each ILR basis ===\n")
cci_results <- data.frame(Gene = ferroptosis_genes, stringsAsFactors=FALSE)
for(basis_name in names(bases)) {
  cat(sprintf("  Basis: %s ... ", basis_name))

  basis <- bases[[basis_name]]
  # Use the helper function with a custom basis
  closed <- close_composition(prop$proportions, pseudocount=1e-6)
  coords <- t(log(t(closed)) %*% basis)
  rownames(coords) <- colnames(basis)
  colnames(coords) <- colnames(closed)

  # Compute top-variance ILR coordinates
  coord_var <- apply(coords, 1, var)
  top_k <- min(4, nrow(coords))
  top_idx <- order(coord_var, decreasing=TRUE)[1:top_k]
  selected_coords <- t(coords[top_idx, , drop=FALSE])

  # Run DE with ILR covariates using limma
  # selected_coords is currently 44 samples × 4 ILR coords (samples in rows)
  # We need to attach the group label to the SAMPLE dimension
  group_vec <- as.factor(true_group)
  # Build design as data.frame with sample IDs as rownames
  design_data <- data.frame(
    row.names = rownames(selected_coords),
    group = group_vec
  )
  # Add each ILR coordinate explicitly (4 of them)
  for(j in 1:min(4, ncol(selected_coords))) {
    design_data[[paste0("ILR", j)]] <- selected_coords[, j]
  }

  # Adjusted model
  design_full <- model.matrix(~ group + ., data = design_data)
  fit <- suppressWarnings(lmFit(exprs_mat[ferroptosis_genes, ], design_full))
  logfc_adj_ilr <- coef(fit)[, "group1"]

  # Unadjusted model
  design_unadj <- model.matrix(~ group, data = design_data)
  fit_unadj <- suppressWarnings(lmFit(exprs_mat[ferroptosis_genes, ], design_unadj))
  logfc_unadj_vec <- coef(fit_unadj)[, "group1"]

  cci_ilr <- abs(logfc_unadj_vec - logfc_adj_ilr) / pmax(0.1, abs(logfc_unadj_vec))
  cci_ilr <- pmin(1, cci_ilr)
  cci_results[[paste0("CCI_", basis_name)]] <- cci_ilr
  cat("done\n")
}

# ── Compute rank stability across bases ──
cat("\n=== Rank stability analysis ===\n")
rank_results <- data.frame(Gene = ferroptosis_genes, stringsAsFactors=FALSE)
for(basis_name in names(bases)) {
  col <- paste0("CCI_", basis_name)
  rank_results[[paste0("Rank_", basis_name)]] <- rank(-cci_results[[col]])  # Negative for descending
}
# Compute Spearman correlation between each pair of bases (use matrix form to avoid NA issue)
rank_mat <- as.matrix(rank_results[, -1])
cat("\nSpearman correlations between bases:\n")
rho_mat <- cor(rank_mat, method="spearman", use="complete.obs")
for(i in 1:(ncol(rank_mat)-1)) {
  for(j in (i+1):ncol(rank_mat)) {
    cat(sprintf("  %s vs %s: rho = %.3f\n", colnames(rank_mat)[i], colnames(rank_mat)[j], rho_mat[i,j]))
  }
}

# ── Mean and SD across bases ──
cci_cols <- grep("^CCI_", names(cci_results), value=TRUE)
cci_results$CCI_mean <- rowMeans(cci_results[, cci_cols])
cci_results$CCI_sd <- apply(cci_results[, cci_cols], 1, sd)
cci_results$CCI_min <- apply(cci_results[, cci_cols], 1, min)
cci_results$CCI_max <- apply(cci_results[, cci_cols], 1, max)
cci_results$CCI_range <- cci_results$CCI_max - cci_results$CCI_min

# Save outputs
out_file <- file.path(out_dir, "Table_S_ILR_BasisStability.csv")
write.csv(cci_results, out_file, row.names=FALSE)
cat(sprintf("\nSaved: %s\n", out_file))

# Summary
cat("\n=== Summary ===\n")
cat(sprintf("Mean coefficient of variation (SD/mean) across bases: %.3f\n",
            mean(cci_results$CCI_sd / pmax(0.01, cci_results$CCI_mean), na.rm=TRUE)))
cat(sprintf("Mean CCI range (max-min) across bases: %.3f\n", mean(cci_results$CCI_range)))
cat(sprintf("Spearman correlation with default: all > 0.85 indicates basis stability\n"))

cat("\n=== ILR basis stability analysis complete ===\n")
