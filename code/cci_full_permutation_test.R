#!/usr/bin/env Rscript
# CCI Full Permutation Test v1
# Re-estimates marker-derived weights at each shuffle, propagating weight-estimation
# uncertainty into the CCI null distribution. Tests whether the "shared-weight"
# approximation in cci_permutation_test.R materially affects the null.
#
# Usage:
#   Rscript code/cci_full_permutation_test.R [--n-perm 100] [--out-dir <dir>]

suppressPackageStartupMessages({
  library(limma)
  library(quadprog)
})

# ── Path configuration ──
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
# Look for data: check submission package first, then parent
DATA_ROOT <- normalizePath(file.path(PROJECT_ROOT, ".."), mustWork = FALSE)
if(!dir.exists(file.path(DATA_ROOT, "geo_analysis_output"))) {
  DATA_ROOT <- PROJECT_ROOT
}
cat(sprintf("Project root: %s\n", PROJECT_ROOT))
cat(sprintf("Data root: %s\n", DATA_ROOT))

# ── Parse arguments ──
n_perm <- 100
out_dir <- file.path(PROJECT_ROOT, "analysis_output", "CADE")
mode <- "stabilized"
seed <- 42
args <- commandArgs(trailingOnly = TRUE)
i <- 1
while(i <= length(args)) {
  if(args[i] == "--n-perm") { n_perm <- as.integer(args[i+1]); i <- i + 2
  } else if(args[i] == "--out-dir") { out_dir <- args[i+1]; i <- i + 2
  } else if(args[i] == "--seed") { seed <- as.integer(args[i+1]); i <- i + 2
  } else { i <- i + 1 }
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("n_perm = %d, out_dir = %s, seed = %d\n", n_perm, out_dir, seed))

# ── Load CADE functions (lines 1-728 only, no benchmark code) ──
cat("Loading CADE functions...\n")
cm_content <- readLines(file.path(PROJECT_ROOT, "code", "cade_method.R"))
fn_content <- cm_content[1:728]
clean_path <- tempfile(fileext = ".R")
writeLines(fn_content, clean_path)
source(clean_path, echo=FALSE)
unlink(clean_path)
cat("  Functions loaded\n")

# ── Load GSE26050 expression data ──
cat("Loading GSE26050 expression data...\n")
expr_file <- file.path(DATA_ROOT, "geo_analysis_output", "GSE26050_expression_corrected.csv")
if(!file.exists(expr_file)) {
  stop("Cannot find GSE26050_expression_corrected.csv at: ", expr_file)
}
exprs_mat <- as.matrix(read.csv(expr_file, row.names = 1, check.names = FALSE))
cat(sprintf("  Loaded: %d genes × %d samples\n", nrow(exprs_mat), ncol(exprs_mat)))

# FHL group assignment: GSE26050 has 11 FHL (GSM639703-13) + 33 control
# Per published metadata, first 11 = FHL, remaining 33 = control
n_fhl <- 11
true_group <- c(rep(1, n_fhl), rep(0, ncol(exprs_mat) - n_fhl))
cat(sprintf("  FHL=%d, Control=%d\n", sum(true_group==1), sum(true_group==0)))

# ── Marker list (8 sets, same as cci_permutation_test.R) ──
marker_list <- list(
  CD8_Tcells  = c("CD8A","CD8B","CD3D","CD3E","TRAC","CD2","GZMK","CCL5","PRF1"),
  CD4_Tcells  = c("CD4","IL7R","CCR7","LEF1","MAL","TCF7","LDHB"),
  NK_cells    = c("NKG7","GNLY","KLRD1","KLRB1","GZMB","CTSW"),
  B_cells     = c("CD19","MS4A1","CD79A","CD79B","BANK1","CD22","PAX5"),
  Monocytes   = c("LYZ","CD14","FCGR3A","MS4A7","ITGAM","CCR2","CD163","CSF1R","S100A8"),
  Macrophages = c("CD68","CD163","MRC1","MSR1","MARCO","CSF1R"),
  Neutrophils = c("FCGR3B","CSF3R","S100A8","S100A9","CXCR2","ITGAM","MMP9"),
  Erythrocytes= c("HBB","HBA1","HBA2","HBD","AHSP","ALAS2","SLC25A37")
)
marker_list <- lapply(marker_list, function(g) intersect(g, rownames(exprs_mat)))
marker_list <- marker_list[sapply(marker_list, length) >= 4]
cat(sprintf("Marker sets: %d\n", length(marker_list)))

# Ferroptosis panel
ferroptosis_genes <- c("SLC7A11","IFNG","TFRC","FTH1","TNF","IL1B","NFKB1","FTL",
                        "STAT3","JAK2","NFE2L2","SLC25A37","CXCL8","SLC40A1","IL6",
                        "GPX4","GCLM","HMOX1","STAT1","NCOA4")
ferroptosis_genes <- intersect(ferroptosis_genes, rownames(exprs_mat))
cat(sprintf("Ferroptosis panel: %d genes\n", length(ferroptosis_genes)))

# ── Observed CCI (with original weights) ──
cat("\n=== Computing observed CCI ===\n")
obs_prop <- estimate_proportions_cade(exprs_mat, marker_list, max_iter=30, tol=1e-5, verbose=FALSE)
cat(sprintf("  Converged at iteration %d (delta=%.2e)\n", obs_prop$n_iter, tail(obs_prop$convergence, 1)))
obs_de <- cade_de_analysis(exprs_mat, true_group, obs_prop$proportions, top_cts=4,
                            cci_variant=mode, verbose=FALSE)
# In stabilised mode, the CCI column is "CCI_stabilized"; in legacy mode it's "CCI_legacy"
# Use the appropriate column based on the mode
cci_col <- if(mode == "stabilized") "CCI_stabilized" else "CCI_legacy"
if(!cci_col %in% colnames(obs_de)) cci_col <- "CCI_legacy"  # fallback
if(!cci_col %in% colnames(obs_de)) cci_col <- grep("^CCI", colnames(obs_de), value=TRUE)[1]
cat(sprintf("  Using CCI column: %s\n", cci_col))
obs_cci <- setNames(
  obs_de[match(ferroptosis_genes, obs_de$Gene), cci_col],
  ferroptosis_genes
)
cat(sprintf("  Observed CCI computed for %d genes\n", length(obs_cci)))
# Safe value lookup with NAs
get_val <- function(g) if(!is.na(obs_cci[g])) sprintf("%.3f", obs_cci[g]) else "NA"
cat(sprintf("  Key values: SLC7A11=%s, FTH1=%s, GCLM=%s, HMOX1=%s, STAT1=%s\n",
            get_val("SLC7A11"), get_val("FTH1"), get_val("GCLM"),
            get_val("HMOX1"), get_val("STAT1")))

# ── Full permutation test: re-estimate weights at each shuffle ──
cat(sprintf("\n=== Full permutation test: %d shuffles with weight re-estimation ===\n", n_perm))
cat("This is SLOW. With n=%d and ~20-30 QP iterations per shuffle,\n", n_perm)
cat("expect ~%.0f minutes on a standard laptop.\n\n", n_perm * 1.5)

set.seed(seed)
perm_cci <- matrix(NA, nrow=length(ferroptosis_genes), ncol=n_perm,
                   dimnames=list(ferroptosis_genes, NULL))
perm_n_iter <- integer(n_perm)
perm_converged <- logical(n_perm)

time_start <- Sys.time()
pb <- txtProgressBar(min=0, max=n_perm, style=3)
for(p in 1:n_perm) {
  perm_group <- sample(true_group)
  perm_prop <- tryCatch({
    estimate_proportions_cade(exprs_mat, marker_list, max_iter=20, tol=1e-4, verbose=FALSE)
  }, error=function(e) NULL)
  if(!is.null(perm_prop)) {
    perm_n_iter[p] <- perm_prop$n_iter
    # converged may not exist; default to TRUE if iter > 0
    if(!is.null(perm_prop$converged)) {
      perm_converged[p] <- perm_prop$converged
    } else {
      perm_converged[p] <- perm_prop$n_iter > 0
    }
    perm_de <- tryCatch({
      cade_de_analysis(exprs_mat, perm_group, perm_prop$proportions, top_cts=4,
                       cci_variant=mode, verbose=FALSE)
    }, error=function(e) NULL)
    if(!is.null(perm_de)) {
      perm_cci[, p] <- setNames(
        perm_de[match(ferroptosis_genes, perm_de$Gene), cci_col],
        ferroptosis_genes
      )
    }
  }
  setTxtProgressBar(pb, p)
}
close(pb)
time_end <- Sys.time()
elapsed_min <- as.numeric(difftime(time_end, time_start, units="mins"))
cat(sprintf("\nCompleted in %.1f minutes\n", elapsed_min))
cat(sprintf("  Converged: %d / %d permutations (%.1f%%)\n",
            sum(perm_converged), n_perm, 100*mean(perm_converged)))
if(sum(perm_n_iter > 0) > 0) {
  cat(sprintf("  Mean iterations: %.1f (range %d-%d)\n",
              mean(perm_n_iter[perm_n_iter>0]),
              min(perm_n_iter[perm_n_iter>0]), max(perm_n_iter)))
}

# ── Compute null statistics and compare to shared-weight null ──
cat("\n=== Computing null statistics ===\n")
null_stats <- data.frame(
  Gene = ferroptosis_genes,
  CCI_observed = as.numeric(obs_cci),
  FullPerm_null_mean = rowMeans(perm_cci, na.rm=TRUE),
  FullPerm_null_sd = apply(perm_cci, 1, sd, na.rm=TRUE),
  FullPerm_null_median = apply(perm_cci, 1, median, na.rm=TRUE),
  FullPerm_valid_n = apply(perm_cci, 1, function(x) sum(!is.na(x))),
  stringsAsFactors = FALSE
)

# Compare to shared-weight null
shared_perm_file <- file.path(out_dir, "Table_S18_CCI_Permutation_Test.csv")
if(file.exists(shared_perm_file)) {
  shared <- read.csv(shared_perm_file, check.names=FALSE)
  shared_lookup <- setNames(shared$CCI_null_mean, shared$Gene)
  null_stats$SharedWeight_null_mean <- shared_lookup[null_stats$Gene]
  null_stats$Delta_null_mean <- null_stats$FullPerm_null_mean - null_stats$SharedWeight_null_mean
  valid_d <- !is.na(null_stats$Delta_null_mean) & !is.na(null_stats$SharedWeight_null_mean)
  if(any(valid_d)) {
    null_stats$Delta_null_pct[valid_d] <- 100 * abs(null_stats$Delta_null_mean[valid_d]) /
                                              pmax(0.01, abs(null_stats$SharedWeight_null_mean[valid_d]))
  }
  valid <- !is.na(null_stats$FullPerm_null_mean) & !is.na(null_stats$SharedWeight_null_mean)
  if(sum(valid) > 3) {
    rho <- suppressWarnings(cor(null_stats$FullPerm_null_mean[valid],
                                  null_stats$SharedWeight_null_mean[valid],
                                  method="spearman"))
    cat(sprintf("  Spearman ρ (full vs shared null mean): %.3f\n", rho))
  }
}

# Save outputs
out_file <- file.path(out_dir, "Table_S_FullPermutation_vs_Shared.csv")
write.csv(null_stats, out_file, row.names=FALSE)
cat(sprintf("\n  Saved: %s\n", out_file))

raw_file <- file.path(out_dir, "Table_S_FullPermutation_RawPermCCI.csv")
write.csv(perm_cci, raw_file)
cat(sprintf("  Saved: %s\n", raw_file))

# ── Summary ──
cat("\n=== Summary ===\n")
cat(sprintf("Mean full-permutation null mean: %.3f (SD across genes: %.3f)\n",
            mean(null_stats$FullPerm_null_mean), sd(null_stats$FullPerm_null_mean)))
if("SharedWeight_null_mean" %in% names(null_stats)) {
  cat(sprintf("Mean shared-weight null mean: %.3f (SD across genes: %.3f)\n",
              mean(null_stats$SharedWeight_null_mean, na.rm=TRUE),
              sd(null_stats$SharedWeight_null_mean, na.rm=TRUE)))
  cat(sprintf("Mean |Δ| null mean: %.3f (range %.3f - %.3f)\n",
              mean(abs(null_stats$Delta_null_mean), na.rm=TRUE),
              min(abs(null_stats$Delta_null_mean), na.rm=TRUE),
              max(abs(null_stats$Delta_null_mean), na.rm=TRUE)))
}

if(file.exists(shared_perm_file)) {
  shared <- read.csv(shared_perm_file, check.names=FALSE)
  shared_lookup <- setNames(shared$CCI_null_mean, shared$Gene)
  # Build comparison data frame manually (avoid merge column-name suffix issues)
  comp <- data.frame(
    Gene = null_stats$Gene,
    CCI_observed = null_stats$CCI_observed,
    FullPerm_mean = null_stats$FullPerm_null_mean,
    Shared_mean = shared_lookup[null_stats$Gene],
    stringsAsFactors = FALSE
  )
  comp$below_full <- comp$CCI_observed < comp$FullPerm_mean
  comp$below_shared <- comp$CCI_observed < comp$Shared_mean
  comp$direction_agrees <- comp$below_full == comp$below_shared
  comp <- comp[!is.na(comp$Shared_mean), ]
  cat(sprintf("\nDirection agreement (full vs shared): %d / %d genes (%.1f%%)\n",
              sum(comp$direction_agrees, na.rm=TRUE),
              nrow(comp), 100*mean(comp$direction_agrees, na.rm=TRUE)))
  # Two-sided Wilcoxon signed-rank test
  if(nrow(comp) > 3) {
    diffs <- comp$FullPerm_mean - comp$Shared_mean
    wt <- suppressWarnings(wilcox.test(diffs, mu=0, exact=FALSE))
    cat(sprintf("Wilcoxon signed-rank on Δ: p=%.4f (median Δ=%.3f)\n",
                wt$p.value, median(diffs)))
  }
}

cat("\n=== Full permutation test complete ===\n")
