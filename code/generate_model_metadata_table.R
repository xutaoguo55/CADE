#!/usr/bin/env Rscript
# Generate reviewer-facing model formula, parameter, metadata and collinearity table.

script_path <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(script_path)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_path <- if (length(file_arg)) normalizePath(sub("^--file=", "", file_arg[1])) else normalizePath("code/generate_model_metadata_table.R")
}
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

raw_dir <- file.path(project_root, "supplementary", "raw_csv_components")
table_dir <- file.path(project_root, "tables")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

weight_file <- file.path(raw_dir, "Table_S24_FHL_CADE_Weights.csv")
if (!file.exists(weight_file)) {
  stop(sprintf("Missing FHL CADE weights: %s", weight_file))
}

weights <- read.csv(weight_file, row.names = 1, check.names = FALSE)
weights <- as.matrix(weights)
storage.mode(weights) <- "numeric"
weights[!is.finite(weights)] <- 0

raw_selected <- c("CD4_Tcells", "Erythrocytes", "Neutrophils", "Macrophages")
missing_raw <- setdiff(raw_selected, colnames(weights))
if (length(missing_raw)) {
  stop(sprintf("Missing selected raw covariates: %s", paste(missing_raw, collapse = ", ")))
}

compute_vif <- function(mat) {
  out <- lapply(seq_len(ncol(mat)), function(j) {
    response <- mat[, j]
    predictors <- mat[, -j, drop = FALSE]
    fit <- lm(response ~ predictors)
    r2 <- summary(fit)$r.squared
    data.frame(
      covariate = colnames(mat)[j],
      r_squared = r2,
      vif = 1 / (1 - r2),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

raw_mat <- weights[, raw_selected, drop = FALSE]
raw_vif <- compute_vif(raw_mat)
raw_condition <- kappa(scale(raw_mat, center = TRUE, scale = TRUE))

close_composition <- function(x, pseudocount = 1e-6) {
  x <- t(x)
  x[x < 0] <- 0
  x <- x + pseudocount
  sweep(x, 2, colSums(x), "/")
}

closed <- close_composition(weights, pseudocount = 1e-6)
basis <- contr.helmert(nrow(closed))
basis <- sweep(basis, 2, sqrt(colSums(basis^2)), "/")
rownames(basis) <- rownames(closed)
colnames(basis) <- paste0("ILR", seq_len(ncol(basis)))
ilr <- t(log(t(closed)) %*% basis)
rownames(ilr) <- colnames(basis)
colnames(ilr) <- colnames(closed)

ilr_variance <- sort(apply(ilr, 1, var), decreasing = TRUE)
ilr_selected <- names(ilr_variance)[seq_len(4)]
ilr_mat <- t(ilr[ilr_selected, , drop = FALSE])
ilr_vif <- compute_vif(ilr_mat)
ilr_condition <- kappa(scale(ilr_mat, center = TRUE, scale = TRUE))

fmt_num <- function(x, digits = 3) {
  ifelse(
    abs(x) >= 10000,
    formatC(x, digits = digits, format = "e"),
    formatC(x, digits = digits, format = "fg", flag = "#")
  )
}

vif_summary <- paste0(
  "raw selected weights max VIF=", fmt_num(max(raw_vif$vif), 4),
  ", condition number=", fmt_num(raw_condition, 4),
  "; selected ILR coordinates max VIF=", fmt_num(max(ilr_vif$vif), 4),
  ", condition number=", fmt_num(ilr_condition, 4)
)

rows <- rbind(
  data.frame(Category = "Cohort metadata", Item = "GSE26050 sample source",
             Value = "PBMC; Affymetrix GPL570; 11 untreated FHL vs 30 healthy paediatric controls",
             Reviewer_check = "GEO accession and Methods dataset overview", stringsAsFactors = FALSE),
  data.frame(Category = "Cohort metadata", Item = "Unavailable metadata",
             Value = "Individual-level genotype, treatment documentation, age-by-sample and batch labels unavailable in public GEO files; treatment status relies on deposited untreated metadata",
             Reviewer_check = "Limitations and sensitivity simulation", stringsAsFactors = FALSE),
  data.frame(Category = "Model formula", Item = "Unadjusted limma",
             Value = "expr_g ~ group",
             Reviewer_check = "group coefficient is the unadjusted bulk DE coefficient", stringsAsFactors = FALSE),
  data.frame(Category = "Model formula", Item = "Raw-weight adjusted limma",
             Value = "expr_g ~ group + CD4_Tcells + Erythrocytes + Neutrophils + Macrophages",
             Reviewer_check = "lineage-level coefficient-sensitivity model; raw weights are closed relative quantities", stringsAsFactors = FALSE),
  data.frame(Category = "Model formula", Item = "ILR adjusted limma",
             Value = sprintf("expr_g ~ group + %s", paste(ilr_selected, collapse = " + ")),
             Reviewer_check = "orthonormal-balance sensitivity model used to audit closed-composition geometry", stringsAsFactors = FALSE),
  data.frame(Category = "Primary parameters", Item = "Marker sets",
             Value = "8 lineages; 66 marker entries; 58 unique genes; 7-9 markers per lineage",
             Reviewer_check = "Supplementary Table S2", stringsAsFactors = FALSE),
  data.frame(Category = "Primary parameters", Item = "FHL covariate count",
             Value = "top_cts=4 for the FHL worked example; generic benchmark screens use top_cts=3 unless stated otherwise",
             Reviewer_check = "k-sensitivity scan in Supplementary Figure S6", stringsAsFactors = FALSE),
  data.frame(Category = "Primary parameters", Item = "CCI denominator threshold",
             Value = "primary CCI requires |logFC_unadj| > 0.1; stabilised CCI uses max(|logFC_unadj|, 0.1)",
             Reviewer_check = "threshold scan across 0.05-0.50", stringsAsFactors = FALSE),
  data.frame(Category = "Primary parameters", Item = "Bootstrap and permutation",
             Value = "200 marker-dropout bootstraps with 20% marker dropout; 1000 label permutations; seeds 42, 2026 and 31415 for ILR multiseed audit",
             Reviewer_check = "Supplementary Tables S4, S8 and S11", stringsAsFactors = FALSE),
  data.frame(Category = "Collinearity diagnostic", Item = "Raw selected weights",
             Value = paste(sprintf("%s VIF=%s", raw_vif$covariate, fmt_num(raw_vif$vif, 4)), collapse = "; "),
             Reviewer_check = sprintf("max VIF=%s; condition number=%s", fmt_num(max(raw_vif$vif), 4), fmt_num(raw_condition, 4)), stringsAsFactors = FALSE),
  data.frame(Category = "Collinearity diagnostic", Item = "Selected ILR coordinates",
             Value = paste(sprintf("%s VIF=%s", ilr_vif$covariate, fmt_num(ilr_vif$vif, 4)), collapse = "; "),
             Reviewer_check = sprintf("max VIF=%s; condition number=%s", fmt_num(max(ilr_vif$vif), 4), fmt_num(ilr_condition, 4)), stringsAsFactors = FALSE),
  data.frame(Category = "Interpretation guardrail", Item = "Panel status",
             Value = "The 20-gene ferroptosis/iron/immune panel is a targeted diagnostic panel curated before interpreting the final CCI ranks; panel results are exploratory and not confirmatory",
             Reviewer_check = "prevents post-hoc mechanism overclaiming", stringsAsFactors = FALSE),
  data.frame(Category = "Interpretation guardrail", Item = "Collinearity interpretation",
             Value = paste0(vif_summary, "; raw lineage covariates are therefore interpreted with ILR and bootstrap concordance rather than as independent causal effects"),
             Reviewer_check = "addresses small-N adjusted-model overfitting and closed-composition collinearity", stringsAsFactors = FALSE)
)

caption <- data.frame(
  Category = "Table_Caption",
  Item = "Table 1H: Reviewer-facing model metadata, parameters and collinearity diagnostics",
  Value = "",
  Reviewer_check = "",
  stringsAsFactors = FALSE
)

write.csv(rbind(caption, rows),
          file.path(table_dir, "Table_01H_Model_Formula_Parameter_Metadata.csv"),
          row.names = FALSE)
write.csv(rows,
          file.path(raw_dir, "Table_S48_ModelParameterMetadata.csv"),
          row.names = FALSE)

cat("Wrote Table_01H_Model_Formula_Parameter_Metadata.csv\n")
cat("Wrote Table_S48_ModelParameterMetadata.csv\n")
cat(sprintf("Raw selected max VIF: %.3f; raw condition number: %.3f\n", max(raw_vif$vif), raw_condition))
cat(sprintf("ILR selected max VIF: %.3f; ILR condition number: %.3f\n", max(ilr_vif$vif), ilr_condition))
