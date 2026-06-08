#!/usr/bin/env Rscript
# Quick regression checks for CADE-ILR and uncertainty-calibrated CCI.

get_script_dir_test <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args_full, value = TRUE)
  if(length(hit) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  }
  try_frame <- tryCatch(sys.frame(1)$ofile, error=function(e) NULL)
  if(!is.null(try_frame) && nzchar(try_frame)) {
    return(dirname(normalizePath(try_frame)))
  }
  getwd()
}

old_skip_main <- Sys.getenv("CADE_SKIP_MAIN", unset=NA)
Sys.setenv(CADE_SKIP_MAIN="1")
source(file.path(get_script_dir_test(), "cade_method.R"), local=FALSE)
if(is.na(old_skip_main)) {
  Sys.unsetenv("CADE_SKIP_MAIN")
} else {
  Sys.setenv(CADE_SKIP_MAIN=old_skip_main)
}

stop_if_false <- function(x, msg) {
  if(!isTRUE(x)) stop(msg, call.=FALSE)
}

set.seed(20260601)
synth <- generate_synthetic_data(n_genes=300, n_samples=24, n_ct=5, de_genes=20)
prop <- estimate_proportions_cade(
  synth$bulk_expr, synth$marker_list,
  max_iter=8, tol=1e-4, verbose=FALSE
)

ilr <- cade_ilr_transform(prop$proportions)
stop_if_false(nrow(ilr$coordinates) == nrow(prop$proportions) - 1,
              "ILR coordinates should have D-1 rows.")
stop_if_false(ncol(ilr$coordinates) == ncol(prop$proportions),
              "ILR coordinates should preserve sample count.")
stop_if_false(all(is.finite(ilr$coordinates)),
              "ILR coordinates should be finite.")
stop_if_false(max(abs(colSums(ilr$closed_composition) - 1)) < 1e-8,
              "Closed compositions should sum to 1 per sample.")

de_ilr <- cade_de_analysis(
  synth$bulk_expr, synth$group, prop$proportions,
  top_cts=4,
  cci_variant="stabilized",
  composition_transform="ilr",
  verbose=FALSE
)
required_cols <- c("CCI", "CCI_lower_approx", "CCI_upper_approx",
                   "CCI_SE_approx", "SE_adj", "CompositionTransform",
                   "CompositionCovariates")
stop_if_false(all(required_cols %in% names(de_ilr)),
              "CADE-ILR output is missing uncertainty or transform columns.")
evaluable <- !is.na(de_ilr$CCI)
stop_if_false(all(de_ilr$CCI[evaluable] >= 0 & de_ilr$CCI[evaluable] <= 1),
              "CCI values should be clipped to [0,1].")
stop_if_false(all(de_ilr$CCI_lower_approx[evaluable] <= de_ilr$CCI_upper_approx[evaluable]),
              "CCI interval lower bounds should not exceed upper bounds.")
stop_if_false(all(de_ilr$CompositionTransform == "ilr"),
              "CompositionTransform should be ilr for ILR run.")

key_genes <- intersect(synth$de_genes[1:8], rownames(synth$bulk_expr))
boot <- cade_bootstrap(
  synth$bulk_expr, synth$marker_list, synth$group,
  n_bootstrap=5,
  top_cts=4,
  cci_variant="stabilized",
  composition_transform="ilr",
  key_genes=key_genes,
  verbose=FALSE
)
boot_cols <- c("CCI_rank_median", "CCI_rank_lower", "CCI_rank_upper",
               "Prob_CCI_lt_0_2", "Prob_CCI_gt_0_5", "Prob_top5_low_CCI")
stop_if_false(all(boot_cols %in% names(boot$summary)),
              "Bootstrap summary is missing rank-stability columns.")
prob_cols <- c("Prob_CCI_lt_0_2", "Prob_CCI_gt_0_5", "Prob_top5_low_CCI")
for(col in prob_cols) {
  x <- boot$summary[[col]]
  x <- x[!is.na(x)]
  stop_if_false(all(x >= 0 & x <= 1),
                sprintf("%s should be a probability in [0,1].", col))
}

cat("CADE-ILR uncertainty regression checks passed.\n")
