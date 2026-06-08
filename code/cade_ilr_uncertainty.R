#!/usr/bin/env Rscript
# ============================================================
# CADE-ILR + uncertainty-calibrated CCI
# Adds compositional-geometry covariates and bootstrap rank stability
# without replacing the legacy CADE outputs.
# ============================================================

get_script_dir_local <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args_full, value = TRUE)
  if(length(hit) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  }
  try_frame <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if(!is.null(try_frame) && nzchar(try_frame)) {
    return(dirname(normalizePath(try_frame)))
  }
  getwd()
}

get_arg_value <- function(args, flag, default) {
  idx <- which(args == flag)
  if(length(idx) == 0 || idx[1] >= length(args)) {
    return(default)
  }
  args[idx[1] + 1]
}

SCRIPT_DIR <- get_script_dir_local()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork=TRUE)
CODE_DIR <- file.path(PROJECT_ROOT, "code")

args <- commandArgs(trailingOnly=TRUE)
n_bootstrap <- as.integer(get_arg_value(args, "--n-bootstrap", "200"))
top_cts <- as.integer(get_arg_value(args, "--top-cts", "4"))
bootstrap_seed <- as.integer(get_arg_value(args, "--seed", "42"))
out_dir_arg <- get_arg_value(args, "--out-dir", file.path(PROJECT_ROOT, "analysis_output", "CADE"))
ILR_OUT_DIR <- if(grepl("^/", out_dir_arg)) out_dir_arg else file.path(PROJECT_ROOT, out_dir_arg)
dir.create(ILR_OUT_DIR, showWarnings=FALSE, recursive=TRUE)

if(is.na(n_bootstrap) || n_bootstrap < 5) {
  stop("--n-bootstrap must be an integer >= 5.")
}
if(is.na(top_cts) || top_cts < 1) {
  stop("--top-cts must be an integer >= 1.")
}
if(is.na(bootstrap_seed)) {
  stop("--seed must be an integer.")
}

old_skip_main <- Sys.getenv("CADE_SKIP_MAIN", unset=NA)
Sys.setenv(CADE_SKIP_MAIN="1")
source(file.path(CODE_DIR, "cade_method.R"), local=FALSE)
if(is.na(old_skip_main)) {
  Sys.unsetenv("CADE_SKIP_MAIN")
} else {
  Sys.setenv(CADE_SKIP_MAIN=old_skip_main)
}
set.seed(bootstrap_seed)

cat("=== CADE-ILR + uncertainty-calibrated CCI ===\n")
cat(sprintf("Project root: %s\n", PROJECT_ROOT))
cat(sprintf("Bootstrap iterations: %d\n", n_bootstrap))
cat(sprintf("Top composition covariates: %d\n", top_cts))
cat(sprintf("Bootstrap seed: %d\n", bootstrap_seed))

expr_path <- file.path(PROJECT_ROOT, "geo_analysis_output", "GSE26050_expression_corrected.csv")
if(!file.exists(expr_path)) {
  stop(sprintf("Corrected expression matrix not found: %s", expr_path))
}
exprs_mat <- as.matrix(read.csv(expr_path, row.names=1, check.names=FALSE))
group <- get_fhl_group(colnames(exprs_mat))
marker_list <- get_fhl_marker_list()

ferroptosis_genes <- c("SLC7A11", "SLC25A37", "FTH1", "FTL", "GPX4",
                       "GCLM", "TFRC", "HMOX1", "SLC40A1", "NCOA4",
                       "STAT3", "JAK2", "IFNG", "STAT1", "NFE2L2",
                       "IL1B", "TNF", "IL6", "CXCL8", "NFKB1")
ferroptosis_genes <- intersect(ferroptosis_genes, rownames(exprs_mat))

cat("Estimating marker-derived weights...\n")
prop_cade <- estimate_proportions_cade(exprs_mat, marker_list, max_iter=50, tol=1e-6)

cat("Running raw-weight stabilized CADE for comparison...\n")
de_raw <- cade_de_analysis(
  exprs_mat, group, prop_cade$proportions,
  top_cts=top_cts,
  cci_variant="stabilized",
  composition_transform="raw"
)

cat("Running ILR-weight stabilized CADE...\n")
de_ilr <- cade_de_analysis(
  exprs_mat, group, prop_cade$proportions,
  top_cts=top_cts,
  cci_variant="stabilized",
  composition_transform="ilr"
)

raw_info <- attr(de_raw, "cade_covariates")
ilr_info <- attr(de_ilr, "cade_covariates")

covariate_table <- rbind(
  data.frame(
    Transform="raw",
    Covariate=names(raw_info$variance),
    Variance=as.numeric(raw_info$variance),
    Selected=names(raw_info$variance) %in% raw_info$selected,
    stringsAsFactors=FALSE
  ),
  data.frame(
    Transform="ilr",
    Covariate=names(ilr_info$variance),
    Variance=as.numeric(ilr_info$variance),
    Selected=names(ilr_info$variance) %in% ilr_info$selected,
    stringsAsFactors=FALSE
  )
)
covariate_table <- covariate_table[order(covariate_table$Transform, -covariate_table$Variance), ]

cat("Running ILR marker-dropout bootstrap and rank-stability analysis...\n")
boot_ilr <- cade_bootstrap(
  exprs_mat, marker_list, group,
  n_bootstrap=n_bootstrap,
  drop_fraction=0.2,
  top_cts=top_cts,
  cci_variant="stabilized",
  composition_transform="ilr",
  key_genes=ferroptosis_genes
)

de_raw_small <- de_raw[, c("Gene", "logFC.adj", "Delta_logFC", "CCI",
                           "CCI_lower_approx", "CCI_upper_approx",
                           "CompositionCovariates")]
names(de_raw_small) <- c("Gene", "logFC_adj_raw", "Delta_logFC_raw",
                         "CCI_raw_stabilized", "CCI_raw_lower_approx",
                         "CCI_raw_upper_approx", "RawCovariates")

de_ilr_small <- de_ilr[, c("Gene", "logFC.adj", "Delta_logFC", "CCI",
                           "CCI_lower_approx", "CCI_upper_approx",
                           "CCI_interval_method", "CompositionCovariates")]
names(de_ilr_small) <- c("Gene", "logFC_adj_ilr", "Delta_logFC_ilr",
                         "CCI_ilr_stabilized", "CCI_ilr_lower_approx",
                         "CCI_ilr_upper_approx", "CCI_ilr_interval_method",
                         "ILRCovariates")

comparison <- merge(de_raw_small, de_ilr_small, by="Gene")
comparison$CCI_ilr_minus_raw <- comparison$CCI_ilr_stabilized - comparison$CCI_raw_stabilized
comparison$Abs_CCI_ilr_minus_raw <- abs(comparison$CCI_ilr_minus_raw)
comparison$logFC_adj_ilr_minus_raw <- comparison$logFC_adj_ilr - comparison$logFC_adj_raw

de_ilr_ferr <- de_ilr[de_ilr$Gene %in% ferroptosis_genes, ]
de_ilr_ferr <- de_ilr_ferr[order(de_ilr_ferr$CCI, de_ilr_ferr$adj.P.Val.adj), ]
boot_summary <- boot_ilr$summary
names(boot_summary)[names(boot_summary) == "CCI_median"] <- "Boot_CCI_median"
names(boot_summary)[names(boot_summary) == "CCI_lower"] <- "Boot_CCI_lower"
names(boot_summary)[names(boot_summary) == "CCI_upper"] <- "Boot_CCI_upper"
names(boot_summary)[names(boot_summary) == "CCI_rank_median"] <- "Boot_CCI_rank_median"
names(boot_summary)[names(boot_summary) == "CCI_rank_lower"] <- "Boot_CCI_rank_lower"
names(boot_summary)[names(boot_summary) == "CCI_rank_upper"] <- "Boot_CCI_rank_upper"

boot_panel_cols <- c("Gene", "logFC_adj_median", "logFC_adj_lower", "logFC_adj_upper",
                     "Boot_CCI_median", "Boot_CCI_lower", "Boot_CCI_upper",
                     "Boot_CCI_rank_median", "Boot_CCI_rank_lower", "Boot_CCI_rank_upper",
                     "Prob_CCI_lt_0_2", "Prob_CCI_gt_0_5", "Prob_top5_low_CCI",
                     "Bootstrap_Evaluable")
de_ilr_ferr <- merge(de_ilr_ferr, boot_summary[, boot_panel_cols], by="Gene", all.x=TRUE)
de_ilr_ferr <- de_ilr_ferr[order(de_ilr_ferr$CCI, de_ilr_ferr$adj.P.Val.adj), ]
de_ilr_ferr$CCI_Rank_Tier <- with(de_ilr_ferr, ifelse(
  is.na(CCI), "Undefined",
  ifelse(CCI < 0.2, "Lowest",
         ifelse(CCI <= 0.5, "Low-moderate", "High"))
))

write.csv(de_ilr, file.path(ILR_OUT_DIR, "Table_CADE_ILR_Full_DE_Results.csv"), row.names=FALSE)
write.csv(de_ilr_ferr, file.path(ILR_OUT_DIR, "Table_CADE_ILR_Ferroptosis_Genes.csv"), row.names=FALSE)
write.csv(boot_summary, file.path(ILR_OUT_DIR, "Table_CADE_ILR_Bootstrap_RankStability.csv"), row.names=FALSE)
write.csv(comparison, file.path(ILR_OUT_DIR, "Table_CADE_ILR_vs_Raw_Comparison.csv"), row.names=FALSE)
write.csv(covariate_table, file.path(ILR_OUT_DIR, "Table_CADE_ILR_Covariates.csv"), row.names=FALSE)

if(requireNamespace("ggplot2", quietly=TRUE)) {
  plot_df <- de_ilr_ferr[!is.na(de_ilr_ferr$Boot_CCI_median), ]
  if(nrow(plot_df) > 0) {
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x=reorder(Gene, Boot_CCI_median), y=Boot_CCI_median)
    ) +
      ggplot2::geom_pointrange(
        ggplot2::aes(ymin=Boot_CCI_lower, ymax=Boot_CCI_upper,
                     color=Prob_top5_low_CCI),
        size=0.8
      ) +
      ggplot2::coord_flip(ylim=c(0, 1)) +
      ggplot2::scale_color_gradient(low="#d73027", high="#1a9850", limits=c(0, 1)) +
      ggplot2::theme_bw(base_size=10) +
      ggplot2::labs(
        title="CADE-ILR bootstrap CCI and rank stability",
        subtitle="Point/range: marker-dropout bootstrap CCI; color: probability of top-5 lowest CCI",
        x="", y="ILR-stabilized CCI", color="Pr(top-5 stable)"
      )
    ggplot2::ggsave(file.path(ILR_OUT_DIR, "Figure_CADE_ILR_RankStability.png"),
                    p, width=8, height=6, dpi=300)
    ggplot2::ggsave(file.path(ILR_OUT_DIR, "Figure_CADE_ILR_RankStability.tif"),
                    p, width=8, height=6, dpi=300)
  }
}

cat("\n=== CADE-ILR summary ===\n")
cat(sprintf("Selected raw covariates: %s\n", paste(raw_info$selected, collapse=", ")))
cat(sprintf("Selected ILR covariates: %s\n", paste(ilr_info$selected, collapse=", ")))
cat(sprintf("Median absolute CCI shift, ILR vs raw: %.3f\n",
            median(comparison$Abs_CCI_ilr_minus_raw, na.rm=TRUE)))
cat("Top low-CCI ILR panel genes:\n")
print(head(de_ilr_ferr[, c("Gene", "logFC.unadj", "logFC.adj", "CCI",
                           "CCI_lower_approx", "CCI_upper_approx",
                           "Boot_CCI_median", "Boot_CCI_lower", "Boot_CCI_upper",
                           "Prob_top5_low_CCI")], 8), row.names=FALSE)
cat(sprintf("\nOutputs saved to: %s\n", ILR_OUT_DIR))
