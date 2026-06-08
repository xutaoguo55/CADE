#!/usr/bin/env Rscript
# ============================================================
# Multi-seed CADE-ILR bootstrap stability validation
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

split_ints <- function(x) {
  vals <- trimws(unlist(strsplit(x, ",")))
  vals <- vals[nzchar(vals)]
  as.integer(vals)
}

safe_range <- function(x, fn) {
  x <- x[is.finite(x)]
  if(length(x) == 0) return(NA_real_)
  fn(x)
}

args <- commandArgs(trailingOnly=TRUE)
SCRIPT_DIR <- get_script_dir_local()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork=TRUE)
CODE_DIR <- file.path(PROJECT_ROOT, "code")
OUT_DIR <- file.path(PROJECT_ROOT, "analysis_output", "CADE")
MULTI_DIR <- file.path(OUT_DIR, "multiseed_ilr_bootstrap")
dir.create(MULTI_DIR, showWarnings=FALSE, recursive=TRUE)

n_bootstrap <- as.integer(get_arg_value(args, "--n-bootstrap", "200"))
top_cts <- as.integer(get_arg_value(args, "--top-cts", "4"))
seeds <- split_ints(get_arg_value(args, "--seeds", "42,2026,31415"))
priority_genes <- c("SLC7A11", "STAT1", "TNF", "NFKB1", "IL1B")

if(is.na(n_bootstrap) || n_bootstrap < 5) {
  stop("--n-bootstrap must be an integer >= 5.")
}
if(is.na(top_cts) || top_cts < 1) {
  stop("--top-cts must be an integer >= 1.")
}
if(length(seeds) < 2 || any(is.na(seeds))) {
  stop("--seeds must contain at least two comma-separated integer seeds.")
}

cat("=== CADE-ILR multi-seed bootstrap stability validation ===\n")
cat(sprintf("Project root: %s\n", PROJECT_ROOT))
cat(sprintf("Bootstrap iterations per seed: %d\n", n_bootstrap))
cat(sprintf("Seeds: %s\n", paste(seeds, collapse=", ")))

all_seed_rows <- list()
rscript <- file.path(R.home("bin"), "Rscript")
runner <- file.path(CODE_DIR, "cade_ilr_uncertainty.R")

for(seed in seeds) {
  seed_dir <- file.path(MULTI_DIR, paste0("seed_", seed))
  dir.create(seed_dir, showWarnings=FALSE, recursive=TRUE)

  cat(sprintf("\n--- Running seed %d ---\n", seed))
  cmd_args <- c(
    runner,
    "--n-bootstrap", as.character(n_bootstrap),
    "--top-cts", as.character(top_cts),
    "--seed", as.character(seed),
    "--out-dir", seed_dir
  )
  status <- system2(rscript, cmd_args)
  if(!identical(status, 0L)) {
    stop(sprintf("CADE-ILR bootstrap run failed for seed %d.", seed))
  }

  ferr_path <- file.path(seed_dir, "Table_CADE_ILR_Ferroptosis_Genes.csv")
  boot_path <- file.path(seed_dir, "Table_CADE_ILR_Bootstrap_RankStability.csv")
  if(!file.exists(ferr_path) || !file.exists(boot_path)) {
    stop(sprintf("Expected CADE-ILR outputs missing for seed %d.", seed))
  }

  ferr <- read.csv(ferr_path, check.names=FALSE)
  boot <- read.csv(boot_path, check.names=FALSE)
  keep_boot <- c("Gene", "Boot_CCI_median", "Boot_CCI_lower", "Boot_CCI_upper",
                 "Boot_CCI_rank_median", "Boot_CCI_rank_lower", "Boot_CCI_rank_upper",
                 "Prob_CCI_lt_0_2", "Prob_CCI_gt_0_5", "Prob_top5_low_CCI",
                 "Bootstrap_Evaluable", "CCI_variant", "CompositionTransform")
  keep_ferr <- c("Gene", "logFC.unadj", "logFC.adj", "CCI", "CCI_lower_approx",
                 "CCI_upper_approx", "adj.P.Val.adj", "CCI_Rank_Tier")
  seed_rows <- merge(
    ferr[, intersect(keep_ferr, names(ferr)), drop=FALSE],
    boot[, intersect(keep_boot, names(boot)), drop=FALSE],
    by="Gene",
    all.x=TRUE
  )
  seed_rows$Seed <- seed
  seed_rows$N_Bootstrap_Per_Seed <- n_bootstrap
  seed_rows$PriorityGene <- seed_rows$Gene %in% priority_genes
  all_seed_rows[[as.character(seed)]] <- seed_rows
}

detailed <- do.call(rbind, all_seed_rows)
rownames(detailed) <- NULL
detailed <- detailed[order(detailed$PriorityGene, detailed$Gene, detailed$Seed,
                           decreasing=FALSE), ]

summarize_gene <- function(df) {
  data.frame(
    Gene=df$Gene[1],
    Seeds=paste(sort(unique(df$Seed)), collapse=","),
    N_Seeds=length(unique(df$Seed)),
    N_Bootstrap_Per_Seed=unique(df$N_Bootstrap_Per_Seed)[1],
    Point_CCI=if("CCI" %in% names(df)) safe_range(df$CCI, median) else NA_real_,
    Boot_CCI_median_min=safe_range(df$Boot_CCI_median, min),
    Boot_CCI_median_max=safe_range(df$Boot_CCI_median, max),
    Boot_CCI_rank_median_min=safe_range(df$Boot_CCI_rank_median, min),
    Boot_CCI_rank_median_max=safe_range(df$Boot_CCI_rank_median, max),
    Prob_top5_low_CCI_min=safe_range(df$Prob_top5_low_CCI, min),
    Prob_top5_low_CCI_max=safe_range(df$Prob_top5_low_CCI, max),
    Prob_CCI_lt_0_2_min=safe_range(df$Prob_CCI_lt_0_2, min),
    Prob_CCI_lt_0_2_max=safe_range(df$Prob_CCI_lt_0_2, max),
    Prob_CCI_gt_0_5_min=safe_range(df$Prob_CCI_gt_0_5, min),
    Prob_CCI_gt_0_5_max=safe_range(df$Prob_CCI_gt_0_5, max),
    Bootstrap_Evaluable_min=safe_range(df$Bootstrap_Evaluable, min),
    Interpretation=if(all(df$Prob_top5_low_CCI >= 0.5, na.rm=TRUE)) {
      "Stable low-CCI rank across seeds"
    } else if(all(df$Prob_CCI_gt_0_5 >= 0.5, na.rm=TRUE)) {
      "Stable high-CCI across seeds"
    } else {
      "Intermediate or seed-sensitive"
    },
    PriorityGene=df$Gene[1] %in% priority_genes,
    stringsAsFactors=FALSE
  )
}

summary_rows <- do.call(rbind, lapply(split(detailed, detailed$Gene), summarize_gene))
summary_rows <- summary_rows[order(!summary_rows$PriorityGene,
                                   summary_rows$Boot_CCI_rank_median_min,
                                   summary_rows$Gene), ]

detail_path <- file.path(OUT_DIR, "Table_CADE_ILR_Multiseed_Stability.csv")
summary_path <- file.path(OUT_DIR, "Table_CADE_ILR_Multiseed_Stability_Summary.csv")
write.csv(detailed, detail_path, row.names=FALSE)
write.csv(summary_rows, summary_path, row.names=FALSE)

cat("\n=== Multi-seed priority-gene summary ===\n")
print(summary_rows[summary_rows$PriorityGene,
                   c("Gene", "Boot_CCI_median_min", "Boot_CCI_median_max",
                     "Boot_CCI_rank_median_min", "Boot_CCI_rank_median_max",
                     "Prob_top5_low_CCI_min", "Prob_top5_low_CCI_max",
                     "Prob_CCI_gt_0_5_min", "Prob_CCI_gt_0_5_max",
                     "Interpretation")],
      row.names=FALSE)
cat(sprintf("\nDetailed multi-seed output: %s\n", detail_path))
cat(sprintf("Summary multi-seed output:  %s\n", summary_path))
