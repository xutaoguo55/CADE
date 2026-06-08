# C2: Treatment sensitivity simulation
# Sensitivity analysis: simulate hypothetical dexamethasone/etoposide exposure effects
# on FHL PBMC DE estimates if treatment-naïve status were misclassified
# Scenarios: 25% (3/11), 50% (6/11), 100% (11/11) treated
set.seed(42)

suppressPackageStartupMessages({
  library(limma)
})

PR <- "/Users/guoxutao/.openclaw/workspace/HLH_Research/CADE_Submission_Package"
SRC <- file.path(PR, "submission_upload_nargab_2026-06-04")

exprs <- as.matrix(read.csv(
  file.path(PR, "geo_analysis_output/GSE26050_expression_corrected.csv"),
  row.names = 1, check.names = FALSE
))

group <- factor(c(rep("FHL", 11), rep("HC", 33)))
design <- model.matrix(~ group)

# --- Dexamethasone-responsive NF-kB targets ---
# Literature: dex 1 mg/kg downregulates ~30% of NF-kB targets within 24h
# Effect sizes from published PBMC/monocyte dex stimulation studies
dex_targets <- list(
  IL1B  = c(-1.2, -1.8),  # strong suppression
  CXCL8 = c(-1.0, -1.6),
  TNF   = c(-0.8, -1.4),
  IL6   = c(-1.0, -1.5),
  CCL2  = c(-0.7, -1.2),
  CXCL10 = c(-0.6, -1.0),
  NFKB1 = c(-0.3, -0.6),  # moderate
  RELA  = c(-0.2, -0.5),
  PTGS2 = c(-0.8, -1.3),
  ICAM1 = c(-0.5, -0.9),
  VCAM1 = c(-0.4, -0.8),
  SELE  = c(-0.5, -0.9),
  MMP9  = c(-0.6, -1.0),
  CCL5  = c(-0.5, -0.9),
  CXCL9 = c(-0.4, -0.8),
  IL18  = c(-0.3, -0.6),
  IL1A  = c(-0.6, -1.0),
  CASP1 = c(-0.2, -0.4),
  NLRP3 = c(-0.3, -0.6)
)

# --- Etoposide-responsive genes ---
# Etoposide induces DNA damage response, apoptosis, cell-cycle arrest
etoposide_targets <- list(
  BAX   = c(0.3, 0.7),
  BBC3  = c(0.4, 0.9),   # PUMA
  PMAIP1 = c(0.3, 0.8),  # NOXA
  CDKN1A = c(0.5, 1.2),  # p21
  TP53  = c(0.2, 0.5),
  GADD45A = c(0.3, 0.7),
  FAS   = c(0.2, 0.5),
  TNFRSF10B = c(0.3, 0.7),  # DR5
  BCL2  = c(-0.3, -0.6),
  MYC   = c(-0.4, -0.8),
  TOP2A = c(-0.5, -1.0),
  CCND1 = c(-0.3, -0.7),
  CCNB1 = c(-0.4, -0.8)
)

all_targets <- c(names(dex_targets), names(etoposide_targets))
common <- intersect(all_targets, rownames(exprs))
cat(sprintf("Treatment targets in expression matrix: %d/%d\n",
  length(common), length(all_targets)))

# Missing genes
missing <- setdiff(all_targets, rownames(exprs))
if (length(missing) > 0) {
  cat("Missing:", paste(missing, collapse = ", "), "\n")
}

# --- Simulation function ---
run_scenario <- function(n_treated, effect_type = "moderate", n_reps = 100) {
  # effect_type: "moderate" = midpoint, "extreme" = worst-case (max suppression/induction)
  results <- matrix(NA, nrow = n_reps, ncol = length(common))
  colnames(results) <- common

  for (rep_i in seq_len(n_reps)) {
    exprs_mod <- exprs
    treated_idx <- sample(1:11, n_treated)

    for (gene in common) {
      vals <- exprs[gene, ]
      if (gene %in% names(dex_targets)) {
        eff_range <- dex_targets[[gene]]
      } else {
        eff_range <- etoposide_targets[[gene]]
      }

      if (effect_type == "extreme") {
        eff <- eff_range[2]  # max effect
      } else {
        eff <- mean(eff_range)
      }

      # Add treatment noise (±20% of effect size)
      treatment_effects <- rnorm(n_treated, mean = eff, sd = abs(eff) * 0.2)
      vals[treated_idx] <- vals[treated_idx] + treatment_effects
      exprs_mod[gene, ] <- vals
    }

    fit <- lmFit(exprs_mod, design)
    fit <- eBayes(fit)
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    # coef 2 is "groupHC"; negate to get FHL - HC
    results[rep_i, ] <- -tt[common, "logFC"]
  }
  results
}

# --- Run scenarios ---
scenarios <- list(
  "3/11 (25%%)"  = 3,
  "6/11 (55%%)"  = 6,
  "11/11 (100%%)" = 11
)

cat("\n=== Conservative treatment-confound simulation ===\n")
cat(sprintf("Replicates per scenario: %d\n", 100))

all_outputs <- list()

for (sc_name in names(scenarios)) {
  n_treated <- scenarios[[sc_name]]
  cat(sprintf("\n--- %s treated ---\n", sc_name))

  # Moderate effect
  mod <- run_scenario(n_treated, "moderate")
  # Extreme effect
  ext <- run_scenario(n_treated, "extreme")

  # Baseline: no treatment
  base <- run_scenario(0, "moderate")

  # Compute shifts
  delta_mod <- colMeans(mod) - colMeans(base)
  delta_ext <- colMeans(ext) - colMeans(base)

  out <- data.frame(
    gene = common,
    delta_moderate = round(delta_mod, 4),
    delta_extreme = round(delta_ext, 4),
    sd_moderate = round(apply(mod, 2, sd), 4),
    sd_extreme = round(apply(ext, 2, sd), 4),
    row.names = NULL
  )
  out <- out[order(out$delta_extreme), ]  # sort by extreme effect

  all_outputs[[sc_name]] <- out

  cat(sprintf("  Mean |delta| (moderate): %.3f\n", mean(abs(delta_mod))))
  cat(sprintf("  Mean |delta| (extreme): %.3f\n", mean(abs(delta_ext))))

  # Top 5 most affected
  cat("  Top 5 extreme-affected:\n")
  print(head(out, 5)[, c("gene", "delta_moderate", "delta_extreme")])

  # Highlight IL1B, CXCL8, TNF
  for (g in c("IL1B", "CXCL8", "TNF")) {
    if (g %in% common) {
      row <- out[out$gene == g, ]
      cat(sprintf("  %s: moderate Δ=%.3f, extreme Δ=%.3f\n",
        g, row$delta_moderate, row$delta_extreme))
    }
  }
}

# --- Save outputs ---
out_dir <- file.path(SRC, "analysis_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Summary table for manuscript
summary_tab <- data.frame(
  Scenario = names(scenarios),
  n_treated = unlist(scenarios),
  mean_abs_delta_moderate = sapply(all_outputs, function(x) mean(abs(x$delta_moderate))),
  mean_abs_delta_extreme = sapply(all_outputs, function(x) mean(abs(x$delta_extreme))),
  max_delta_moderate = sapply(all_outputs, function(x) min(x$delta_moderate)),
  max_delta_extreme = sapply(all_outputs, function(x) min(x$delta_extreme)),
  IL1B_delta_moderate = sapply(all_outputs, function(x) x$delta_moderate[x$gene == "IL1B"]),
  IL1B_delta_extreme = sapply(all_outputs, function(x) x$delta_extreme[x$gene == "IL1B"]),
  CXCL8_delta_moderate = sapply(all_outputs, function(x) x$delta_moderate[x$gene == "CXCL8"]),
  CXCL8_delta_extreme = sapply(all_outputs, function(x) x$delta_extreme[x$gene == "CXCL8"]),
  TNF_delta_moderate = sapply(all_outputs, function(x) x$delta_moderate[x$gene == "TNF"]),
  TNF_delta_extreme = sapply(all_outputs, function(x) x$delta_extreme[x$gene == "TNF"]),
  row.names = NULL
)

write.csv(summary_tab,
  file.path(out_dir, "C2_treatment_confound_summary.csv"),
  row.names = FALSE)
cat(sprintf("\nSaved summary to C2_treatment_confound_summary.csv\n"))

# Full per-gene tables per scenario
for (sc_name in names(all_outputs)) {
  fname <- gsub("[ /()%]+", "_", sc_name)
  fname <- gsub("_+$", "", fname)
  write.csv(all_outputs[[sc_name]],
    file.path(out_dir, sprintf("C2_treatment_confound_%s.csv", fname)),
    row.names = FALSE)
}
cat("Saved per-scenario tables.\n")
cat("\n=== DONE ===\n")
