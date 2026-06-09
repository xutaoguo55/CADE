#!/usr/bin/env Rscript
# ============================================================
# CADE: End-to-End Reproducibility Script
# ============================================================
# Run this script from any directory to reproduce the full
# CADE analysis pipeline. All paths are derived from the
# location of this script.
#
# Usage:
#   Rscript run_all.R
#   Rscript run_all.R --skip-benchmarks  # Skip synthetic and standalone benchmarks
#   Rscript run_all.R --step 1           # Run only step 1
#
# Output:
#   analysis_output/     — package-local CADE results, DE, GSVA, WGCNA
#   geo_analysis_output/ — package-local GEO downloads, expression matrices
# ============================================================

# ── Parse command-line arguments ───────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
skip_benchmarks <- "--skip-benchmarks" %in% args
single_step <- NULL
for (i in seq_along(args)) {
  if (args[i] == "--step" && i < length(args)) {
    single_step <- as.integer(args[i + 1])
  }
}
if ("--step" %in% args && is.null(single_step)) {
  stop("--step requires an integer value from 1 to 11.")
}
if (!is.null(single_step) && (is.na(single_step) || single_step < 1 || single_step > 12)) {
  stop("--step must be an integer from 1 to 12.")
}

# ── Path configuration ─────────────────────────────────────────
SCRIPT_DIR <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile)),
  error = function(e) {
    # Fallback for Rscript: extract --file path from commandArgs
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      script_path <- sub("^--file=", "", file_arg[1])
      return(dirname(normalizePath(script_path)))
    }
    getwd()
  }
)
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
CODE_DIR <- file.path(PROJECT_ROOT, "code")

if (!file.exists(file.path(PROJECT_ROOT, "README.md")) ||
    !dir.exists(CODE_DIR) ||
    !dir.exists(file.path(PROJECT_ROOT, "manuscript"))) {
  stop(sprintf("Could not resolve CADE package root from script directory: %s", SCRIPT_DIR))
}

# Create output directories
dir.create(file.path(PROJECT_ROOT, "analysis_output", "CADE"),
           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PROJECT_ROOT, "geo_analysis_output"),
           showWarnings = FALSE, recursive = TRUE)

cat(sprintf("CADE Reproducibility Run — %s\n", Sys.time()))
cat(sprintf("Project root: %s\n", PROJECT_ROOT))
cat(sprintf("Code dir:     %s\n", CODE_DIR))
cat(sprintf("R version:    %s\n", R.version.string))
cat(sprintf("Platform:     %s\n", R.version$platform))
cat(strrep("=", 60), "\n\n")

# ── Helper: run a script and report status ─────────────────────
run_script <- function(script_name, step_num, description) {
  script_path <- file.path(CODE_DIR, script_name)
  cat(sprintf("\n%s\n", strrep("=", 60)))
  cat(sprintf("STEP %d: %s\n", step_num, description))
  cat(sprintf("Script: %s\n", script_name))
  cat(sprintf("%s\n\n", strrep("=", 60)))

  if (!file.exists(script_path)) {
    cat(sprintf("WARNING: Script not found: %s\n", script_path))
    return(invisible(FALSE))
  }

  start_time <- Sys.time()
  old_skip_synth <- Sys.getenv("CADE_SKIP_SYNTHETIC_BENCHMARKS", unset = NA)
  if (skip_benchmarks && identical(script_name, "cade_method.R")) {
    Sys.setenv(CADE_SKIP_SYNTHETIC_BENCHMARKS = "1")
  }
  on.exit({
    if (is.na(old_skip_synth)) {
      Sys.unsetenv("CADE_SKIP_SYNTHETIC_BENCHMARKS")
    } else {
      Sys.setenv(CADE_SKIP_SYNTHETIC_BENCHMARKS = old_skip_synth)
    }
  }, add = TRUE)
  exit_code <- tryCatch({
    source(script_path, echo = FALSE, local = FALSE)
    0L
  }, error = function(e) {
    cat(sprintf("\nERROR: %s\n", e$message))
    return(1L)
  })

  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  if (exit_code == 0) {
    cat(sprintf("\n--- Step %d COMPLETE (%.1f sec) ---\n", step_num, elapsed))
  } else {
    cat(sprintf("\n--- Step %d FAILED (%.1f sec) ---\n", step_num, elapsed))
    stop(sprintf("Step %d failed. Fix errors before continuing.", step_num))
  }
  return(invisible(exit_code == 0))
}

# ── Execution plan ─────────────────────────────────────────────
steps <- list(
  # step, script, description
  list(1, "geo_de_analysis.R",            "Download GSE26050 + standard DE"),
  list(2, "fix_cade_expression.R",        "Correct expression matrix for CADE"),
  list(3, "cade_method.R",                "CADE: weight estimation + adjusted DE + CCI"),
  list(4, "cade_ilr_uncertainty.R",       "CADE-ILR compositional covariates + CCI uncertainty/rank stability"),
  list(5, "cci_permutation_test.R",       "Permutation calibration of CCI"),
  list(6, "wgcna_standalone.R",           "WGCNA co-expression network"),
  list(7, "deeper_analysis_v3.R",         "GSVA, immune scoring, targeted enrichment with WGCNA integration"),
  list(8, "additional_analysis.R",        "TF scoring, cross-disease comparison"),
  list(9, "cade_real_scRNA_benchmark_v2.R","Real scRNA-seq benchmark (standalone)"),
  list(10, "cade_pbmc_benchmark.R",       "PBMC benchmark (standalone)"),
  list(11, "cade_external_validation_v4.R","Sepsis external validation (standalone)"),
  list(12, "empirical_comparator_runtime_benchmark.R","Empirical comparator + runtime/scalability benchmark")
)

# ── Run steps ──────────────────────────────────────────────────
for (s in steps) {
  if (!is.null(single_step) && s[[1]] != single_step) {
    cat(sprintf("Skipping step %d (--step %d)\n", s[[1]], single_step))
    next
  }
  if (skip_benchmarks && s[[1]] >= 9) {
    cat(sprintf("Skipping step %d (benchmarks skipped)\n", s[[1]]))
    next
  }
  run_script(s[[2]], s[[1]], s[[3]])
}

# ── Final report ────────────────────────────────────────────────
cat(sprintf("\n%s\n", strrep("=", 60)))
cat(sprintf("CADE REPRODUCIBILITY RUN COMPLETE — %s\n", Sys.time()))
cat(sprintf("Output directories:\n"))
cat(sprintf("  %s\n", file.path(PROJECT_ROOT, "analysis_output")))
cat(sprintf("  %s\n", file.path(PROJECT_ROOT, "geo_analysis_output")))
cat(sprintf("%s\n", strrep("=", 60)))
