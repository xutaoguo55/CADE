#!/usr/bin/env Rscript
# ============================================================
# CADE: Cell-Aware Differential Expression
# A method for composition-corrected DE analysis in small-N
# bulk transcriptomic studies
# ============================================================
library(limma)
library(quadprog)
library(dplyr)
library(ggplot2)
library(patchwork)

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

set.seed(42)
OUT_DIR <- file.path(PROJECT_ROOT, "analysis_output", "CADE")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

cat("=== CADE: Cell-Aware Differential Expression ===\n")

# ============================================================
# Part 1: Core Algorithm
# ============================================================

safe_softmax <- function(x) {
  x <- as.numeric(x)
  if(all(!is.finite(x))) {
    return(rep(1 / length(x), length(x)))
  }
  x[!is.finite(x)] <- min(x[is.finite(x)], na.rm = TRUE)
  e_x <- exp(x - max(x, na.rm = TRUE))
  s <- sum(e_x)
  if(!is.finite(s) || s <= 0) {
    return(rep(1 / length(x), length(x)))
  }
  e_x / s
}

close_composition <- function(prop_mat, pseudocount=1e-6) {
  prop_mat <- as.matrix(prop_mat)
  storage.mode(prop_mat) <- "numeric"
  if(is.null(rownames(prop_mat))) {
    rownames(prop_mat) <- paste0("Component", seq_len(nrow(prop_mat)))
  }
  if(any(!is.finite(prop_mat))) {
    prop_mat[!is.finite(prop_mat)] <- 0
  }
  prop_mat[prop_mat < 0] <- 0
  prop_mat <- prop_mat + pseudocount
  col_sums <- colSums(prop_mat)
  bad <- !is.finite(col_sums) | col_sums <= 0
  if(any(bad)) {
    prop_mat[, bad] <- 1 / nrow(prop_mat)
    col_sums <- colSums(prop_mat)
  }
  sweep(prop_mat, 2, col_sums, "/")
}

cade_ilr_transform <- function(prop_mat, pseudocount=1e-6) {
  closed <- close_composition(prop_mat, pseudocount=pseudocount)
  n_comp <- nrow(closed)
  if(n_comp < 2) {
    stop("ILR transform requires at least two composition components.")
  }
  basis <- contr.helmert(n_comp)
  basis <- sweep(basis, 2, sqrt(colSums(basis^2)), "/")
  rownames(basis) <- rownames(closed)
  colnames(basis) <- paste0("ILR", seq_len(ncol(basis)))

  coords <- t(log(t(closed)) %*% basis)
  rownames(coords) <- colnames(basis)
  colnames(coords) <- colnames(closed)
  attr(coords, "basis") <- basis
  attr(coords, "pseudocount") <- pseudocount

  list(
    coordinates = coords,
    closed_composition = closed,
    basis = basis,
    pseudocount = pseudocount
  )
}

select_cade_covariates <- function(prop_mat, top_cts=3,
                                   composition_transform=c("raw", "ilr"),
                                   ilr_pseudocount=1e-6) {
  composition_transform <- match.arg(composition_transform)
  if(composition_transform == "ilr") {
    transformed <- cade_ilr_transform(prop_mat, pseudocount=ilr_pseudocount)
    covar_mat <- transformed$coordinates
  } else {
    transformed <- NULL
    covar_mat <- as.matrix(prop_mat)
  }
  if(nrow(covar_mat) < 1) {
    stop("No composition covariates available.")
  }
  covar_vars <- apply(covar_mat, 1, var, na.rm=TRUE)
  covar_vars[!is.finite(covar_vars)] <- 0
  n_use <- min(top_cts, nrow(covar_mat))
  covar_use <- names(sort(covar_vars, decreasing=TRUE))[seq_len(n_use)]

  list(
    covariates = covar_mat,
    selected = covar_use,
    variance = covar_vars,
    transform = composition_transform,
    transformed = transformed
  )
}

safe_quantile <- function(x, prob) {
  if(all(is.na(x))) return(NA_real_)
  unname(quantile(x, prob, na.rm=TRUE, names=FALSE))
}

safe_mean <- function(x) {
  if(all(is.na(x))) return(NA_real_)
  mean(x, na.rm=TRUE)
}

#' Iterative Marker-Based Weight Estimation
estimate_proportions_cade <- function(expr_mat, marker_list, max_iter=50,
                                       tol=1e-6, min_markers=3, ridge=1e-8,
                                       verbose=TRUE) {
  n_ct <- length(marker_list)
  ct_names <- names(marker_list)
  n_samples <- ncol(expr_mat)

  # Filter to markers present in data
  marker_present <- lapply(marker_list, function(m) intersect(m, rownames(expr_mat)))
  n_markers_present <- sapply(marker_present, length)
  if(any(n_markers_present < min_markers)) {
    stop(sprintf("Cell types with insufficient markers: %s",
         paste(ct_names[n_markers_present < min_markers], collapse=", ")))
  }
  if(isTRUE(verbose)) {
    cat(sprintf("  Markers per cell type: %s\n",
        paste(sprintf("%s=%d", ct_names, n_markers_present), collapse=", ")))
  }

  # Step 1: Initial marker-derived weight estimates via mean marker expression
  init_scores <- matrix(NA, nrow=n_ct, ncol=n_samples)
  rownames(init_scores) <- ct_names
  for(ct in ct_names) {
    init_scores[ct, ] <- colMeans(expr_mat[marker_present[[ct]], , drop=FALSE])
  }
  # Softmax normalization per sample
  prop_mat <- apply(init_scores, 2, safe_softmax)
  rownames(prop_mat) <- ct_names

  # Step 2: Iterative refinement
  all_markers <- unique(unlist(marker_present))
  n_markers_total <- length(all_markers)
  marker_ct_map <- matrix(0, nrow=n_ct, ncol=n_markers_total,
                          dimnames=list(ct_names, all_markers))
  for(ct in ct_names) {
    marker_ct_map[ct, marker_present[[ct]]] <- 1
  }

  history <- list()
  for(iter in 1:max_iter) {
    # 2a: Estimate cell-type-specific expression profiles
    ct_profiles <- matrix(0, nrow=n_markers_total, ncol=n_ct,
                          dimnames=list(all_markers, ct_names))
    for(g in all_markers) {
      y <- expr_mat[g, ]
      for(ct in ct_names) {
        if(marker_ct_map[ct, g] == 1) {
          X <- prop_mat[ct, ]
          fit <- lm(y ~ X)
          coef_ct <- coef(fit)
          ct_profiles[g, ct] <- ifelse(length(coef_ct) >= 2 && is.finite(coef_ct[2]), coef_ct[2], 0)
        }
      }
    }
    ct_profiles[is.na(ct_profiles)] <- 0

    # 2b: Update marker-derived weights via constrained quadratic programming
    prop_old <- prop_mat
    for(i in 1:n_samples) {
      y_i <- expr_mat[all_markers, i]
      Dmat <- t(ct_profiles) %*% ct_profiles + diag(ridge, n_ct)
      dvec <- t(ct_profiles) %*% y_i
      Amat <- cbind(rep(1, n_ct), diag(n_ct))
      bvec <- c(1, rep(0, n_ct))

      sol <- tryCatch(
        solve.QP(Dmat, dvec, Amat, bvec, meq=1),
        error = function(e) NULL
      )
      if(!is.null(sol)) {
        w <- pmax(sol$solution, 0)
        w_sum <- sum(w)
        if(is.finite(w_sum) && w_sum > 0) {
          prop_mat[, i] <- w / w_sum
        }
      }
    }

    # Check convergence
    delta <- max(abs(prop_mat - prop_old))
    history[[iter]] <- delta
    if(delta < tol) {
      if(isTRUE(verbose)) {
        cat(sprintf("  Converged at iteration %d (delta=%.2e)\n", iter, delta))
      }
      break
    }
    if(iter == max_iter) {
      if(isTRUE(verbose)) {
        cat(sprintf("  Max iterations reached (final delta=%.2e)\n", delta))
      }
    }
  }

  list(
    proportions = prop_mat,  # legacy field name; values are relative marker-derived weights
    ct_profiles = ct_profiles,
    marker_present = marker_present,
    convergence = unlist(history),
    n_iter = iter
  )
}

#' Composition-Aware Differential Expression
#'
#' CCI is an absolute coefficient-change / model-sensitivity metric.
#' It is not a formal mediation proportion or causal attributable fraction.
#' Inspect logFC.unadj, logFC.adj, Delta_logFC, and EffectChange_Direction together.
cade_de_analysis <- function(expr_mat, group, prop_mat, top_cts=3,
                             min_abs_logfc=0.1, se_floor_mult=1.0,
                             cci_variant=c("legacy", "stabilized"),
                             composition_transform=c("raw", "ilr"),
                             ilr_pseudocount=1e-6,
                             case_level=NULL,
                             verbose=TRUE) {
  cci_variant <- match.arg(cci_variant)
  composition_transform <- match.arg(composition_transform)
  if(length(group) != ncol(expr_mat)) {
    stop("Length of group must equal number of columns in expr_mat.")
  }
  if(ncol(prop_mat) != ncol(expr_mat)) {
    stop("prop_mat must have the same number of samples (columns) as expr_mat.")
  }

  if(is.logical(group)) {
    if(length(unique(group)) != 2) {
      stop("group must be binary (two levels).")
    }
    group_binary <- as.numeric(group)
  } else if(is.character(group) || is.factor(group)) {
    group_chr <- as.character(group)
    levels_present <- unique(group_chr)
    if(length(levels_present) != 2) {
      stop("group must be binary (two levels).")
    }
    if(is.null(case_level)) {
      known_case_levels <- c("FHL", "HLH", "Disease", "Diseased", "Case", "Patient",
                             "Sepsis", "SIRS", "MAS", "SJIA_MAS", "Treatment")
      matched <- known_case_levels[tolower(known_case_levels) %in% tolower(levels_present)]
      if(length(matched) == 1) {
        case_level <- levels_present[tolower(levels_present) == tolower(matched)]
      } else {
        case_level <- levels_present[1]
        warning(sprintf(
          "case_level not supplied; using first observed level as case: %s",
          case_level
        ))
      }
    }
    if(!(case_level %in% levels_present)) {
      stop(sprintf("case_level '%s' not found in group values.", case_level))
    }
    group_binary <- as.numeric(group_chr == case_level)
  } else {
    group_binary <- as.numeric(group)
    if(length(unique(group_binary)) != 2) {
      stop("group must be binary (two unique values).")
    }
  }

  # Select top-varying raw composition weights or ILR coordinates as covariates.
  covar_info <- select_cade_covariates(
    prop_mat,
    top_cts=top_cts,
    composition_transform=composition_transform,
    ilr_pseudocount=ilr_pseudocount
  )
  covar_mat <- covar_info$covariates
  ct_use <- covar_info$selected
  if(isTRUE(verbose)) {
    cat(sprintf("  Composition transform: %s\n", composition_transform))
    cat(sprintf("  Composition covariates used: %s\n", paste(ct_use, collapse=", ")))
  }

  # Design matrix: group + top marker-derived composition covariates
  design_adj <- model.matrix(~ group_binary + t(covar_mat[ct_use, , drop=FALSE]))
  design_unadj <- model.matrix(~ group_binary)

  # Fit limma models
  fit_adj <- lmFit(expr_mat, design_adj)
  fit_adj <- eBayes(fit_adj, trend=TRUE)
  de_adj <- topTable(fit_adj, coef=2, number=Inf, adjust.method="BH")
  de_adj$Gene <- rownames(de_adj)

  fit_unadj <- lmFit(expr_mat, design_unadj)
  fit_unadj <- eBayes(fit_unadj, trend=TRUE)
  de_unadj <- topTable(fit_unadj, coef=2, number=Inf, adjust.method="BH")
  de_unadj$Gene <- rownames(de_unadj)

  de_adj_sub <- de_adj[, c("Gene", "logFC", "t", "P.Value", "adj.P.Val")]
  names(de_adj_sub)[names(de_adj_sub) == "t"] <- "t.adj"

  # Merge results
  results <- merge(
    de_unadj[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")],
    de_adj_sub,
    by="Gene", suffixes=c(".unadj", ".adj")
  )

  # Compute Composition Confounding Index (CCI)
  # CCI is an absolute coefficient-change metric. The signed fields below
  # distinguish attenuation from amplification or direction reversal.
  results$Delta_logFC <- results$logFC.adj - results$logFC.unadj
  results$EffectChange_Direction <- with(results, ifelse(
    is.na(logFC.unadj) | abs(logFC.unadj) <= min_abs_logfc, "undefined",
    ifelse(sign(logFC.unadj) != sign(logFC.adj), "direction_reversal",
      ifelse(abs(logFC.adj) < abs(logFC.unadj), "attenuation", "amplification")
    )
  ))

  # Legacy CCI (for backward compatibility with previous results/manuscript tables)
  results$CCI_legacy <- with(results,
    ifelse(abs(logFC.unadj) > min_abs_logfc,
           abs(logFC.unadj - logFC.adj) / abs(logFC.unadj),
           NA)
  )
  results$CCI_legacy <- pmax(0, pmin(1, results$CCI_legacy))

  # Stabilized CCI: denominator is max(|logFC_unadj|, adaptive floor)
  # Floor uses both fixed min_abs_logfc and limma-based SE proxy.
  results$SE_unadj <- with(results, ifelse(
    is.na(t) | abs(t) < 1e-8, NA, abs(logFC.unadj / t)
  ))
  results$SE_adj <- with(results, ifelse(
    is.na(t.adj) | abs(t.adj) < 1e-8, NA, abs(logFC.adj / t.adj)
  ))
  results$CCI_floor <- with(results, ifelse(
    is.na(SE_unadj), min_abs_logfc, pmax(min_abs_logfc, se_floor_mult * SE_unadj)
  ))
  results$CCI_denominator <- with(results, pmax(abs(logFC.unadj), CCI_floor))
  results$CCI_stabilized_raw <- with(results,
    ifelse(is.na(logFC.unadj) | is.na(logFC.adj), NA,
           abs(logFC.unadj - logFC.adj) / CCI_denominator)
  )
  results$CCI_stabilized <- pmax(0, pmin(1, results$CCI_stabilized_raw))
  results$BaselineReliability <- with(results, ifelse(
    is.na(logFC.unadj), "missing",
    ifelse(abs(logFC.unadj) < CCI_floor, "small_baseline", "stable_baseline")
  ))

  # Active CCI column used downstream; legacy can still be selected explicitly.
  results$CCI <- if(cci_variant == "legacy") results$CCI_legacy else results$CCI_stabilized

  raw_den <- pmax(abs(results$logFC.unadj), 1e-8)
  stable_ratio_se <- sqrt(
    (results$logFC.adj^2 / raw_den^4) * results$SE_unadj^2 +
      (1 / raw_den^2) * results$SE_adj^2
  )
  floor_delta_se <- with(results, sqrt(SE_unadj^2 + SE_adj^2) / CCI_denominator)
  results$CCI_SE_approx <- ifelse(
    results$BaselineReliability == "stable_baseline",
    stable_ratio_se,
    floor_delta_se
  )
  results$CCI_lower_approx <- pmax(0, results$CCI - 1.96 * results$CCI_SE_approx)
  results$CCI_upper_approx <- pmin(1, results$CCI + 1.96 * results$CCI_SE_approx)
  results$CCI_interval_method <- ifelse(
    results$BaselineReliability == "stable_baseline",
    "delta_ratio_independent_cov0",
    ifelse(results$BaselineReliability == "small_baseline",
           "stabilized_delta_numerator_only",
           "not_available")
  )
  results$CompositionTransform <- composition_transform
  results$CompositionCovariates <- paste(ct_use, collapse=";")

  attr(results, "cade_covariates") <- covar_info
  results
}

#' Bootstrap Uncertainty Quantification for CADE
cade_bootstrap <- function(expr_mat, marker_list, group, n_bootstrap=100,
                           drop_fraction=0.2,
                           top_cts=3,
                           cci_variant=c("legacy", "stabilized"),
                           composition_transform=c("raw", "ilr"),
                           ilr_pseudocount=1e-6,
                           key_genes=c("SLC7A11", "SLC25A37", "FTH1", "FTL", "GPX4",
                                       "GCLM", "TFRC", "HMOX1", "SLC40A1", "NCOA4"),
                           verbose=TRUE) {
  cci_variant <- match.arg(cci_variant)
  composition_transform <- match.arg(composition_transform)
  ct_names <- names(marker_list)
  n_ct <- length(ct_names)

  key_genes <- intersect(key_genes, rownames(expr_mat))

  empty_summary <- function() {
    data.frame(
      Gene = character(),
      logFC_adj_median = numeric(),
      logFC_adj_lower = numeric(),
      logFC_adj_upper = numeric(),
      CCI_median = numeric(),
      CCI_lower = numeric(),
      CCI_upper = numeric(),
      CCI_rank_median = numeric(),
      CCI_rank_lower = numeric(),
      CCI_rank_upper = numeric(),
      Prob_CCI_lt_0_2 = numeric(),
      Prob_CCI_gt_0_5 = numeric(),
      Prob_top5_low_CCI = numeric(),
      Bootstrap_Evaluable = integer(),
      CCI_variant = character(),
      CompositionTransform = character(),
      stringsAsFactors = FALSE
    )
  }

  # Storage
  boot_props <- array(NA, dim=c(n_ct, ncol(expr_mat), n_bootstrap),
                      dimnames=list(ct_names, colnames(expr_mat), NULL))
  boot_logFC <- matrix(NA, nrow=length(key_genes), ncol=n_bootstrap,
                       dimnames=list(key_genes, NULL))
  boot_CCI <- matrix(NA, nrow=length(key_genes), ncol=n_bootstrap,
                     dimnames=list(key_genes, NULL))
  boot_CCI_rank <- matrix(NA, nrow=length(key_genes), ncol=n_bootstrap,
                          dimnames=list(key_genes, NULL))

  if(isTRUE(verbose)) {
    cat(sprintf("  Bootstrap: %d iterations", n_bootstrap))
    pb <- txtProgressBar(min=0, max=n_bootstrap, style=3)
  }

  for(b in 1:n_bootstrap) {
    marker_sub <- lapply(marker_list, function(m) {
      m_present <- intersect(m, rownames(expr_mat))
      if(length(m_present) < 3) {
        return(m_present)
      }
      n_keep <- max(3, round(length(m_present) * (1 - drop_fraction)))
      sample(m_present, n_keep)
    })

    prop_est <- tryCatch({
      estimate_proportions_cade(expr_mat, marker_sub, max_iter=30, tol=1e-4,
                                verbose=FALSE)
    }, error = function(e) NULL)

    if(!is.null(prop_est)) {
      boot_props[, , b] <- prop_est$proportions

      de_boot <- tryCatch({
        cade_de_analysis(
          expr_mat, group, prop_est$proportions, top_cts=top_cts,
          cci_variant=cci_variant,
          composition_transform=composition_transform,
          ilr_pseudocount=ilr_pseudocount,
          verbose=FALSE
        )
      }, error = function(e) NULL)

      if(!is.null(de_boot)) {
        rownames(de_boot) <- de_boot$Gene
        boot_logFC[, b] <- de_boot[key_genes, "logFC.adj"]
        boot_CCI[, b] <- de_boot[key_genes, "CCI"]
        boot_CCI_rank[, b] <- rank(boot_CCI[, b], na.last="keep", ties.method="average")
      }
    }
    if(isTRUE(verbose)) setTxtProgressBar(pb, b)
  }
  if(isTRUE(verbose)) close(pb)

  top_k <- min(5, length(key_genes))

  if(length(key_genes) == 0) {
    warning("No requested key genes were present in the expression matrix; returning empty bootstrap summary.")
    return(list(
      prop_dist = boot_props,
      logFC_dist = boot_logFC,
      CCI_dist = boot_CCI,
      CCI_rank_dist = boot_CCI_rank,
      summary = empty_summary()
    ))
  }

  # Summarize bootstrap distributions
  boot_summary <- data.frame(
    Gene = key_genes,
    logFC_adj_median = apply(boot_logFC, 1, safe_quantile, 0.5),
    logFC_adj_lower = apply(boot_logFC, 1, safe_quantile, 0.025),
    logFC_adj_upper = apply(boot_logFC, 1, safe_quantile, 0.975),
    CCI_median = apply(boot_CCI, 1, safe_quantile, 0.5),
    CCI_lower = apply(boot_CCI, 1, safe_quantile, 0.025),
    CCI_upper = apply(boot_CCI, 1, safe_quantile, 0.975),
    CCI_rank_median = apply(boot_CCI_rank, 1, safe_quantile, 0.5),
    CCI_rank_lower = apply(boot_CCI_rank, 1, safe_quantile, 0.025),
    CCI_rank_upper = apply(boot_CCI_rank, 1, safe_quantile, 0.975),
    Prob_CCI_lt_0_2 = apply(boot_CCI, 1, function(x) safe_mean(x < 0.2)),
    Prob_CCI_gt_0_5 = apply(boot_CCI, 1, function(x) safe_mean(x > 0.5)),
    Prob_top5_low_CCI = apply(boot_CCI_rank, 1, function(x) safe_mean(x <= top_k)),
    Bootstrap_Evaluable = apply(!is.na(boot_CCI), 1, sum),
    CCI_variant = cci_variant,
    CompositionTransform = composition_transform,
    stringsAsFactors = FALSE
  )

  list(
    prop_dist = boot_props,
    logFC_dist = boot_logFC,
    CCI_dist = boot_CCI,
    CCI_rank_dist = boot_CCI_rank,
    summary = boot_summary
  )
}

# ============================================================
# Part 2: Benchmarking Framework
# ============================================================

# Helper: random Dirichlet samples
rdirichlet <- function(n, alpha) {
  x <- matrix(rgamma(n * length(alpha), shape=alpha, scale=1),
              nrow=n, ncol=length(alpha), byrow=TRUE)
  x / rowSums(x)
}

generate_synthetic_data <- function(n_genes=5000, n_samples=44,
                                     n_ct=5, ct_bias=TRUE, de_genes=200) {
  ct_names <- c("Neutrophils", "Bcells", "Tcells", "Monocytes", "Erythrocytes")
  ct_names <- ct_names[1:n_ct]

  # Cell-type-specific expression profiles
  ct_baseline <- matrix(rnorm(n_genes * n_ct, mean=0, sd=1.5), nrow=n_genes, ncol=n_ct)
  rownames(ct_baseline) <- paste0("Gene", 1:n_genes)
  colnames(ct_baseline) <- ct_names

  # Make some genes cell-type-specific (markers)
  n_markers_per_ct <- 8
  for(i in 1:n_ct) {
    start_idx <- (i-1) * n_markers_per_ct + 1
    end_idx <- i * n_markers_per_ct
    ct_baseline[start_idx:end_idx, i] <- ct_baseline[start_idx:end_idx, i] + 4
  }

  # True proportions
  n_case <- max(1, round(n_samples * 11 / 44))
  group <- c(rep(TRUE, n_case), rep(FALSE, n_samples - n_case))
  true_props <- matrix(0, nrow=n_ct, ncol=n_samples)
  rownames(true_props) <- ct_names
  colnames(true_props) <- paste0("Sample", 1:n_samples)

  if(ct_bias) {
    for(i in 1:n_samples) {
      if(group[i]) {
        true_props[, i] <- rdirichlet(1, c(0.35, 0.05, 0.05, 0.15, 0.40) * 50)
      } else {
        true_props[, i] <- rdirichlet(1, c(0.15, 0.25, 0.25, 0.20, 0.15) * 50)
      }
    }
  } else {
    base <- rep(1/n_ct, n_ct) * 50
    for(i in 1:n_samples) true_props[, i] <- rdirichlet(1, base)
  }

  # Generate bulk expression
  bulk_expr <- ct_baseline %*% true_props + matrix(rnorm(n_genes * n_samples, 0, 0.3),
                                                     nrow=n_genes, ncol=n_samples)

  # Add DE signal to non-marker genes
  de_start <- n_markers_per_ct * n_ct + 1
  de_end <- de_start + de_genes - 1
  for(g in de_start:de_end) {
    bulk_expr[g, group] <- bulk_expr[g, group] + rnorm(1, mean=2, sd=0.5)
  }

  rownames(bulk_expr) <- paste0("Gene", 1:n_genes)
  colnames(bulk_expr) <- paste0("Sample", 1:n_samples)

  # Marker gene sets
  marker_list <- list()
  for(i in 1:n_ct) {
    start_idx <- (i-1) * n_markers_per_ct + 1
    end_idx <- i * n_markers_per_ct
    marker_list[[ct_names[i]]] <- paste0("Gene", start_idx:end_idx)
  }

  # True DE status
  true_de <- rep(FALSE, n_genes)
  true_de[de_start:de_end] <- TRUE
  names(true_de) <- paste0("Gene", 1:n_genes)

  list(
    bulk_expr = bulk_expr,
    true_props = true_props,
    group = group,
    marker_list = marker_list,
    true_de = true_de,
    de_genes = paste0("Gene", de_start:de_end),
    ct_baseline = ct_baseline
  )
}

benchmark_cade <- function(synth, n_bootstrap=100) {
  cat("\n--- Benchmarking CADE ---\n")

  # 1. CADE
  cat("  Running CADE...\n")
  prop_cade <- estimate_proportions_cade(synth$bulk_expr, synth$marker_list)
  de_cade <- cade_de_analysis(synth$bulk_expr, synth$group, prop_cade$proportions)

  # 2. Standard limma
  cat("  Running standard limma...\n")
  design_unadj <- model.matrix(~ synth$group)
  fit_unadj <- lmFit(synth$bulk_expr, design_unadj)
  fit_unadj <- eBayes(fit_unadj, trend=TRUE)
  de_unadj <- topTable(fit_unadj, coef=2, number=Inf, adjust.method="BH")

  # 3. Simple marker-mean weight estimation
  cat("  Running simple marker-mean...\n")
  prop_simple <- matrix(NA, nrow=length(synth$marker_list), ncol=ncol(synth$bulk_expr))
  rownames(prop_simple) <- names(synth$marker_list)
  for(ct in names(synth$marker_list)) {
    m <- intersect(synth$marker_list[[ct]], rownames(synth$bulk_expr))
    prop_simple[ct, ] <- colMeans(synth$bulk_expr[m, , drop=FALSE])
  }
  prop_simple_norm <- apply(prop_simple, 2, function(x) {
    e_x <- exp(x - max(x))
    e_x / sum(e_x)
  })

  # 4. CADE with bootstrap
  cat("  Running CADE bootstrap...\n")
  synthetic_key_genes <- synth$de_genes[seq_len(min(10, length(synth$de_genes)))]
  boot_res <- cade_bootstrap(synth$bulk_expr, synth$marker_list,
                             synth$group, n_bootstrap=n_bootstrap,
                             key_genes=synthetic_key_genes)

  # --- Evaluation ---
  # A. Proportion estimation accuracy
  prop_cor_cade <- sapply(rownames(synth$true_props), function(ct) {
    cor(synth$true_props[ct, ], prop_cade$proportions[ct, ], method="pearson")
  })
  prop_cor_simple <- sapply(rownames(synth$true_props), function(ct) {
    cor(synth$true_props[ct, ], prop_simple_norm[ct, ], method="pearson")
  })

  # B. DE gene detection (AUROC)
  scores_cade <- abs(de_cade$logFC.adj)
  names(scores_cade) <- de_cade$Gene
  scores_unadj <- abs(de_unadj$logFC)
  names(scores_unadj) <- rownames(de_unadj)

  library(pROC)
  auroc_cade <- auc(roc(synth$true_de[names(scores_cade)], scores_cade,
                        quiet=TRUE))
  auroc_unadj <- auc(roc(synth$true_de[names(scores_unadj)], scores_unadj,
                         quiet=TRUE))

  # C. CCI for true DE genes vs background
  de_cade_annotated <- de_cade
  de_cade_annotated$TrueDE <- de_cade_annotated$Gene %in% synth$de_genes
  cci_de_true <- mean(de_cade_annotated$CCI[de_cade_annotated$TrueDE & !is.na(de_cade_annotated$CCI)])
  cci_de_false <- mean(de_cade_annotated$CCI[!de_cade_annotated$TrueDE & !is.na(de_cade_annotated$CCI)])

  # Summary
  cat("\n  === Benchmark Results ===\n")
  cat(sprintf("  Proportion correlation (CADE):    %s\n",
      paste(sprintf("%.3f", prop_cor_cade), collapse=", ")))
  cat(sprintf("  Proportion correlation (Simple):   %s\n",
      paste(sprintf("%.3f", prop_cor_simple), collapse=", ")))
  cat(sprintf("  DE detection AUROC (CADE):         %.3f\n", auroc_cade))
  cat(sprintf("  DE detection AUROC (Unadjusted):    %.3f\n", auroc_unadj))
  cat(sprintf("  CCI for true DE genes:             %.3f\n", cci_de_true))
  cat(sprintf("  CCI for background genes:          %.3f\n", cci_de_false))

  list(
    prop_cor_cade = prop_cor_cade,
    prop_cor_simple = prop_cor_simple,
    auroc_cade = auroc_cade,
    auroc_unadj = auroc_unadj,
    cci_de_true = cci_de_true,
    cci_de_false = cci_de_false,
    de_cade = de_cade,
    de_unadj = de_unadj,
    prop_cade = prop_cade,
    prop_simple = prop_simple_norm,
    boot_summary = boot_res$summary
  )
}

# ============================================================
# Part 3: Application to FHL Data (GSE26050)
# ============================================================

get_fhl_case_samples <- function() {
  c(
    "GSM639703", "GSM639704", "GSM639705", "GSM639706", "GSM639707",
    "GSM639708", "GSM639709", "GSM639710", "GSM639711", "GSM639712",
    "GSM639713"
  )
}

get_fhl_group <- function(sample_names) {
  sample_names %in% get_fhl_case_samples()
}

get_fhl_marker_list <- function() {
  list(
    CD8_Tcells = c("CD8A", "CD8B", "PRF1", "GZMB", "GZMA", "GNLY", "NKG7", "CD3E", "CD3D"),
    CD4_Tcells = c("CD4", "CD3E", "CD3D", "IL7R", "CCR7", "LEF1", "TCF7", "SELL"),
    NK_cells = c("NKG7", "GNLY", "PRF1", "GZMB", "KLRB1", "KLRD1", "KLRF1", "NCR1", "CD160"),
    B_cells = c("CD19", "CD79A", "CD79B", "MS4A1", "PAX5", "BLK", "BANK1", "CD22"),
    Monocytes = c("CD14", "FCGR3A", "CSF1R", "ITGAM", "LYZ", "S100A8", "S100A9", "VCAN"),
    Macrophages = c("CD68", "CD163", "MRC1", "MSR1", "CSF1R", "ITGAM", "TLR2", "TLR4"),
    Neutrophils = c("FCGR3B", "CXCR2", "CXCL8", "CSF3R", "MMP8", "MMP9", "ELANE", "MPO"),
    Erythrocytes = c("HBB", "HBA1", "HBA2", "GYPA", "ALAS2", "CA1", "SLC4A1", "AHSP")
  )
}

run_cade_fhl <- function() {
  cat("\n=== Applying CADE to FHL data ===\n")

  # Load real data
  exprs_mat <- as.matrix(read.csv(
    file.path(PROJECT_ROOT, "geo_analysis_output", "GSE26050_expression_corrected.csv"),
    row.names=1, check.names=FALSE
  ))
  cat(sprintf("  Loaded: %d genes x %d samples\n", nrow(exprs_mat), ncol(exprs_mat)))

  # Group assignment
  group <- get_fhl_group(colnames(exprs_mat))

  # Curated marker gene sets
  marker_list <- get_fhl_marker_list()

  # Run CADE
  cat("  Estimating marker-derived weights...\n")
  prop_cade <- estimate_proportions_cade(exprs_mat, marker_list, max_iter=50, tol=1e-6)

  # Compare with simple marker-mean
  cat("  Comparing with simple marker-mean...\n")
  prop_simple <- matrix(NA, nrow=length(marker_list), ncol=ncol(exprs_mat))
  rownames(prop_simple) <- names(marker_list)
  for(ct in names(marker_list)) {
    m <- intersect(marker_list[[ct]], rownames(exprs_mat))
    prop_simple[ct, ] <- colMeans(exprs_mat[m, , drop=FALSE])
  }
  prop_simple_norm <- apply(prop_simple, 2, function(x) {
    e_x <- exp(x - max(x))
    e_x / sum(e_x)
  })

  # DE analysis
  cat("  Running composition-adjusted DE...\n")
  de_cade <- cade_de_analysis(exprs_mat, group, prop_cade$proportions, top_cts=4)

  # Bootstrap
  cat("  Bootstrap uncertainty estimation...\n")
  boot_res <- cade_bootstrap(exprs_mat, marker_list, group, n_bootstrap=200)

  # Ferroptosis gene results
  ferroptosis_genes <- c("SLC7A11", "SLC25A37", "FTH1", "FTL", "GPX4",
                          "GCLM", "TFRC", "HMOX1", "SLC40A1", "NCOA4",
                          "STAT3", "JAK2", "IFNG", "STAT1", "NFE2L2",
                          "IL1B", "TNF", "IL6", "CXCL8", "NFKB1")

  de_cade_ferr <- de_cade[de_cade$Gene %in% ferroptosis_genes, ]
  de_cade_ferr <- de_cade_ferr[order(de_cade_ferr$adj.P.Val.adj), ]

  cat("\n  === CADE Results for Ferroptosis Genes ===\n")
  cat(sprintf("  %-12s %8s %8s %8s %6s %6s %s\n",
              "Gene", "logFC.raw", "logFC.adj", "ΔlogFC", "CCI", "FDR.adj", "Interpretation"))
  for(i in 1:nrow(de_cade_ferr)) {
    r <- de_cade_ferr[i,]
    interp <- if(!is.na(r$CCI)) {
      if(r$CCI < 0.2) "Mainly genuine" else if(r$CCI < 0.5) "Mixed" else "Mainly composition"
    } else "Low expr"
    cat(sprintf("  %-12s %+8.3f %+8.3f %+8.3f %6.3f %6.1e %s\n",
        r$Gene, r$logFC.unadj, r$logFC.adj,
        r$logFC.unadj - r$logFC.adj,
        ifelse(is.na(r$CCI), 0, r$CCI),
        r$adj.P.Val.adj, interp))
  }

  # ============================================================
  # Save results FIRST (before figure generation)
  # ============================================================
  write.csv(de_cade, file.path(OUT_DIR, "Table_CADE_Full_DE_Results.csv"), row.names=FALSE)
  write.csv(de_cade_ferr, file.path(OUT_DIR, "Table_CADE_Ferroptosis_Genes.csv"), row.names=FALSE)
  write.csv(boot_res$summary, file.path(OUT_DIR, "Table_CADE_Bootstrap_Summary.csv"), row.names=FALSE)
  write.csv(t(prop_cade$proportions), file.path(OUT_DIR, "Table_CADE_Proportions.csv"))
  cat("  CSV results saved\n")

  # ============================================================
  # Figures
  # ============================================================

  # Prep data
  prop_comparison <- data.frame()
  ct_names <- names(marker_list)
  for(ct in ct_names) {
    prop_comparison <- rbind(prop_comparison, data.frame(
      Sample = colnames(exprs_mat),
      CellType = ct,
      CADE = prop_cade$proportions[ct, ],
      Simple = prop_simple_norm[ct, ],
      Group = ifelse(group, "FHL", "HC"),
      stringsAsFactors = FALSE
    ))
  }

  # --- Figure 1: Proportion comparison ---
  p1 <- ggplot(prop_comparison, aes(x=Simple, y=CADE, color=Group)) +
    geom_point(alpha=0.6, size=1.5) +
    geom_abline(slope=1, intercept=0, lty=2, alpha=0.3) +
    facet_wrap(~CellType, scales="free") +
    scale_color_manual(values=c("FHL"="#E74C3C", "HC"="#3498DB")) +
    theme_bw(base_size=10) +
    labs(title="CADE vs Simple Marker-Mean Proportion Estimates",
         subtitle="Diagonal = agreement | GSE26050 (FHL n=11, HC n=33)",
         x="Simple Marker-Mean Proportion", y="CADE Estimated Proportion")

  ggsave(file.path(OUT_DIR, "Figure_CADE_Proportion_Comparison.png"), p1,
         width=10, height=8, dpi=300)
  ggsave(file.path(OUT_DIR, "Figure_CADE_Proportion_Comparison.tif"), p1,
         width=10, height=8, dpi=300)

  # --- Figure 2: CCI waterfall ---
  de_ferr_plot <- de_cade_ferr[!is.na(de_cade_ferr$CCI), ]

  p2 <- ggplot(de_ferr_plot, aes(x=reorder(Gene, -CCI), y=CCI)) +
    geom_bar(stat="identity", aes(fill=CCI), color="black", lwd=0.3) +
    geom_text(aes(label=sprintf("%.2f", CCI)), vjust=-0.5, size=3) +
    scale_fill_gradient2(low="#2ECC71", mid="#F39C12", high="#E74C3C",
                         midpoint=0.5, limits=c(0,1)) +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(title="Composition Confounding Index (CCI) for Ferroptosis Genes",
         subtitle="CCI ~ 0: stable after adjustment | CCI ~ 1: highly model-sensitive",
         x="", y="CCI") +
    coord_cartesian(ylim=c(0, 1.2))

  ggsave(file.path(OUT_DIR, "Figure_CADE_CCI_Ferroptosis.png"), p2,
         width=10, height=5, dpi=300)
  ggsave(file.path(OUT_DIR, "Figure_CADE_CCI_Ferroptosis.tif"), p2,
         width=10, height=5, dpi=300)

  # --- Figure 3: logFC adjustment scatter ---
  # Delta_logFC already computed in cade_de_analysis (adj - unadj)
  p3 <- ggplot(de_ferr_plot, aes(x=logFC.unadj, y=logFC.adj, label=Gene)) +
    geom_abline(slope=1, intercept=0, lty=2, alpha=0.3) +
    geom_point(aes(color=CCI), size=3) +
    geom_text(vjust=-0.7, size=2.8) +
    scale_color_gradient2(low="#3498DB", mid="#F39C12", high="#E74C3C",
                          midpoint=0.5, limits=c(0,1)) +
    theme_bw(base_size=10) +
    labs(title="Composition-Adjusted vs Raw logFC for Ferroptosis Genes",
         subtitle="Points below diagonal: raw logFC overestimates due to composition",
         x="Raw logFC (unadjusted)", y="Composition-Adjusted logFC (CADE)")

  ggsave(file.path(OUT_DIR, "Figure_CADE_logFC_Adjustment.png"), p3,
         width=7, height=6, dpi=300)
  ggsave(file.path(OUT_DIR, "Figure_CADE_logFC_Adjustment.tif"), p3,
         width=7, height=6, dpi=300)

  # --- Figure 4: Bootstrap uncertainty ---
  boot_df <- boot_res$summary[!is.na(boot_res$summary$CCI_median), ]
  n_boot_omitted <- nrow(boot_res$summary) - nrow(boot_df)
  if(n_boot_omitted > 0) {
    cat(sprintf("  Note: omitted %d gene(s) with undefined bootstrap CCI from Figure_CADE_Bootstrap_CCI.\n",
                n_boot_omitted))
  }
  p4 <- ggplot(boot_df, aes(x=reorder(Gene, -CCI_median), y=CCI_median)) +
    geom_pointrange(aes(ymin=CCI_lower, ymax=CCI_upper, color=CCI_median),
                    size=0.8) +
    scale_color_gradient2(low="#2ECC71", mid="#F39C12", high="#E74C3C",
                          midpoint=0.5, limits=c(0,1)) +
    theme_bw(base_size=10) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(title="Bootstrap CCI Estimates (200 iterations)",
         subtitle="Error bars: 95% CI from marker gene resampling",
         x="", y="CCI (median +- 95% CI)") +
    coord_cartesian(ylim=c(0, 1.2))

  ggsave(file.path(OUT_DIR, "Figure_CADE_Bootstrap_CCI.png"), p4,
         width=10, height=5, dpi=300)
  ggsave(file.path(OUT_DIR, "Figure_CADE_Bootstrap_CCI.tif"), p4,
         width=10, height=5, dpi=300)

  # --- Figure 5: Marker-derived weight heatmap ---
  prop_heat <- prop_cade$proportions
  # Pre-scale manually (handles zero-variance rows safely)
  prop_heat_scaled <- t(scale(t(prop_heat)))
  prop_heat_scaled[is.na(prop_heat_scaled)] <- 0
  anno_col <- data.frame(Group=ifelse(group, "FHL", "HC"), row.names=colnames(prop_heat_scaled))
  anno_colors <- list(Group=c("FHL"="#E74C3C", "HC"="#3498DB"))

  # Use gplots heatmap.2 as more robust alternative
  suppressPackageStartupMessages(library(gplots))
  tryCatch({
    png(file.path(OUT_DIR, "Figure_CADE_Proportion_Heatmap.png"), width=10, height=5,
        units="in", res=300)
    heatmap.2(prop_heat_scaled, col=colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
              ColSideColors=ifelse(group, "#E74C3C", "#3498DB"),
              Rowv=FALSE, Colv=TRUE, dendrogram="column",
              scale="none", trace="none", key=TRUE, density.info="none",
              margins=c(6, 10), cexRow=0.9, cexCol=0.5,
              main="CADE Marker-Derived Cell-Type Weights\nGSE26050 (row z-score)")
    legend("topright", legend=c("FHL", "HC"), fill=c("#E74C3C", "#3498DB"),
           border=NA, bty="n", cex=0.8)
    dev.off()
  }, error=function(e) cat(sprintf("  Note: heatmap PNG failed: %s\n", e$message)))

  tryCatch({
    tiff(file.path(OUT_DIR, "Figure_CADE_Proportion_Heatmap.tif"), width=10, height=5,
         units="in", res=300)
    heatmap.2(prop_heat_scaled, col=colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
              ColSideColors=ifelse(group, "#E74C3C", "#3498DB"),
              Rowv=FALSE, Colv=TRUE, dendrogram="column",
              scale="none", trace="none", key=TRUE, density.info="none",
              margins=c(6, 10), cexRow=0.9, cexCol=0.5,
              main="CADE Marker-Derived Cell-Type Weights")
    legend("topright", legend=c("FHL", "HC"), fill=c("#E74C3C", "#3498DB"),
           border=NA, bty="n", cex=0.8)
    dev.off()
  }, error=function(e) cat(sprintf("  Note: heatmap TIF failed: %s\n", e$message)))

  # Save marker-derived weight comparison data
  write.csv(prop_comparison, file.path(OUT_DIR, "Table_CADE_Prop_Comparison.csv"), row.names=FALSE)

  cat(sprintf("\n  Output saved to: %s\n", OUT_DIR))

  list(
    de_cade = de_cade,
    de_ferr = de_cade_ferr,
    prop_cade = prop_cade,
    boot_res = boot_res,
    prop_comparison = prop_comparison
  )
}

# ============================================================
# Part 4: Run Everything
# ============================================================

if(identical(Sys.getenv("CADE_SKIP_MAIN", "0"), "1")) {
  cat("CADE_SKIP_MAIN=1 -> core functions loaded; main benchmark/application skipped.\n")
} else {
  if(identical(Sys.getenv("CADE_SKIP_SYNTHETIC_BENCHMARKS", "0"), "1")) {
    cat("\nCADE_SKIP_SYNTHETIC_BENCHMARKS=1 -> synthetic benchmark skipped.\n")
  } else {
    cat("\n========================================\n")
    cat("Part A: Synthetic data benchmark\n")
    cat("========================================\n")

    benchmark_scenarios <- list(
      "Strong_DE_large_bias" = list(de=200, bias=TRUE),
      "Strong_DE_no_bias"    = list(de=200, bias=FALSE),
      "Weak_DE_large_bias"   = list(de=50, bias=TRUE)
    )

    benchmark_results <- list()
    for(scenario in names(benchmark_scenarios)) {
      cat(sprintf("\n--- Scenario: %s ---\n", scenario))
      params <- benchmark_scenarios[[scenario]]

      for(rep in 1:3) {
        cat(sprintf("  Replicate %d/3...\n", rep))
        synth <- generate_synthetic_data(n_genes=5000, n_samples=44, n_ct=5,
                                         ct_bias=params$bias, de_genes=params$de)
        bench <- benchmark_cade(synth, n_bootstrap=50)
        bench$scenario <- scenario
        bench$replicate <- rep
        bench$de_genes <- params$de
        bench$bias <- params$bias
        benchmark_results[[paste(scenario, rep)]] <- bench
      }
    }

    # Aggregate benchmark results
    bench_summary <- do.call(rbind, lapply(benchmark_results, function(b) {
      data.frame(
        Scenario = b$scenario,
        Replicate = b$replicate,
        PropCor_CADE_mean = mean(b$prop_cor_cade),
        PropCor_Simple_mean = mean(b$prop_cor_simple),
        AUROC_CADE = b$auroc_cade,
        AUROC_Unadj = b$auroc_unadj,
        CCI_TrueDE = b$cci_de_true,
        CCI_Background = b$cci_de_false,
        stringsAsFactors = FALSE
      )
    }))

    cat("\n========================================\n")
    cat("Benchmark Summary\n")
    cat("========================================\n")
    print(bench_summary, row.names=FALSE)
    write.csv(bench_summary, file.path(OUT_DIR, "Table_Benchmark_Summary.csv"), row.names=FALSE)

    # Aggregate by scenario
    cat("\n  Average by scenario:\n")
    agg <- bench_summary %>% group_by(Scenario) %>%
      summarise(
        PropCor_CADE = mean(PropCor_CADE_mean),
        PropCor_Simple = mean(PropCor_Simple_mean),
        AUROC_CADE = mean(AUROC_CADE),
        AUROC_Unadj = mean(AUROC_Unadj),
        CCI_TrueDE = mean(CCI_TrueDE),
        CCI_Background = mean(CCI_Background),
        .groups="drop"
      )
    print(agg, row.names=FALSE)
    write.csv(agg, file.path(OUT_DIR, "Table_Benchmark_Aggregated.csv"), row.names=FALSE)
  }

  cat("\n========================================\n")
  cat("Part B: FHL data application\n")
  cat("========================================\n")

  fhl_results <- run_cade_fhl()

  cat("\n=== CADE Analysis Complete ===\n")
  cat("Key finding: CCI ranks genes by composition-adjustment sensitivity\n")
}
