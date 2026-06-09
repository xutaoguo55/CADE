#!/usr/bin/env Rscript
# Empirical comparator and runtime/scalability benchmark for CADE.
#
# This script is self-contained and uses simulated ground truth so every method
# is evaluated against the same known cell-intrinsic DE and composition-shift
# labels. It does not report unavailable external tools as if they were run.

suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(RUVSeq)
  library(pROC)
  library(quadprog)
})

get_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args_full, value = TRUE)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", hit[1]))))
  }
  tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
}

SCRIPT_DIR <- get_script_dir()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
OUT_DIR <- file.path(PROJECT_ROOT, "analysis_output", "benchmark_comparison")
TABLE_DIR <- file.path(PROJECT_ROOT, "tables")
RAW_DIR <- file.path(PROJECT_ROOT, "supplementary", "raw_csv_components")
FIG_DIR <- file.path(PROJECT_ROOT, "figures")
PNG_DIR <- file.path(PROJECT_ROOT, "manuscript", "docx_embedded_figures_png")
for (d in c(OUT_DIR, TABLE_DIR, RAW_DIR, FIG_DIR, PNG_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

set.seed(20260609)

safe_auc <- function(truth, score) {
  ok <- is.finite(score) & !is.na(truth)
  truth <- as.integer(truth[ok])
  score <- as.numeric(score[ok])
  if (length(unique(truth)) < 2 || length(unique(score)) < 2) return(NA_real_)
  as.numeric(pROC::roc(truth, score, quiet = TRUE, direction = "auto")$auc)
}

safe_pr_auc <- function(truth, score) {
  ok <- is.finite(score) & !is.na(truth)
  truth <- as.integer(truth[ok])
  score <- as.numeric(score[ok])
  if (sum(truth == 1) == 0 || sum(truth == 0) == 0) return(NA_real_)
  ord <- order(score, decreasing = TRUE)
  truth <- truth[ord]
  tp <- cumsum(truth == 1)
  fp <- cumsum(truth == 0)
  recall <- tp / sum(truth == 1)
  precision <- tp / pmax(tp + fp, 1)
  recall <- c(0, recall)
  precision <- c(precision[1], precision)
  sum(diff(recall) * precision[-1])
}

softmax_rows <- function(scores) {
  t(apply(scores, 1, function(x) {
    e <- exp(x - max(x, na.rm = TRUE))
    e / sum(e)
  }))
}

close_rows <- function(x, eps = 1e-6) {
  x <- pmax(as.matrix(x), 0) + eps
  x / rowSums(x)
}

ilr_coords <- function(prop_samples_by_ct, eps = 1e-6) {
  closed <- close_rows(prop_samples_by_ct, eps)
  n_comp <- ncol(closed)
  basis <- contr.helmert(n_comp)
  basis <- sweep(basis, 2, sqrt(colSums(basis^2)), "/")
  colnames(basis) <- paste0("ILR", seq_len(ncol(basis)))
  log(closed) %*% basis
}

simulate_ground_truth <- function(n_genes = 3000, n_samples = 48, n_ct = 5,
                                  n_markers_per_ct = 12, n_true_de = 80,
                                  bias = 0.5, marker_noise = 0,
                                  seed = 1) {
  set.seed(seed)
  n_half <- n_samples / 2
  genes <- paste0("Gene", seq_len(n_genes))
  ct_names <- paste0("CT", seq_len(n_ct))

  marker_count <- n_ct * n_markers_per_ct
  marker_idx <- split(seq_len(marker_count), rep(seq_len(n_ct), each = n_markers_per_ct))
  marker_genes <- lapply(marker_idx, function(ii) genes[ii])
  names(marker_genes) <- ct_names

  ct_specific <- seq_len(min(500, n_genes - n_true_de - 20))
  available_de <- setdiff(seq_len(n_genes), ct_specific)
  true_de <- sample(available_de, n_true_de)

  profile <- matrix(rnorm(n_genes * n_ct, mean = 8, sd = 0.55), n_genes, n_ct,
                    dimnames = list(genes, ct_names))
  for (k in seq_len(n_ct)) {
    profile[marker_idx[[k]], k] <- profile[marker_idx[[k]], k] + 4.0
    profile[marker_idx[[k]], setdiff(seq_len(n_ct), k)] <-
      profile[marker_idx[[k]], setdiff(seq_len(n_ct), k)] - 1.0
  }
  ct_specific_non_marker <- setdiff(ct_specific, unlist(marker_idx))
  ct_bins <- split(ct_specific_non_marker, rep(seq_len(n_ct), length.out = length(ct_specific_non_marker)))
  for (k in seq_len(n_ct)) {
    ii <- ct_bins[[k]]
    profile[ii, k] <- profile[ii, k] + 2.8
    profile[ii, setdiff(seq_len(n_ct), k)] <- profile[ii, setdiff(seq_len(n_ct), k)] - 1.2
  }

  p_ctrl <- rep(1 / n_ct, n_ct)
  p_case_max <- c(0.46, 0.08, 0.07, 0.27, 0.12)
  p_case_max <- p_case_max / sum(p_case_max)
  p_case <- p_ctrl + bias * (p_case_max - p_ctrl)
  p_case <- p_case / sum(p_case)

  props <- matrix(NA_real_, n_samples, n_ct, dimnames = list(paste0("S", seq_len(n_samples)), ct_names))
  for (i in seq_len(n_samples)) {
    center <- if (i <= n_half) p_ctrl else p_case
    jitter <- rgamma(n_ct, shape = pmax(center * 80, 1), rate = 80)
    props[i, ] <- jitter / sum(jitter)
  }

  bulk <- t(profile %*% t(props))
  disease_effect <- rep(0, n_genes)
  disease_effect[true_de] <- rnorm(n_true_de, mean = 1.3, sd = 0.25)
  bulk[(n_half + 1):n_samples, ] <- sweep(
    bulk[(n_half + 1):n_samples, , drop = FALSE], 2, disease_effect, "+"
  )
  bulk <- bulk + matrix(rnorm(n_samples * n_genes, 0, 0.20), n_samples, n_genes)
  bulk <- t(bulk)
  rownames(bulk) <- genes
  colnames(bulk) <- paste0("S", seq_len(n_samples))

  marker_list <- marker_genes
  if (marker_noise > 0) {
    n_swap <- max(1, floor(length(unlist(marker_list)) * marker_noise))
    random_genes <- sample(setdiff(genes, unlist(marker_list)), n_swap)
    flat <- unlist(marker_list)
    flat[sample(seq_along(flat), n_swap)] <- random_genes
    marker_list <- split(flat, rep(ct_names, each = n_markers_per_ct))
  }

  mismatch_profile <- profile + matrix(rnorm(length(profile), 0, 1.1), nrow(profile), ncol(profile))
  mismatch_profile[, sample(seq_len(n_ct))] <- mismatch_profile

  list(
    bulk = bulk,
    group = factor(rep(c("Control", "Case"), each = n_half), levels = c("Control", "Case")),
    marker_list = marker_list,
    true_props = props,
    matched_profile = profile,
    mismatched_profile = mismatch_profile,
    true_de = genes[true_de],
    composition_sensitive = genes[setdiff(ct_specific, true_de)]
  )
}

marker_weights <- function(bulk, marker_list) {
  scores <- sapply(marker_list, function(gs) colMeans(bulk[intersect(gs, rownames(bulk)), , drop = FALSE]))
  colnames(scores) <- names(marker_list)
  softmax_rows(scores)
}

qp_weights_from_reference <- function(bulk, reference_profile) {
  genes <- intersect(rownames(bulk), rownames(reference_profile))
  Y <- bulk[genes, , drop = FALSE]
  S <- reference_profile[genes, , drop = FALSE]
  n_ct <- ncol(S)
  D <- t(S) %*% S
  diag(D) <- diag(D) + 1e-6
  out <- matrix(NA_real_, ncol(Y), n_ct, dimnames = list(colnames(Y), colnames(S)))
  for (i in seq_len(ncol(Y))) {
    sol <- tryCatch(
      quadprog::solve.QP(Dmat = D, dvec = as.vector(t(S) %*% Y[, i]),
                         Amat = cbind(rep(1, n_ct), diag(n_ct)),
                         bvec = c(1, rep(0, n_ct)), meq = 1)$solution,
      error = function(e) rep(1 / n_ct, n_ct)
    )
    out[i, ] <- pmax(sol, 0) / sum(pmax(sol, 0))
  }
  out
}

fit_limma_score <- function(bulk, group, covars = NULL) {
  if (is.null(covars) || ncol(as.matrix(covars)) == 0) {
    design <- model.matrix(~ group)
  } else {
    covars <- as.matrix(covars)
    covars <- covars[, apply(covars, 2, var, na.rm = TRUE) > 1e-8, drop = FALSE]
    if (ncol(covars) > 0) {
      keep <- seq_len(ncol(covars))
      repeat {
        design <- model.matrix(~ group + covars[, keep, drop = FALSE])
        if (qr(design)$rank == ncol(design) || length(keep) <= 1) break
        keep <- keep[-length(keep)]
      }
    } else {
      design <- model.matrix(~ group)
    }
  }
  tt <- limma::topTable(limma::eBayes(limma::lmFit(bulk, design), trend = TRUE),
                        coef = 2, number = Inf, sort.by = "none")
  list(t = abs(tt$t), logFC = tt$logFC)
}

run_sva_score <- function(bulk, group) {
  mod <- model.matrix(~ group)
  mod0 <- model.matrix(~ 1, data.frame(group = group))
  sv <- tryCatch(sva::sva(bulk, mod, mod0, n.sv = 2)$sv, error = function(e) NULL)
  fit_limma_score(bulk, group, sv)
}

run_ruv_score <- function(bulk, group) {
  base <- fit_limma_score(bulk, group)
  ctrl <- order(base$t, decreasing = FALSE)[seq_len(min(500, length(base$t)))]
  W <- tryCatch(RUVSeq::RUVg(bulk, ctrl, k = 2)$W, error = function(e) NULL)
  fit_limma_score(bulk, group, W)
}

cade_lite_score <- function(bulk, group, marker_list, top_cts = 4, transform = "raw") {
  w <- marker_weights(bulk, marker_list)
  cov <- if (identical(transform, "ilr")) ilr_coords(w) else w
  vars <- apply(cov, 2, var, na.rm = TRUE)
  cov <- cov[, names(sort(vars, decreasing = TRUE))[seq_len(min(top_cts, ncol(cov)))], drop = FALSE]
  unadj <- fit_limma_score(bulk, group)
  adj <- fit_limma_score(bulk, group, cov)
  cci <- pmin(pmax(abs(unadj$logFC - adj$logFC) / pmax(abs(unadj$logFC), 0.1), 0), 1)
  list(t = adj$t, logFC = adj$logFC, CCI = cci, weights = w)
}

cade_full_weights <- function(bulk, marker_list, max_iter = 6, tol = 1e-5) {
  w <- t(marker_weights(bulk, marker_list))
  ct_names <- rownames(w)
  marker_present <- lapply(marker_list, intersect, rownames(bulk))
  all_markers <- unique(unlist(marker_present))
  marker_ct_map <- matrix(0, length(ct_names), length(all_markers),
                          dimnames = list(ct_names, all_markers))
  for (ct in ct_names) marker_ct_map[ct, marker_present[[ct]]] <- 1
  for (iter in seq_len(max_iter)) {
    old <- w
    profiles <- matrix(0, length(all_markers), length(ct_names),
                       dimnames = list(all_markers, ct_names))
    for (g in all_markers) {
      y <- bulk[g, ]
      for (ct in ct_names) {
        if (marker_ct_map[ct, g] == 1 && sd(w[ct, ]) > 1e-8) {
          profiles[g, ct] <- coef(lm(y ~ w[ct, ]))[2]
        }
      }
    }
    profiles[!is.finite(profiles)] <- 0
    D <- t(profiles) %*% profiles
    diag(D) <- diag(D) + 1e-6
    for (i in seq_len(ncol(bulk))) {
      sol <- tryCatch(
        quadprog::solve.QP(D, as.vector(t(profiles) %*% bulk[all_markers, i]),
                           cbind(rep(1, length(ct_names)), diag(length(ct_names))),
                           c(1, rep(0, length(ct_names))), meq = 1)$solution,
        error = function(e) old[, i]
      )
      w[, i] <- pmax(sol, 0) / sum(pmax(sol, 0))
    }
    if (max(abs(w - old)) < tol) break
  }
  t(w)
}

timed <- function(expr) {
  gc()
  start <- proc.time()[["elapsed"]]
  value <- force(expr)
  elapsed <- proc.time()[["elapsed"]] - start
  list(value = value, elapsed = elapsed)
}

benchmark_comparators <- function() {
  bias_levels <- c(0, 25, 50, 100)
  reps <- 3
  rows <- list()
  i <- 0
  for (bias in bias_levels) {
    for (rep_i in seq_len(reps)) {
      dat <- simulate_ground_truth(bias = bias / 100, seed = 9000 + bias * 10 + rep_i)
      truth_de <- rownames(dat$bulk) %in% dat$true_de
      truth_comp <- rownames(dat$bulk) %in% dat$composition_sensitive

      methods <- list()
      methods[["limma"]] <- timed(fit_limma_score(dat$bulk, dat$group))
      mw <- timed(marker_weights(dat$bulk, dat$marker_list))
      methods[["limma+marker"]] <- timed(fit_limma_score(dat$bulk, dat$group, mw$value))
      methods[["limma+SVA"]] <- timed(run_sva_score(dat$bulk, dat$group))
      methods[["limma+RUVg"]] <- timed(run_ruv_score(dat$bulk, dat$group))
      matched <- timed(qp_weights_from_reference(dat$bulk, dat$matched_profile))
      mismatched <- timed(qp_weights_from_reference(dat$bulk, dat$mismatched_profile))
      methods[["matched-NNLS+limma"]] <- timed(fit_limma_score(dat$bulk, dat$group, matched$value))
      methods[["mismatched-NNLS+limma"]] <- timed(fit_limma_score(dat$bulk, dat$group, mismatched$value))
      methods[["oracle-proportion+limma"]] <- timed(fit_limma_score(dat$bulk, dat$group, dat$true_props))
      methods[["CADE-lite"]] <- timed(cade_lite_score(dat$bulk, dat$group, dat$marker_list, transform = "raw"))
      full_w <- timed(cade_full_weights(dat$bulk, dat$marker_list, max_iter = 6))
      methods[["CADE-full"]] <- timed({
        cov <- full_w$value
        vars <- apply(cov, 2, var)
        cov <- cov[, names(sort(vars, decreasing = TRUE))[1:min(4, ncol(cov))], drop = FALSE]
        unadj <- fit_limma_score(dat$bulk, dat$group)
        adj <- fit_limma_score(dat$bulk, dat$group, cov)
        cci <- pmin(pmax(abs(unadj$logFC - adj$logFC) / pmax(abs(unadj$logFC), 0.1), 0), 1)
        list(t = adj$t, logFC = adj$logFC, CCI = cci)
      })
      methods[["CADE-ILR"]] <- timed(cade_lite_score(dat$bulk, dat$group, dat$marker_list, transform = "ilr"))

      for (method in names(methods)) {
        result <- methods[[method]]$value
        tscore <- result$t
        cci <- result$CCI
        top100 <- order(tscore, decreasing = TRUE)[seq_len(min(100, length(tscore)))]
        comp_fp_top100 <- mean(truth_comp[top100])
        cci_auc <- if (!is.null(cci)) {
          keep <- truth_de | truth_comp
          safe_auc(truth_de[keep], 1 - cci[keep])
        } else {
          NA_real_
        }
        i <- i + 1
        rows[[i]] <- data.frame(
          bias_percent = bias,
          replicate = rep_i,
          method = method,
          auroc_true_de = safe_auc(truth_de, tscore),
          auprc_true_de = safe_pr_auc(truth_de, tscore),
          composition_fp_fraction_top100 = comp_fp_top100,
          cci_stable_vs_composition_auc = cci_auc,
          runtime_sec = methods[[method]]$elapsed,
          stringsAsFactors = FALSE
        )
      }
      cat(sprintf("Comparator benchmark bias=%d%% replicate=%d complete\n", bias, rep_i))
    }
  }
  do.call(rbind, rows)
}

runtime_scalability <- function() {
  scenarios <- data.frame(
    n_genes = c(1000, 3000, 6000, 3000, 3000),
    n_samples = c(24, 48, 72, 48, 48),
    bootstrap_iterations = c(0, 0, 0, 25, 50),
    stringsAsFactors = FALSE
  )
  rows <- list()
  j <- 0
  for (s in seq_len(nrow(scenarios))) {
    dat <- simulate_ground_truth(
      n_genes = scenarios$n_genes[s],
      n_samples = scenarios$n_samples[s],
      bias = 0.5,
      seed = 7000 + s
    )
    for (mode in c("CADE-lite", "CADE-full", "CADE-ILR")) {
      res <- timed({
        if (mode == "CADE-lite") {
          fit <- cade_lite_score(dat$bulk, dat$group, dat$marker_list, transform = "raw")
        } else if (mode == "CADE-ILR") {
          fit <- cade_lite_score(dat$bulk, dat$group, dat$marker_list, transform = "ilr")
        } else {
          fw <- cade_full_weights(dat$bulk, dat$marker_list, max_iter = 6)
          fit <- fit_limma_score(dat$bulk, dat$group, fw)
        }
        if (scenarios$bootstrap_iterations[s] > 0 && mode == "CADE-lite") {
          key <- rownames(dat$bulk)[seq_len(min(30, nrow(dat$bulk)))]
          boot_cci <- replicate(scenarios$bootstrap_iterations[s], {
            sub_markers <- lapply(dat$marker_list, function(gs) {
              n_keep <- max(3, floor(length(gs) * 0.8))
              sample(gs, n_keep)
            })
            tmp <- cade_lite_score(dat$bulk, dat$group, sub_markers, transform = "raw")
            tmp$CCI[key]
          })
          fit$bootstrap_summary <- rowMeans(boot_cci, na.rm = TRUE)
        }
        fit
      })
      j <- j + 1
      rows[[j]] <- data.frame(
        n_genes = scenarios$n_genes[s],
        n_samples = scenarios$n_samples[s],
        bootstrap_iterations = ifelse(mode == "CADE-lite", scenarios$bootstrap_iterations[s], 0),
        mode = mode,
        runtime_sec = res$elapsed,
        object_footprint_mb = as.numeric(object.size(list(dat$bulk, res$value))) / 1024^2,
        stringsAsFactors = FALSE
      )
      cat(sprintf("Runtime benchmark genes=%d samples=%d boot=%d mode=%s complete\n",
                  scenarios$n_genes[s], scenarios$n_samples[s],
                  ifelse(mode == "CADE-lite", scenarios$bootstrap_iterations[s], 0), mode))
    }
  }
  do.call(rbind, rows)
}

write_captioned_csv <- function(df, path, caption) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste0("Table_Caption,", caption), con)
  utils::write.table(df, con, sep = ",", row.names = FALSE, col.names = TRUE, qmethod = "double")
}

plot_summary <- function(comp_summary, runtime_df) {
  tif_path <- file.path(FIG_DIR, "SuppFigure_S11_Comparator_Runtime_Benchmark.tif")
  png_path <- file.path(PNG_DIR, "SuppFigure_S11_Comparator_Runtime_Benchmark.png")
  draw <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(mfrow = c(2, 2), mar = c(4.6, 4.6, 2.5, 1.2), oma = c(0, 0, 2, 0))
    methods <- unique(comp_summary$method)
    cols <- setNames(rainbow(length(methods), s = 0.65, v = 0.75), methods)
    plot(NA, xlim = range(comp_summary$bias_percent), ylim = c(0.45, 1.02),
         xlab = "Composition bias (%)", ylab = "AUROC for true DE",
         main = "A. Comparator DE recovery")
    for (m in methods) {
      d <- comp_summary[comp_summary$method == m, ]
      lines(d$bias_percent, d$auroc_true_de_mean, type = "b", pch = 16, col = cols[m], lwd = 1.8)
    }
    legend("bottomleft", legend = methods, col = cols[methods], lty = 1, pch = 16, cex = 0.58, bty = "n")

    plot(NA, xlim = range(comp_summary$bias_percent), ylim = c(0, 1),
         xlab = "Composition bias (%)", ylab = "Composition-driven genes in top 100",
         main = "B. False composition prioritisation")
    for (m in methods) {
      d <- comp_summary[comp_summary$method == m, ]
      lines(d$bias_percent, d$composition_fp_fraction_top100_mean, type = "b", pch = 16, col = cols[m], lwd = 1.8)
    }

    cci <- comp_summary[is.finite(comp_summary$cci_stable_vs_composition_auc_mean), ]
    plot(NA, xlim = range(cci$bias_percent), ylim = c(0.45, 1.02),
         xlab = "Composition bias (%)", ylab = "CCI AUROC",
         main = "C. CADE CCI stable-vs-composition ranking")
    for (m in unique(cci$method)) {
      d <- cci[cci$method == m, ]
      lines(d$bias_percent, d$cci_stable_vs_composition_auc_mean, type = "b", pch = 16, col = cols[m], lwd = 1.8)
    }

    rt <- runtime_df[runtime_df$bootstrap_iterations == 0, ]
    rt_modes <- unique(rt$mode)
    rt_cols <- setNames(c("#2b8cbe", "#b2182b", "#4daf4a")[seq_along(rt_modes)], rt_modes)
    plot(NA, xlim = range(rt$n_genes), ylim = c(0, max(rt$runtime_sec) * 1.15),
         xlab = "Genes", ylab = "Runtime (sec)", main = "D. Runtime scalability")
    for (m in rt_modes) {
      d <- rt[rt$mode == m, ]
      d <- d[order(d$n_genes), ]
      lines(d$n_genes, d$runtime_sec, type = "b", pch = 16, col = rt_cols[m], lwd = 1.8)
    }
    legend("topleft", legend = rt_modes, col = rt_cols[rt_modes],
           lty = 1, pch = 16, cex = 0.7, bty = "n")
    mtext("Supplementary Figure S11 | Empirical comparator and runtime/scalability benchmark",
          outer = TRUE, cex = 1.1, font = 2)
  }
  tiff(tif_path, width = 8.5, height = 7.0, units = "in", res = 300, compression = "lzw")
  draw()
  dev.off()
  png(png_path, width = 8.5, height = 7.0, units = "in", res = 300)
  draw()
  dev.off()
  cat(sprintf("Wrote %s\n", tif_path))
  cat(sprintf("Wrote %s\n", png_path))
}

cat("=== CADE empirical comparator and runtime benchmark ===\n")
cat(sprintf("Project root: %s\n", PROJECT_ROOT))
comp_raw <- benchmark_comparators()
runtime_raw <- runtime_scalability()

summarise_metric <- function(metric) {
  out <- aggregate(comp_raw[[metric]], comp_raw[, c("bias_percent", "method")],
                   function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE))
  names(out)[3] <- paste0(metric, "_mean")
  out
}
metric_tables <- lapply(
  c("auroc_true_de", "auprc_true_de", "composition_fp_fraction_top100",
    "cci_stable_vs_composition_auc", "runtime_sec"),
  summarise_metric
)
comp_summary <- Reduce(
  function(x, y) merge(x, y, by = c("bias_percent", "method"), all = TRUE),
  metric_tables
)
comp_summary <- comp_summary[order(comp_summary$bias_percent, comp_summary$method), ]

write.csv(comp_raw, file.path(OUT_DIR, "empirical_comparator_benchmark_raw.csv"), row.names = FALSE)
write.csv(comp_summary, file.path(OUT_DIR, "empirical_comparator_benchmark_summary.csv"), row.names = FALSE)
write.csv(runtime_raw, file.path(OUT_DIR, "runtime_scalability_benchmark.csv"), row.names = FALSE)
write.csv(comp_raw, file.path(RAW_DIR, "Table_S46_EmpiricalComparator_Raw.csv"), row.names = FALSE)
write.csv(runtime_raw, file.path(RAW_DIR, "Table_S47_RuntimeScalability_Raw.csv"), row.names = FALSE)
write_captioned_csv(
  comp_summary,
  file.path(TABLE_DIR, "Table_01F_Empirical_Comparator_Benchmark.csv"),
  "Table 1F: Empirical comparator benchmark on simulated ground truth (mean of 3 replicates)"
)
write_captioned_csv(
  runtime_raw,
  file.path(TABLE_DIR, "Table_01G_Runtime_Scalability.csv"),
  "Table 1G: Runtime and object-footprint scalability benchmark"
)
plot_summary(comp_summary, runtime_raw)

cat("=== DONE ===\n")
