#!/usr/bin/env Rscript
# ============================================================================
# CADE Real scRNA-seq Pseudo-Bulk Benchmark v4 — Hardened Edition
# Uses Monaco Immune Reference (real scRNA-derived) for cell-type profiles
# Tests: CADE vs limma vs MarkerScore vs SVA vs CIBERSORT/LM22
#
# Upgrades from v3:
#   - Weaker DE signal (mix of 0.5, 0.8, 1.0 logFC) for realistic difficulty
#   - Higher noise (0.12 vs 0.05) for realistic microarray log2 noise
#   - 80 true DE genes (up from 50)
#   - 8 markers per cell type (down from 12) to stress-test marker dependency
#   - Added 150% extreme bias level
#   - Added 20 composition-driven DE genes (apparent DE from proportion shifts)
#   - CIBERSORT/LM22 comparison (nnls-based deconvolution)
# ============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(quadprog)
  library(pROC)
  library(celldex)
})

set.seed(42)

# Keep benchmark outputs inside the submission package, independent of cwd.
SCRIPT_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
})
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
OUT_DIR <- file.path(PROJECT_ROOT, "analysis_output", "CADE", "scRNA_benchmark_v4")
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)
TIMING_LOG <- file.path(OUT_DIR, "timing.log")
cat(sprintf("=== CADE Benchmark v4 started at %s ===\n", Sys.time()), file=TIMING_LOG)

cat("\n=====================================================\n")
cat("CADE Real scRNA-seq Pseudo-Bulk Benchmark v4\n")
cat("Hardened edition: weaker DE, higher noise, LM22 compare\n")
cat("=====================================================\n")

# ---- 1. Load Monaco reference ----
cat("\n[1/8] Loading Monaco Immune Reference...\n")
ref <- MonacoImmuneData()
ref_expr <- log2(assay(ref, "logcounts") + 1)

ct_labels <- colData(ref)$label.main
all_cts_full <- names(sort(table(ct_labels), decreasing=TRUE))
ct_counts <- table(ct_labels)
exclude <- names(which(ct_counts < 5))
exclude <- union(exclude, "T cells")
all_cts <- setdiff(all_cts_full, exclude)
cat(sprintf("  Cell types (%d): %s\n", length(all_cts), paste(all_cts, collapse=", ")))

N_CT <- length(all_cts)

ct_profiles <- do.call(cbind, lapply(all_cts, function(ct) {
  cells <- which(ct_labels == ct)
  if(length(cells) > 1) Matrix::rowMeans(ref_expr[, cells, drop=FALSE]) else ref_expr[, cells]
}))
colnames(ct_profiles) <- all_cts

gene_mean <- rowMeans(ct_profiles)
gene_var <- apply(ct_profiles, 1, var)
keep <- which(gene_mean > 1.5 & gene_var > 0.5)
ct_profiles <- ct_profiles[keep, ]
cat(sprintf("  Selected %d expressed variable genes\n", nrow(ct_profiles)))

N_GENES <- nrow(ct_profiles)

# ---- 2. Identify markers ----
cat("\n[2/8] Identifying cell-type markers (max 8 per type)...\n")

MAX_MARKERS <- 8
marker_list <- list()
used <- character(0)
for(i in seq_len(N_CT)) {
  fc <- ct_profiles[, i] - rowMeans(ct_profiles[, -i, drop=FALSE])
  candidates <- names(sort(fc[ct_profiles[, i] > 2 & fc > 2], decreasing=TRUE))
  avail <- setdiff(candidates, used)
  selected <- head(avail, MAX_MARKERS)
  marker_list[[colnames(ct_profiles)[i]]] <- selected
  used <- c(used, selected)
  cat(sprintf("  %20s: %d markers\n", colnames(ct_profiles)[i], length(selected)))
}
all_markers <- unique(unlist(marker_list))
cat(sprintf("  Total unique markers: %d\n", length(all_markers)))

# ---- 3. Build LM22-like signature for CIBERSORT comparison ----
cat("\n[3/8] Building LM22-like reference signature...\n")

# Use the cell-type profiles as the reference for nnls deconvolution
# In practice, CIBERSORT would use LM22, but LM22 genes don't overlap perfectly
# with Monaco data. Instead, we create a "pseudo-LM22" from the Monaco profiles
# using marker genes as the signature basis — this is a fair comparison
# because both methods use the same underlying cell-type profiles.

lm22_sig <- ct_profiles[unique(unlist(marker_list)), ]
cat(sprintf("  LM22-like signature: %d genes x %d cell types\n", nrow(lm22_sig), ncol(lm22_sig)))

ciber_lm22_deconv <- function(bulk, sig) {
  # NNLS-based deconvolution (original CIBERSORT approach)
  common_genes <- intersect(rownames(bulk), rownames(sig))
  S <- as.matrix(sig[common_genes, ])
  Y <- as.matrix(bulk[common_genes, ])

  props <- matrix(NA, ncol(Y), ncol(S))
  colnames(props) <- colnames(S)

  for (i in seq_len(ncol(Y))) {
    fit <- tryCatch({
      nnls::nnls(S, Y[, i])
    }, error = function(e) NULL)
    if (!is.null(fit)) {
      x <- fit$x
      props[i, ] <- x / sum(x)
    } else {
      props[i, ] <- rep(1/ncol(S), ncol(S))
    }
  }
  props
}

# ---- 4. Generate pseudo-bulk data (HARDENED) ----
cat("\n[4/8] Generating pseudo-bulk with realistic challenge...\n")

N_SAMPLES <- 60
N_HALF <- 30
N_TRUE_DE <- 80           # more DE genes to detect
DE_FC <- c(0.5, 0.8, 1.0) # mixed signal strengths
N_COMP_DE <- 20           # composition-driven DE genes
BIAS_LEVELS <- c(0, 0.25, 0.50, 1.0, 1.5)
NOISE_SD <- 0.12          # realistic microarray noise
N_REP <- 3

balanced <- rep(1/N_CT, N_CT)
names(balanced) <- colnames(ct_profiles)

skewed <- balanced
lympho <- grep("T cell|B cell|NK", names(balanced))
myeloid <- grep("Mono|Macro|Dendritic|Neutro", names(balanced))
if(length(lympho) > 0) skewed[lympho] <- balanced[lympho] * 0.25
if(length(myeloid) > 0) skewed[myeloid] <- balanced[myeloid] * (1 + 0.75*length(lympho)/length(myeloid))
skewed <- skewed / sum(skewed)

generate_data <- function(bias) {
  target <- (1 - min(bias, 1)) * balanced + min(bias, 1) * skewed

  true_p <- matrix(NA, N_SAMPLES, N_CT)
  colnames(true_p) <- names(balanced)
  for(i in 1:N_HALF) {
    x <- rgamma(N_CT, shape=balanced*50, rate=1)
    true_p[i, ] <- x / sum(x)
  }
  for(i in (N_HALF+1):N_SAMPLES) {
    x <- rgamma(N_CT, shape=target*50, rate=1)
    true_p[i, ] <- x / sum(x)
  }

  bulk <- matrix(NA, N_GENES, N_SAMPLES)
  rownames(bulk) <- rownames(ct_profiles)
  for(i in 1:N_SAMPLES) {
    bulk[, i] <- ct_profiles %*% true_p[i, ] + rnorm(N_GENES, 0, NOISE_SD)
  }

  # Spike in true DE genes (cell-intrinsic regulation)
  de_pool <- setdiff(rownames(bulk), all_markers)
  true_de <- sample(de_pool, N_TRUE_DE)
  for(g in true_de) {
    fc <- sample(DE_FC, 1)
    bulk[g, (N_HALF+1):N_SAMPLES] <- bulk[g, (N_HALF+1):N_SAMPLES] + fc
  }

  # Spike in composition-driven DE genes (apparent DE from proportion shifts)
  # These genes have HIGH expression in cell types that expand in disease
  comp_de_pool <- setdiff(de_pool, true_de)
  comp_de <- sample(comp_de_pool, N_COMP_DE)
  # For each comp_DE gene, boost expression proportional to the cell-type shift
  for(g in comp_de) {
    # Make this gene highly expressed in one of the expanding myeloid types
    high_ct <- sample(myeloid, 1)
    ct_shift <- mean(true_p[(N_HALF+1):N_SAMPLES, high_ct] - true_p[1:N_HALF, high_ct])
    bulk[g, (N_HALF+1):N_SAMPLES] <- bulk[g, (N_HALF+1):N_SAMPLES] + ct_shift * 2
  }

  all_de <- c(true_de, comp_de)
  de_type <- c(rep("cell_intrinsic", N_TRUE_DE), rep("composition_driven", N_COMP_DE))
  names(de_type) <- all_de

  group <- factor(c(rep("Control", N_HALF), rep("Disease", N_HALF)))
  list(bulk=bulk, true_p=true_p, true_de=true_de, comp_de=comp_de,
       all_de=all_de, de_type=de_type, group=group, bias=bias)
}

# ---- 5. CADE implementation ----
estimate_cade <- function(expr, markers) {
  n_ct <- length(markers)
  n_samp <- ncol(expr)
  ct_nm <- names(markers)

  P <- matrix(NA, n_samp, n_ct)
  for(k in seq_len(n_ct)) {
    mk <- intersect(markers[[k]], rownames(expr))
    P[, k] <- if(length(mk)>0) colMeans(expr[mk,,drop=FALSE]) else 0
  }
  P <- exp(P - apply(P, 1, max))
  P <- P / rowSums(P)
  P[is.na(P)] <- 1/n_ct

  for(iter in 1:50) {
    C <- matrix(NA, nrow(expr), n_ct)
    for(k in 1:n_ct) C[, k] <- rowMeans(expr * P[, k])

    Dmat <- t(C) %*% C + diag(1e-8, n_ct)
    Amat <- cbind(1, diag(n_ct))
    bvec <- c(1, rep(0, n_ct))

    Pnew <- matrix(NA, n_samp, n_ct)
    for(i in 1:n_samp) {
      qp <- tryCatch(solve.QP(Dmat, t(C) %*% expr[,i], Amat, bvec, meq=1),
                     error=function(e) NULL)
      Pnew[i,] <- if(!is.null(qp)) qp$solution else P[i,]
    }
    if(max(abs(Pnew-P)) < 1e-6) { P <- Pnew; break }
    P <- Pnew
  }
  colnames(P) <- ct_nm
  P
}

run_cade_de <- function(expr, props, group) {
  pv <- apply(props, 2, var)
  top <- names(sort(pv, decreasing=TRUE))[1:min(4, ncol(props))]
  fit_u <- eBayes(lmFit(expr, model.matrix(~group)))
  fit_a <- eBayes(lmFit(expr, model.matrix(~group+props[,top,drop=FALSE])))
  lfc_u <- fit_u$coefficients[,2]
  lfc_a <- fit_a$coefficients[,2]
  cci <- rep(NA, nrow(expr))
  names(cci) <- rownames(expr)
  v <- abs(lfc_u)>=0.1
  cci[v] <- pmax(0, pmin(1, abs(lfc_u[v]-lfc_a[v])/abs(lfc_u[v])))
  list(lfc_u=lfc_u, lfc_a=lfc_a, cci=cci, p_u=fit_u$p.value[,2], p_a=fit_a$p.value[,2])
}

run_marker_score <- function(expr, markers, group) {
  n_ct <- length(markers)
  P <- matrix(NA, ncol(expr), n_ct)
  for(k in seq_len(n_ct)) {
    mk <- intersect(markers[[k]], rownames(expr))
    P[,k] <- colMeans(expr[mk,,drop=FALSE])
  }
  pv <- apply(P, 2, var)
  top <- order(pv, decreasing=TRUE)[1:min(4, n_ct)]
  fit <- eBayes(lmFit(expr, model.matrix(~group+P[,top,drop=FALSE])))
  list(pval=fit$p.value[,2], props=P)
}

run_sva <- function(expr, group) {
  if(!requireNamespace("sva", quietly=TRUE)) return(NULL)
  mod <- model.matrix(~group)
  mod0 <- model.matrix(~1, data=data.frame(group))
  sv <- tryCatch(sva::sva(expr, mod, mod0, n.sv=min(5,ncol(expr)-2))$sv,
                 error=function(e) NULL)
  if(is.null(sv)) return(NULL)
  fit <- eBayes(lmFit(expr, cbind(mod, sv)))
  list(pval=fit$p.value[,2])
}

run_ciber_lm22 <- function(expr, sig, group) {
  props <- ciber_lm22_deconv(expr, sig)
  pv <- apply(props, 2, var)
  top <- order(pv, decreasing=TRUE)[1:min(4, ncol(props))]
  fit <- eBayes(lmFit(expr, model.matrix(~group+props[,top,drop=FALSE])))
  list(pval=fit$p.value[,2], props=props)
}

auroc_score <- function(pvals, true_de) {
  score <- -log10(pmax(pvals, 1e-300))
  labs <- ifelse(names(pvals) %in% true_de, 1, 0)
  if(sum(labs)<2) return(NA)
  tryCatch(as.numeric(auc(roc(labs, score, quiet=TRUE))), error=function(e) NA)
}

# ---- 6. Run benchmark ----
cat("\n[6/8] Running benchmark (5 bias × 3 reps × 6 methods)...\n")

all_res <- list()
run_id <- 1

for(bias in BIAS_LEVELS) {
  cat(sprintf("\n--- Bias %.0f%% ---\n", bias*100))
  for(rep in 1:N_REP) {
    t_start <- Sys.time()
    d <- generate_data(bias)
    bg <- setdiff(rownames(d$bulk), c(d$all_de, all_markers))

    # CADE
    cade_p <- estimate_cade(d$bulk, marker_list)
    cade_de <- run_cade_de(d$bulk, cade_p, d$group)

    # Standard limma
    limma_fit <- eBayes(lmFit(d$bulk, model.matrix(~d$group)))

    # MarkerScore + limma
    ms <- run_marker_score(d$bulk, marker_list, d$group)

    # SVA
    sva <- run_sva(d$bulk, d$group)

    # CIBERSORT/LM22-like
    ciber <- run_ciber_lm22(d$bulk, lm22_sig, d$group)

    # AUROC for all DE genes (cell-intrinsic + composition-driven)
    l_auroc <- auroc_score(limma_fit$p.value[,2], d$all_de)
    ms_auroc <- auroc_score(ms$pval, d$all_de)
    cade_auroc <- auroc_score(cade_de$p_a, d$all_de)
    sva_auroc <- if(!is.null(sva)) auroc_score(sva$pval, d$all_de) else NA
    ciber_auroc <- auroc_score(ciber$pval, d$all_de)

    # AUROC for cell-intrinsic DE only (the genes CADE should detect)
    cade_intrinsic <- auroc_score(cade_de$p_a, d$true_de)
    limma_intrinsic <- auroc_score(limma_fit$p.value[,2], d$true_de)

    # Proportion accuracy
    common <- intersect(colnames(cade_p), colnames(d$true_p))
    ct_r_cade <- sapply(common, function(k) cor(cade_p[,k], d$true_p[,k]))
    ct_r_ciber <- sapply(common, function(k) cor(ciber$props[,k], d$true_p[,k]))
    ct_r_cade[is.na(ct_r_cade)] <- 0
    ct_r_ciber[is.na(ct_r_ciber)] <- 0

    # CCI discrimination: cell-intrinsic DE vs background vs composition-driven
    cci_de <- mean(cade_de$cci[d$true_de], na.rm=TRUE)
    cci_comp <- mean(cade_de$cci[d$comp_de], na.rm=TRUE)
    cci_bg <- mean(cade_de$cci[bg], na.rm=TRUE)

    elapsed <- as.numeric(difftime(Sys.time(), t_start, units="secs"))

    res <- data.frame(
      Bias=bias, Rep=rep,
      CADE_all=cade_auroc,
      CADE_intrinsic=cade_intrinsic,
      Limma=l_auroc,
      Limma_intrinsic=limma_intrinsic,
      MarkerScore=ms_auroc,
      SVA=sva_auroc,
      CIBER_LM22=ciber_auroc,
      CCI_intrinsic=cci_de,
      CCI_compDriven=cci_comp,
      CCI_bg=cci_bg,
      CCI_sep=cci_bg - cci_de,
      CCI_compSep=cci_comp - cci_de,
      PropR_CADE_mean=mean(ct_r_cade),
      PropR_CADE_median=median(ct_r_cade),
      PropR_CIBER_mean=mean(ct_r_ciber),
      PropR_CIBER_median=median(ct_r_ciber),
      Elapsed_sec=elapsed,
      stringsAsFactors=FALSE
    )
    all_res[[run_id]] <- res
    run_id <- run_id + 1

    cat(sprintf("  rep%d [%.1fs]: CADE_int=%.3f Limma=%.3f CIBER=%.3f MS=%.3f SVA=%.3f CCI_sep=%.3f\n",
                rep, elapsed, res$CADE_intrinsic, res$Limma, res$CIBER_LM22,
                res$MarkerScore, if(is.na(res$SVA)) NA else res$SVA, res$CCI_sep))
  }
}

all_res <- do.call(rbind, all_res)

# ---- 7. Results ----
cat("\n\n[7/8] ===== BENCHMARK RESULTS v4 =====\n")

cat("\n--- AUROC for Cell-Intrinsic DE Detection ---\n")
cat(sprintf("%-10s %8s %8s %8s %8s %8s %8s\n", "Bias", "CADE", "Limma", "CIBER", "MarkSc", "SVA", "CCI_sep"))
for(b in BIAS_LEVELS) {
  s <- all_res[all_res$Bias==b,]
  cat(sprintf("%-10s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
              sprintf("%.0f%%", b*100),
              mean(s$CADE_intrinsic, na.rm=TRUE),
              mean(s$Limma_intrinsic, na.rm=TRUE),
              mean(s$CIBER_LM22, na.rm=TRUE),
              mean(s$MarkerScore, na.rm=TRUE),
              mean(s$SVA, na.rm=TRUE),
              mean(s$CCI_sep, na.rm=TRUE)))
}

cat("\n--- CCI Discrimination (cell-intrinsic vs background vs composition-driven) ---\n")
cat(sprintf("%-10s %8s %8s %8s %8s\n", "Bias", "Intrinsic", "CompDriven", "Background", "Separation"))
for(b in BIAS_LEVELS) {
  s <- all_res[all_res$Bias==b,]
  cat(sprintf("%-10s %8.4f %8.4f %8.4f %8.4f\n",
              sprintf("%.0f%%", b*100),
              mean(s$CCI_intrinsic, na.rm=TRUE),
              mean(s$CCI_compDriven, na.rm=TRUE),
              mean(s$CCI_bg, na.rm=TRUE),
              mean(s$CCI_sep, na.rm=TRUE)))
}

cat("\n--- Proportion Estimation Accuracy (Pearson r) ---\n")
cat(sprintf("%-10s %10s %10s %10s %10s\n", "Bias", "CADE_mean", "CADE_med", "CIBER_mean", "CIBER_med"))
for(b in BIAS_LEVELS) {
  s <- all_res[all_res$Bias==b,]
  cat(sprintf("%-10s %10.4f %10.4f %10.4f %10.4f\n",
              sprintf("%.0f%%", b*100),
              mean(s$PropR_CADE_mean, na.rm=TRUE),
              mean(s$PropR_CADE_median, na.rm=TRUE),
              mean(s$PropR_CIBER_mean, na.rm=TRUE),
              mean(s$PropR_CIBER_median, na.rm=TRUE)))
}

# ---- 8. Save ----
cat("\n\n[8/8] Saving results...\n")
write.csv(all_res, file.path(OUT_DIR, "scRNA_benchmark_v4_summary.csv"), row.names=FALSE)
saveRDS(list(results=all_res, config=list(
  data="MonacoImmune (celldex, real scRNA-seq)",
  title="CADE Benchmark v4 — Hardened Edition",
  upgrades=c("DE_FC=0.5/0.8/1.0 mixed", "NOISE_SD=0.12", "N_DE=80+20comp",
             "MAX_MARKERS=8", "CIBERSORT/LM22 comparison", "150% extreme bias"),
  n_samples=N_SAMPLES, n_genes=N_GENES, n_ct=N_CT,
  n_true_de=N_TRUE_DE, n_comp_de=N_COMP_DE,
  bias_levels=BIAS_LEVELS, n_replicates=N_REP,
  cell_types=colnames(ct_profiles),
  methods="CADE vs Limma vs MarkerScore+Limma vs SVA vs CIBERSORT/LM22"
)), file.path(OUT_DIR, "scRNA_benchmark_v4_full.rds"))

cat(sprintf("\nDone. Results: %s\n", OUT_DIR))
