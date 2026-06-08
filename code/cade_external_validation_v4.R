#!/usr/bin/env Rscript
# CADE External Validation v4
# Downloads sepsis dataset (GSE28750, GPL570) + maps probes→genes
# Uses INLINE CADE functions
library(limma)
library(quadprog)
library(GEOquery)
library(pROC)

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

# ── Runtime options ──
# Usage:
#   Rscript cade_external_validation_v4.R
#   Rscript cade_external_validation_v4.R --mode legacy
#   Rscript cade_external_validation_v4.R --mode stabilized
args <- commandArgs(trailingOnly = TRUE)
mode <- "stabilized"
if("--mode" %in% args) {
  idx <- which(args == "--mode")
  if(length(idx) > 0 && idx[1] < length(args)) {
    mode_in <- tolower(args[idx[1] + 1])
    if(mode_in %in% c("legacy", "stabilized")) {
      mode <- mode_in
    } else {
      stop("Invalid --mode. Use 'legacy' or 'stabilized'.")
    }
  }
}

# ── CADE core functions (inline) ──

estimate_proportions_cade <- function(expr_mat, marker_list, max_iter=50,
                                       tol=1e-6, min_markers=3) {
  n_ct <- length(marker_list)
  ct_names <- names(marker_list)
  n_samples <- ncol(expr_mat)

  marker_present <- lapply(marker_list, function(m) intersect(m, rownames(expr_mat)))
  n_markers_present <- sapply(marker_present, length)
  if(any(n_markers_present < min_markers)) {
    stop(sprintf("Cell types with insufficient markers: %s",
         paste(ct_names[n_markers_present < min_markers], collapse=", ")))
  }

  init_scores <- matrix(NA, nrow=n_ct, ncol=n_samples)
  rownames(init_scores) <- ct_names
  for(ct in ct_names) {
    init_scores[ct, ] <- colMeans(expr_mat[marker_present[[ct]], , drop=FALSE])
  }
  prop_mat <- apply(init_scores, 2, function(x) {
    e_x <- exp(x - max(x))
    e_x / sum(e_x)
  })
  rownames(prop_mat) <- ct_names

  all_markers <- unique(unlist(marker_present))
  n_markers_total <- length(all_markers)
  marker_ct_map <- matrix(0, nrow=n_ct, ncol=n_markers_total,
                          dimnames=list(ct_names, all_markers))
  for(ct in ct_names) {
    marker_ct_map[ct, marker_present[[ct]]] <- 1
  }

  for(iter in 1:max_iter) {
    ct_profiles <- matrix(0, nrow=n_markers_total, ncol=n_ct,
                          dimnames=list(all_markers, ct_names))
    for(g in all_markers) {
      y <- expr_mat[g, ]
      for(ct in ct_names) {
        if(marker_ct_map[ct, g] == 1) {
          X <- prop_mat[ct, ]
          fit <- lm(y ~ X)
          ct_profiles[g, ct] <- coef(fit)[2]
        }
      }
    }
    ct_profiles[is.na(ct_profiles)] <- 0

    prop_old <- prop_mat
    for(i in 1:n_samples) {
      y_i <- expr_mat[all_markers, i]
      Dmat <- t(ct_profiles) %*% ct_profiles
      dvec <- t(ct_profiles) %*% y_i
      Amat <- cbind(rep(1, n_ct), diag(n_ct))
      bvec <- c(1, rep(0, n_ct))
      qp <- solve.QP(Dmat, dvec, Amat, bvec, meq=1)
      prop_mat[, i] <- qp$solution
    }

    delta <- max(abs(prop_mat - prop_old))
    if(delta < tol) break
  }

  list(proportions=prop_mat, iterations=iter, final_delta=delta)
}

cade_de_analysis <- function(expr_mat, group, prop_mat, top_cts=3,
                             min_abs_logfc=0.1, se_floor_mult=1.0,
                             cci_variant=c("legacy", "stabilized")) {
  cci_variant <- match.arg(cci_variant)
  group_binary <- as.numeric(group)
  ct_vars <- apply(prop_mat, 1, var)
  ct_use <- names(sort(ct_vars, decreasing=TRUE)[1:min(top_cts, nrow(prop_mat))])

  design_adj <- model.matrix(~ group_binary + t(prop_mat[ct_use, , drop=FALSE]))
  design_unadj <- model.matrix(~ group_binary)

  fit_adj <- lmFit(expr_mat, design_adj)
  fit_adj <- eBayes(fit_adj, trend=TRUE)
  de_adj <- topTable(fit_adj, coef=2, number=Inf, adjust.method="BH")
  de_adj$Gene <- rownames(de_adj)

  fit_unadj <- lmFit(expr_mat, design_unadj)
  fit_unadj <- eBayes(fit_unadj, trend=TRUE)
  de_unadj <- topTable(fit_unadj, coef=2, number=Inf, adjust.method="BH")
  de_unadj$Gene <- rownames(de_unadj)

  results <- merge(
    de_unadj[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")],
    de_adj[, c("Gene", "logFC", "P.Value", "adj.P.Val")],
    by="Gene", suffixes=c(".unadj", ".adj")
  )
  results$CCI_legacy <- with(results,
    ifelse(abs(logFC.unadj) > min_abs_logfc,
           abs(logFC.unadj - logFC.adj) / abs(logFC.unadj), NA)
  )
  results$CCI_legacy <- pmax(0, pmin(1, results$CCI_legacy))

  results$SE_unadj <- with(results, ifelse(
    is.na(t) | abs(t) < 1e-8, NA, abs(logFC.unadj / t)
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
  results$CCI <- if(cci_variant == "legacy") results$CCI_legacy else results$CCI_stabilized
  results
}

# ── Probe-to-gene mapping for GPL570 ──

map_gpl570_probes <- function(expr, eset) {
  fd <- fData(eset)
  if(!"Gene Symbol" %in% colnames(fd)) {
    cat("  No 'Gene Symbol' column in feature data\n")
    return(NULL)
  }

  symbols <- as.character(fd[rownames(expr), "Gene Symbol"])
  names(symbols) <- rownames(expr)

  # Remove probes without gene symbol
  has_symbol <- symbols != "" & !is.na(symbols)
  cat(sprintf("  Probes with gene symbol: %d/%d\n", sum(has_symbol), length(has_symbol)))
  expr <- expr[has_symbol, ]
  symbols <- symbols[has_symbol]

  # Handle multi-gene probes (take first symbol after /// split)
  symbols <- sapply(strsplit(symbols, " /// "), `[`, 1)
  names(symbols) <- rownames(expr)

  # Collapse multiple probes per gene: keep probe with highest mean expression
  dup_genes <- names(which(table(symbols) > 1))
  cat(sprintf("  Genes with multiple probes: %d\n", length(dup_genes)))

  # For each gene, keep the probe with highest mean expression
  expr_gene <- matrix(NA, nrow=length(unique(symbols)), ncol=ncol(expr))
  rownames(expr_gene) <- unique(symbols)
  colnames(expr_gene) <- colnames(expr)

  for(g in unique(symbols)) {
    probes <- rownames(expr)[symbols == g]
    if(length(probes) == 1) {
      expr_gene[g, ] <- expr[probes, ]
    } else {
      # Keep probe with highest mean expression
      best <- probes[which.max(rowMeans(expr[probes, , drop=FALSE]))]
      expr_gene[g, ] <- expr[best, ]
    }
  }

  cat(sprintf("  Mapped to genes: %d genes x %d samples\n", nrow(expr_gene), ncol(expr_gene)))
  expr_gene
}

OUT_DIR <- file.path(PROJECT_ROOT, "analysis_output", "CADE", "validation")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

cat(sprintf("=== CADE External Validation v4 (mode=%s) ===\n\n", mode))

# ══════════════════════════════════════════════════════════════════════════════
# Marker genes (same as FHL analysis for consistency)
# ══════════════════════════════════════════════════════════════════════════════

MARKER_LIST <- list(
  CD8_Tcells = c("CD8A", "CD8B", "PRF1", "GZMB", "GZMA", "GNLY", "NKG7", "CD3E", "CD3D"),
  CD4_Tcells = c("CD4", "CD3E", "CD3D", "IL7R", "CCR7", "LEF1", "TCF7", "SELL"),
  NK_cells = c("NKG7", "GNLY", "PRF1", "GZMB", "KLRB1", "KLRD1", "KLRF1", "NCR1", "CD160"),
  B_cells = c("CD19", "CD79A", "CD79B", "MS4A1", "PAX5", "BLK", "BANK1", "CD22"),
  Monocytes = c("CD14", "FCGR3A", "CSF1R", "ITGAM", "LYZ", "S100A8", "S100A9", "VCAN"),
  Neutrophils = c("FCGR3B", "CXCR2", "CXCL8", "CSF3R", "MMP8", "MMP9", "ELANE", "MPO")
)

# Sepsis-relevant genes for CCI analysis
SEPSIS_GENES <- c("IL1B", "TNF", "IL6", "CXCL8", "IL10", "HMOX1",
                  "S100A8", "S100A9", "CD14", "TLR2", "TLR4",
                  "NFKB1", "STAT3", "HIF1A", "MPO", "ELANE",
                  "MMP8", "MMP9", "FCGR3B", "CSF3R")

# ══════════════════════════════════════════════════════════════════════════════
# Download GSE28750 (Sepsis whole blood, GPL570)
# ══════════════════════════════════════════════════════════════════════════════

cat("Downloading GSE28750 (sepsis whole blood, GPL570)...\n")

gse <- tryCatch({
  getGEO("GSE28750", GSEMatrix=TRUE, getGPL=TRUE)
}, error=function(e) {
  cat(sprintf("Download failed: %s\n", e$message))
  return(NULL)
})

if(is.null(gse) || length(gse) == 0) {
  stop("Could not download GSE28750")
}

eset <- gse[[1]]
expr_probe <- exprs(eset)
pheno <- pData(eset)

cat(sprintf("Loaded: %d probes x %d samples\n", nrow(expr_probe), ncol(expr_probe)))

# ── Probe → Gene mapping ──
cat("Mapping probes to genes (GPL570)...\n")
expr_gene <- map_gpl570_probes(expr_probe, eset)

# ── Log2 transform ──
if(max(expr_gene, na.rm=TRUE) > 15) {
  cat("Applying log2 transform...\n")
  expr_gene <- log2(expr_gene + 1)
}
cat(sprintf("Expression range: [%.2f, %.2f]\n", min(expr_gene), max(expr_gene)))

# ── Group assignment ──
cat("\nPhenotype columns:\n")
for(col in colnames(pheno)) {
  cat(sprintf("  %s: %s\n", col, paste(unique(pheno[[col]]), collapse=" | ")))
}

# Use characteristics_ch1.1 for group (health status: sepsis / healthy / post_surgical)
vals <- as.character(pheno[["characteristics_ch1.1"]])
cat(sprintf("\nGroup values: %s\n", paste(unique(vals), collapse=" | ")))

is_sepsis <- grepl("sepsis", vals, ignore.case=TRUE)
is_healthy <- grepl("healthy", vals, ignore.case=TRUE)
is_surgical <- grepl("post_surgical|surgical", vals, ignore.case=TRUE)

cat(sprintf("Sepsis: %d, Healthy: %d, Post-surgical: %d\n",
            sum(is_sepsis), sum(is_healthy), sum(is_surgical)))

# Compare sepsis vs healthy (exclude post_surgical)
keep <- is_sepsis | is_healthy
group <- is_sepsis[keep]
expr_final <- expr_gene[, keep]
colnames(expr_final) <- paste0("Sample", 1:ncol(expr_final))

cat(sprintf("Final: %d sepsis, %d healthy, %d genes\n", sum(group), sum(!group), nrow(expr_final)))

# ── Run CADE ──
cat("\n─────────────── CADE Validation ───────────────\n")

# Check markers
marker_avail <- sapply(MARKER_LIST, function(m) sum(m %in% rownames(expr_final)))
cat("Markers available (gene-level):\n")
for(nm in names(marker_avail)) {
  cat(sprintf("  %-16s: %d/%d\n", nm, marker_avail[nm], length(MARKER_LIST[[nm]])))
}

ct_keep <- names(marker_avail)[marker_avail >= 3]
cat(sprintf("Using %d cell types with >=3 markers\n", length(ct_keep)))
marker_filt <- MARKER_LIST[ct_keep]

cade_prop <- estimate_proportions_cade(expr_final, marker_filt, max_iter=50, tol=1e-6)
cat(sprintf("CADE converged: iter=%d, delta=%.2e\n",
            cade_prop$n_iter, tail(cade_prop$convergence, 1)))

cade_de <- cade_de_analysis(expr_final, group, cade_prop$proportions, top_cts=3, cci_variant=mode)
rownames(cade_de) <- cade_de$Gene

# ── Results ──

# CCI for background genes
all_markers <- unique(unlist(marker_filt))
cci_bg <- cade_de[setdiff(rownames(cade_de), all_markers), "CCI"]
cci_bg <- cci_bg[!is.na(cci_bg)]

cat(sprintf("\nCCI background (n=%d): mean=%.4f, median=%.4f, SD=%.4f\n",
            length(cci_bg), mean(cci_bg), median(cci_bg), sd(cci_bg)))
cat(sprintf("CCI percentiles: 5th=%.4f, 25th=%.4f, 75th=%.4f, 95th=%.4f\n",
            quantile(cci_bg, 0.05), quantile(cci_bg, 0.25),
            quantile(cci_bg, 0.75), quantile(cci_bg, 0.95)))

# Proportion differences
cat("\nCADE-estimated cell-type proportion differences (sepsis - healthy):\n")
for(i in 1:nrow(cade_prop$proportions)) {
  ct <- rownames(cade_prop$proportions)[i]
  d_mean <- mean(cade_prop$proportions[i, group])
  h_mean <- mean(cade_prop$proportions[i, !group])
  cat(sprintf("  %-16s: sepsis=%.4f, healthy=%.4f, Δ=%+.4f\n",
              ct, d_mean, h_mean, d_mean - h_mean))
}

# Known gene CCI
g_avail <- intersect(SEPSIS_GENES, rownames(cade_de))
if(length(g_avail) > 0) {
  gene_df <- cade_de[g_avail, c("Gene", "logFC.unadj", "logFC.adj", "CCI", "adj.P.Val.adj")]
  gene_df <- gene_df[order(gene_df$CCI), ]
  cat(sprintf("\nSepsis signaling genes (n=%d, CCI-ranked):\n", nrow(gene_df)))
  cat(sprintf("  %-12s %8s %8s %8s %8s\n", "Gene", "logFC_raw", "logFC_adj", "CCI", "FDR_adj"))
  for(i in 1:nrow(gene_df)) {
    r <- gene_df[i,]
    cat(sprintf("  %-12s %+8.2f %+8.2f %8.3f %8.4f\n",
                r$Gene, r$logFC.unadj, r$logFC.adj, r$CCI, r$adj.P.Val.adj))
  }
}

# DE counts
cat(sprintf("\nDE genes (FDR<0.05): %d (CADE) vs %d (unadj)\n",
            sum(cade_de$adj.P.Val.adj < 0.05, na.rm=TRUE),
            sum(cade_de$adj.P.Val.unadj < 0.05, na.rm=TRUE)))

# ── Save results ──
# CCI summary
cci_summary <- data.frame(
  Metric = c("Mean", "Median", "SD", "Q05", "Q25", "Q75", "Q95", "N_Genes"),
  Value = c(mean(cci_bg), median(cci_bg), sd(cci_bg),
            as.numeric(quantile(cci_bg, c(0.05, 0.25, 0.75, 0.95))),
            length(cci_bg)),
  stringsAsFactors = FALSE
)
file_suffix <- if(mode == "legacy") "" else paste0("_", mode)
write.csv(cci_summary, file.path(OUT_DIR, paste0("Table_Sepsis_GSE28750_CCI_Distribution", file_suffix, ".csv")),
          row.names=FALSE)

# Proportion differences
prop_df <- data.frame(
  CellType = rownames(cade_prop$proportions),
  Sepsis = rowMeans(cade_prop$proportions[, group, drop=FALSE]),
  Healthy = rowMeans(cade_prop$proportions[, !group, drop=FALSE]),
  Delta = rowMeans(cade_prop$proportions[, group, drop=FALSE]) -
          rowMeans(cade_prop$proportions[, !group, drop=FALSE]),
  stringsAsFactors = FALSE
)
write.csv(prop_df, file.path(OUT_DIR, "Table_Sepsis_GSE28750_Proportions.csv"),
          row.names=FALSE)

# Gene-level CCI
write.csv(gene_df, file.path(OUT_DIR, paste0("Table_Sepsis_GSE28750_GeneCCI", file_suffix, ".csv")),
          row.names=FALSE)

cat(sprintf("\nResults saved to: %s\n", OUT_DIR))
cat("Files:\n")
for(f in list.files(OUT_DIR, pattern="*.csv")) {
  cat(sprintf("  - %s\n", f))
}
cat("\n=== Done ===\n")
