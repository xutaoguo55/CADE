#!/usr/bin/env Rscript
# ============================================================
# Standalone WGCNA analysis (clean session to avoid conflicts)
# ============================================================
library(WGCNA)
library(dplyr)

set.seed(42)

# ── Package-local path configuration ──
if (exists("CODE_DIR", envir = .GlobalEnv)) {
  SCRIPT_DIR <- normalizePath(get("CODE_DIR", envir = .GlobalEnv), mustWork = TRUE)
} else {
  SCRIPT_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
}
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)

OUT_DIR <- file.path(PROJECT_ROOT, "analysis_output", "CADE")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

cat("=== WGCNA Analysis (standalone) ===\n")

# Load expression matrix
exprs_mat <- as.matrix(read.csv(
  file.path(PROJECT_ROOT, "geo_analysis_output", "GSE26050_processed_matrix.csv"),
  row.names=1, check.names=FALSE
))
cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(exprs_mat), ncol(exprs_mat)))

# Group assignment based on GEO metadata
# GSE26050: GSM639703-639713 = FHL (PRF/WT patients), GSM639714-639746 = HC (CTRL)
sample_names <- colnames(exprs_mat)
group <- ifelse(sample_names %in% c(
  "GSM639703", "GSM639704", "GSM639705", "GSM639706", "GSM639707",
  "GSM639708", "GSM639709", "GSM639710", "GSM639711", "GSM639712",
  "GSM639713"
), "FHL", "HC")
cat(sprintf("  FHL: %d, HC: %d\n", sum(group=="FHL"), sum(group=="HC")))

# Focus on most variable genes
n_genes_wgcna <- min(5000, nrow(exprs_mat))
gene_vars <- apply(exprs_mat, 1, var, na.rm=TRUE)
top_var_genes <- names(sort(gene_vars, decreasing=TRUE)[1:n_genes_wgcna])
exprs_wgcna <- exprs_mat[top_var_genes, ]
datExpr <- t(exprs_wgcna)

# Check for missing genes/samples
gsg <- goodSamplesGenes(datExpr, verbose=0)
if(!gsg$allOK) {
  datExpr <- datExpr[, gsg$goodGenes]
  cat(sprintf("  Removed %d problematic genes\n", sum(!gsg$goodGenes)))
}

# Pick soft threshold
cat("  Picking soft threshold...\n")
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=0, networkType="signed")
best_power <- sft$powerEstimate
if(is.na(best_power)) {
  best_power <- which(sft$fitIndices[, "SFT.R.sq"] > 0.8)[1]
  if(is.na(best_power)) best_power <- 6
}
cat(sprintf("  Best soft power: %d (scale-free R=%.3f)\n",
            best_power, sft$fitIndices[best_power, "SFT.R.sq"]))

# Blockwise module detection
cat("  Running blockwiseModules...\n")
enableWGCNAThreads(nThreads=2)
net <- blockwiseModules(datExpr, power=best_power,
                         TOMType="signed", minModuleSize=30,
                         reassignThreshold=0, mergeCutHeight=0.25,
                         numericLabels=TRUE, pamRespectsDendro=FALSE,
                         saveTOMs=FALSE, verbose=0, maxBlockSize=5000)

cat(sprintf("  Detected %d modules\n", length(unique(net$colors))))

# Correlate modules with FHL status
fhl_binary <- as.numeric(group == "FHL")
module_eigengenes <- net$MEs
module_trait_cor <- cor(module_eigengenes, fhl_binary, use="complete.obs")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

# Find modules significantly associated with FHL
module_stats <- data.frame(
  Module = gsub("ME", "", colnames(module_eigengenes)),
  Cor_FHL = as.vector(module_trait_cor),
  P_value = as.vector(module_trait_pval),
  Size = as.vector(table(net$colors)),
  stringsAsFactors = FALSE
)
module_stats <- module_stats[module_stats$Module != "0", ]
module_stats$FDR <- p.adjust(module_stats$P_value, method="BH")

cat("\n  Module association with FHL:\n")
for(i in 1:min(10, nrow(module_stats))) {
  r <- module_stats[i,]
  sig_mark <- if(r$FDR < 0.05) "***" else if(r$FDR < 0.10) "*" else ""
  cat(sprintf("    Module %s: r=%+.3f P=%.4f FDR=%.4f Size=%d %s\n",
              r$Module, r$Cor_FHL, r$P_value, r$FDR, r$Size, sig_mark))
}

# Extract hub genes from FHL-associated modules
sig_modules <- module_stats$Module[module_stats$FDR < 0.05]
ferroptosis_all <- c("SLC7A11", "GPX4", "ACSL4", "LPCAT3", "TFRC",
                      "SLC25A37", "FTH1", "FTL", "HMOX1", "NCOA4",
                      "KEAP1", "NFE2L2", "GCLM", "GCLC", "SLC40A1",
                      "SAT1", "ALOX15", "TXN", "TXNRD1", "SOD1", "SOD2", "CAT")

for(mod in sig_modules) {
  mod_genes <- top_var_genes[net$colors == as.numeric(mod)]
  if(length(mod_genes) > 0) {
    mod_kme <- cor(datExpr[, mod_genes, drop=FALSE],
                   module_eigengenes[, paste0("ME", mod)], use="complete.obs")
    mod_kme <- sort(abs(mod_kme[,1]), decreasing=TRUE)
    top_hubs <- names(mod_kme)[1:min(30, length(mod_kme))]

    cat(sprintf("\n  Module %s (r=%.3f, FDR=%.2e, size=%d)\n",
                mod, module_stats$Cor_FHL[module_stats$Module == mod],
                module_stats$FDR[module_stats$Module == mod], length(mod_genes)))
    cat(sprintf("  Top hubs: %s\n", paste(top_hubs[1:15], collapse=", ")))

    # Ferroptosis genes in this module
    ferr_in_mod <- intersect(mod_genes, ferroptosis_all)
    if(length(ferr_in_mod) > 0) {
      cat(sprintf("  FERROPTOSIS GENES IN MODULE: %s\n", paste(ferr_in_mod, collapse=", ")))
    }
  }
}

# Save results
write.csv(module_stats, file.path(OUT_DIR, "Table_WGCNA_Module_FHL_Correlation.csv"), row.names=FALSE)

# Save module assignments for all top genes
module_assignments <- data.frame(
  Gene = top_var_genes,
  Module = net$colors,
  stringsAsFactors = FALSE
)
write.csv(module_assignments, file.path(OUT_DIR, "Table_WGCNA_Gene_Module_Assignments.csv"), row.names=FALSE)

# Save hub genes from significant modules
hub_genes_out <- data.frame()
for(mod in sig_modules) {
  mod_genes <- top_var_genes[net$colors == as.numeric(mod)]
  if(length(mod_genes) > 0) {
    mod_kme <- cor(datExpr[, mod_genes, drop=FALSE],
                   module_eigengenes[, paste0("ME", mod)], use="complete.obs")
    mod_kme <- sort(abs(mod_kme[,1]), decreasing=TRUE)
    top_hubs <- names(mod_kme)[1:min(50, length(mod_kme))]
    hub_genes_out <- rbind(hub_genes_out, data.frame(
      Module = mod,
      Gene = top_hubs,
      kME = mod_kme[1:min(50, length(mod_kme))],
      Rank = 1:min(50, length(mod_kme)),
      stringsAsFactors = FALSE
    ))
  }
}
write.csv(hub_genes_out, file.path(OUT_DIR, "Table_WGCNA_Hub_Genes.csv"), row.names=FALSE)

cat(sprintf("\n=== WGCNA Complete: %d modules, %d significant ===\n",
            nrow(module_stats), length(sig_modules)))
