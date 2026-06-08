#!/usr/bin/env Rscript
# ============================================================
# Deeper Analysis v3: Cell deconvolution, WGCNA, ferroptosis
# scoring, immune infiltration, and TF activity inference
# All analyses use real data from GEO GSE26050
# ============================================================
library(GEOquery)
library(limma)
library(dplyr)
library(tibble)
library(ggplot2)
library(GSVA)
library(GSEABase)
library(fgsea)
library(msigdbr)

# ── Package-local path configuration ──
if (exists("CODE_DIR", envir = .GlobalEnv)) {
  SCRIPT_DIR <- normalizePath(get("CODE_DIR", envir = .GlobalEnv), mustWork = TRUE)
} else {
  SCRIPT_DIR <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
}
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)

set.seed(42)
OUT_DIR <- file.path(PROJECT_ROOT, "analysis_output", "CADE")
dir.create(OUT_DIR, showWarnings=FALSE, recursive=TRUE)

# ============================================================
# Step 1: Load and prepare data
# ============================================================
cat("=== Step 1: Loading expression data ===\n")

# Load the processed matrix from GEO re-analysis
exprs_mat <- as.matrix(read.csv(
  file.path(PROJECT_ROOT, "geo_analysis_output", "GSE26050_processed_matrix.csv"),
  row.names=1, check.names=FALSE
))
cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(exprs_mat), ncol(exprs_mat)))

# Load group assignments
gse <- getGEO("GSE26050", GSEMatrix=TRUE, getGPL=FALSE,
              destdir=file.path(PROJECT_ROOT, "geo_analysis_output"))
eset <- gse[[1]]
pdata <- pData(eset)
sample_titles <- pdata$title
group <- ifelse(grepl("CTRL", sample_titles), "HC", "FHL")
cat(sprintf("  FHL: %d, HC: %d\n", sum(group=="FHL"), sum(group=="HC")))

# Load group assignments from known GEO metadata
# GSE26050: GSM639703-639713 = FHL, GSM639714-639746 = HC
common_samples <- intersect(colnames(exprs_mat), colnames(eset))
if(length(common_samples) < ncol(exprs_mat)) {
  exprs_mat <- exprs_mat[, common_samples]
}
group <- ifelse(common_samples %in% c(
  "GSM639703", "GSM639704", "GSM639705", "GSM639706", "GSM639707",
  "GSM639708", "GSM639709", "GSM639710", "GSM639711", "GSM639712",
  "GSM639713"
), "FHL", "HC")
cat(sprintf("  Final samples: %d (FHL=%d, HC=%d)\n",
            length(common_samples), sum(group=="FHL"), sum(group=="HC")))

# ============================================================
# Step 2: Per-sample ferroptosis scoring (GSVA)
# ============================================================
cat("\n=== Step 2: Ferroptosis pathway scoring per sample ===\n")

# Define ferroptosis gene sets from curated sources
ferroptosis_sets <- list(
  Ferroptosis_Core = c("SLC7A11", "GPX4", "ACSL4", "LPCAT3", "TFRC",
                         "SLC25A37", "FTH1", "FTL", "HMOX1", "NCOA4",
                         "KEAP1", "NFE2L2", "GCLM", "GCLC", "SLC40A1",
                         "SAT1", "ALOX15"),
  Ferroptosis_Defense = c("SLC7A11", "GPX4", "NFE2L2", "GCLM", "GCLC",
                           "FTH1", "FTL", "SLC40A1", "HMOX1", "TXN",
                           "TXNRD1", "SOD1", "SOD2", "CAT", "PRDX1", "PRDX6"),
  Ferroptosis_Drivers = c("TFRC", "SLC25A37", "ACSL4", "LPCAT3", "ALOX5",
                           "ALOX15", "NCOA4", "SAT1", "LOX", "VDAC2", "VDAC3"),
  Iron_Homeostasis = c("FTH1", "FTL", "TFRC", "SLC40A1", "SLC25A37",
                        "HMOX1", "IREB2", "HAMP", "TFR2", "STEAP3"),
  GSH_Metabolism = c("SLC7A11", "GCLM", "GCLC", "GPX4", "GSR",
                      "GSS", "SLC3A2", "GGT1")
)

# Run GSVA (v2.x API with parameter objects)
cat("  Running GSVA for ferroptosis gene sets...\n")
gsva_param <- gsvaParam(exprs_mat, ferroptosis_sets, minSize=3, maxSize=200)
gsva_ferroptosis <- gsva(gsva_param, verbose=FALSE)

# Test differential enrichment between FHL vs HC
gsva_results <- data.frame()
for(set_name in rownames(gsva_ferroptosis)) {
  fhl_scores <- gsva_ferroptosis[set_name, group == "FHL"]
  hc_scores <- gsva_ferroptosis[set_name, group == "HC"]
  tt <- t.test(fhl_scores, hc_scores)
  gsva_results <- rbind(gsva_results, data.frame(
    GeneSet = set_name,
    Mean_FHL = mean(fhl_scores), Mean_HC = mean(hc_scores),
    Delta = mean(fhl_scores) - mean(hc_scores),
    t_stat = tt$statistic, P_value = tt$p.value,
    stringsAsFactors = FALSE
  ))
}
gsva_results$FDR <- p.adjust(gsva_results$P_value, method="BH")

cat("  Per-sample ferroptosis enrichment (GSVA):\n")
for(i in 1:nrow(gsva_results)) {
  r <- gsva_results[i,]
  cat(sprintf("    %-25s Δ=%+6.3f P=%.4f FDR=%.4f\n",
              r$GeneSet, r$Delta, r$P_value, r$FDR))
}

write.csv(gsva_results, file.path(OUT_DIR, "Table_Ferroptosis_GSVA_Scores.csv"), row.names=FALSE)

# ============================================================
# Step 3: Immune cell type scoring using curated marker genes
# ============================================================
cat("\n=== Step 3: Immune cell type scoring ===\n")

# Curated immune cell marker gene sets (from LM22 / CIBERSORT + literature)
immune_sigs <- list(
  CD8_Tcells = c("CD8A", "CD8B", "PRF1", "GZMB", "GZMA", "GNLY", "NKG7", "CD3E", "CD3D"),
  CD4_Tcells = c("CD4", "CD3E", "CD3D", "IL7R", "CCR7", "LEF1", "TCF7", "SELL"),
  NK_cells = c("NKG7", "GNLY", "PRF1", "GZMB", "KLRB1", "KLRD1", "KLRF1", "NCR1", "CD160"),
  B_cells = c("CD19", "CD79A", "CD79B", "MS4A1", "PAX5", "BLK", "BANK1", "CD22"),
  Monocytes = c("CD14", "FCGR3A", "CSF1R", "ITGAM", "LYZ", "S100A8", "S100A9", "VCAN"),
  Macrophages = c("CD68", "CD163", "MRC1", "MSR1", "CSF1R", "ITGAM", "TLR2", "TLR4"),
  Neutrophils = c("FCGR3B", "CXCR2", "CXCL8", "CSF3R", "MMP8", "MMP9", "ELANE", "MPO"),
  Erythrocytes = c("HBB", "HBA1", "HBA2", "GYPA", "ALAS2", "CA1", "SLC4A1", "AHSP"),
  Hematopoietic_Stem = c("CD34", "KIT", "FLT3", "GATA2", "HOXA9", "MEIS1", "PROM1")
)

# Score each sample using ssGSEA-like approach (mean expression of marker genes)
immune_scores <- matrix(NA, nrow=length(immune_sigs), ncol=ncol(exprs_mat))
rownames(immune_scores) <- names(immune_sigs)
colnames(immune_scores) <- colnames(exprs_mat)

for(ct in names(immune_sigs)) {
  markers_present <- intersect(immune_sigs[[ct]], rownames(exprs_mat))
  if(length(markers_present) >= 3) {
    immune_scores[ct, ] <- colMeans(exprs_mat[markers_present, , drop=FALSE])
  } else {
    cat(sprintf("  WARNING: %s has only %d markers present\n", ct, length(markers_present)))
  }
}

# Differential testing
immune_results <- data.frame()
for(ct in rownames(immune_scores)) {
  fhl_val <- immune_scores[ct, group == "FHL"]
  hc_val <- immune_scores[ct, group == "HC"]
  if(length(fhl_val) > 1 && length(hc_val) > 1) {
    tt <- t.test(fhl_val, hc_val)
    immune_results <- rbind(immune_results, data.frame(
      CellType = ct, Mean_FHL = mean(fhl_val), Mean_HC = mean(hc_val),
      Delta = mean(fhl_val) - mean(hc_val),
      P_value = tt$p.value, stringsAsFactors = FALSE
    ))
  }
}
immune_results$FDR <- p.adjust(immune_results$P_value, method="BH")

cat("  Immune cell type enrichment in FHL vs HC:\n")
for(i in 1:nrow(immune_results)) {
  r <- immune_results[i,]
  sig_mark <- if(r$FDR < 0.05) "***" else if(r$FDR < 0.10) "*" else ""
  cat(sprintf("    %-25s Δ=%+7.3f FDR=%.4f %s\n",
              r$CellType, r$Delta, r$FDR, sig_mark))
}

write.csv(immune_results, file.path(OUT_DIR, "Table_Immune_CellType_Scores.csv"), row.names=FALSE)

# ============================================================
# Step 4: WGCNA Module Integration (from standalone analysis)
# ============================================================
cat("\n=== Step 4: WGCNA Co-expression Network Results ===\n")

# Read WGCNA results from standalone run
wgcna_mod_file <- file.path(OUT_DIR, "Table_WGCNA_Module_FHL_Correlation.csv")
module_stats <- NULL
if(file.exists(wgcna_mod_file)) {
  module_stats <- read.csv(wgcna_mod_file, stringsAsFactors=FALSE)
  cat(sprintf("  Loaded %d WGCNA modules from standalone analysis\n", nrow(module_stats)))

  sig_modules <- module_stats$Module[module_stats$FDR < 0.05]
  cat(sprintf("  FHL-associated modules (FDR<0.05): %d\n", length(sig_modules)))
  for(i in 1:nrow(module_stats)) {
    r <- module_stats[i,]
    sig_mark <- if(r$FDR < 0.05) "***" else if(r$FDR < 0.10) "*" else ""
    cat(sprintf("    Module %s: r=%+.3f FDR=%.2e Size=%d %s\n",
                r$Module, r$Cor_FHL, r$FDR, r$Size, sig_mark))
  }

  # Read hub genes
  hub_file <- file.path(OUT_DIR, "Table_WGCNA_Hub_Genes.csv")
  if(file.exists(hub_file)) {
    hub_genes <- read.csv(hub_file, stringsAsFactors=FALSE)
    for(mod in sig_modules) {
      mod_hubs <- hub_genes[hub_genes$Module == mod, ]
      if(nrow(mod_hubs) > 0) {
        cat(sprintf("    Module %s top hubs: %s\n", mod,
                    paste(head(mod_hubs$Gene, 10), collapse=", ")))
      }
    }
  }
} else {
  cat("  WGCNA results not found. Run wgcna_standalone.R first.\n")
}

# ============================================================
# Step 5: TF Activity Inference (DoRothEA-like)
# ============================================================
cat("\n=== Step 5: TF activity inference ===\n")

# Curated TF-target relationships from DoRothEA + literature
# Using TFs relevant to ferroptosis and HLH
tf_regulons <- list(
  STAT3 = c("FTH1", "FTL", "SLC7A11", "HMOX1", "IL6", "SOCS3", "FOS", "JUNB",
             "BCL2", "CCND1", "MYC", "VEGF", "HIF1A"),
  NFKB1 = c("IL1B", "IL6", "TNF", "CXCL8", "PTGS2", "CCL2", "ICAM1", "VCAM1",
             "BCL2L1", "XIAP", "SOD2", "MMP9"),
  NFE2L2 = c("SLC7A11", "GCLM", "GCLC", "HMOX1", "FTH1", "FTL", "GPX4",
              "NQO1", "TXNRD1", "PRDX1", "SOD1", "CAT", "AKR1C1", "ME1"),
  STAT1 = c("IRF1", "IRF7", "IFNG", "MX1", "MX2", "OAS1", "OAS2", "ISG15",
             "IFIT1", "IFIT3", "GBP1", "GBP2", "CIITA", "TAP1"),
  HIF1A = c("SLC7A11", "TFRC", "HMOX1", "VEGF", "LDHA", "HK2", "PKM",
             "SLC2A1", "BNIP3", "BNIP3L", "EPO", "PGK1", "ALDOA"),
  IRF1 = c("IFNG", "CXCL8", "GBP1", "OAS1", "TAP1", "PSMB9", "CASP1",
            "STAT1", "CIITA", "IL15"),
  TP53 = c("CDKN1A", "BAX", "BBC3", "MDM2", "GADD45A", "RRM2B", "FDXR",
            "SESN1", "SESN2", "TIGAR", "SLC7A11", "GPX4")
)

# Calculate TF activity score = mean expression of target genes
tf_scores <- matrix(NA, nrow=length(tf_regulons), ncol=ncol(exprs_mat))
rownames(tf_scores) <- names(tf_regulons)
colnames(tf_scores) <- colnames(exprs_mat)

for(tf in names(tf_regulons)) {
  targets_present <- intersect(tf_regulons[[tf]], rownames(exprs_mat))
  if(length(targets_present) >= 5) {
    tf_scores[tf, ] <- colMeans(exprs_mat[targets_present, , drop=FALSE])
  }
}

# Differential TF activity
tf_results <- data.frame()
for(tf in rownames(tf_scores)) {
  fhl_val <- tf_scores[tf, group == "FHL"]
  hc_val <- tf_scores[tf, group == "HC"]
  if(length(fhl_val) > 1 && length(hc_val) > 1) {
    tt <- t.test(fhl_val, hc_val)
    tf_results <- rbind(tf_results, data.frame(
      TF = tf, Mean_FHL = mean(fhl_val), Mean_HC = mean(hc_val),
      Delta = mean(fhl_val) - mean(hc_val),
      P_value = tt$p.value, stringsAsFactors = FALSE
    ))
  }
}
tf_results$FDR <- p.adjust(tf_results$P_value, method="BH")

cat("  TF activity inference (target gene mean expression):\n")
for(i in 1:nrow(tf_results)) {
  r <- tf_results[i,]
  direction <- if(r$Delta > 0) "ACTIVATED" else "REPRESSED"
  sig_mark <- if(r$FDR < 0.05) "***" else if(r$FDR < 0.10) "*" else ""
  cat(sprintf("    %-10s Δ=%+7.3f FDR=%.4f [%s] %s\n",
              r$TF, r$Delta, r$FDR, direction, sig_mark))
}

write.csv(tf_results, file.path(OUT_DIR, "Table_TF_Activity_Inference.csv"), row.names=FALSE)

# ============================================================
# Step 6: Published HLH/ferroptosis gene signature validation
# ============================================================
cat("\n=== Step 6: External signature validation ===\n")

# Validate key published gene signatures in our dataset
published_sigs <- list(
  HLH_Inflammatory = c("IFNG", "CXCL9", "CXCL10", "IL6", "TNF", "IL1B",
                         "IL18", "CCL2", "CCL3", "CXCL8"),
  Ferroptosis_Sensitivity = c("SLC7A11", "GPX4", "NFE2L2", "FTH1",
                               "ACSL4", "TFRC", "GCLM", "HMOX1"),
  JAK_STAT_Hyperactivation = c("STAT1", "STAT3", "JAK2", "SOCS1", "SOCS3",
                                "IRF1", "IFNG", "IFNGR1", "IL6ST"),
  NFKB_Activation = c("NFKB1", "RELA", "IKBKB", "TNF", "IL1B", "IL6",
                       "CXCL8", "PTGS2", "CCL2", "ICAM1")
)

# Score each sample for each signature
sig_scores <- matrix(NA, nrow=length(published_sigs), ncol=ncol(exprs_mat))
rownames(sig_scores) <- names(published_sigs)
colnames(sig_scores) <- colnames(exprs_mat)

for(sig_name in names(published_sigs)) {
  genes_present <- intersect(published_sigs[[sig_name]], rownames(exprs_mat))
  if(length(genes_present) >= 3) {
    sig_scores[sig_name, ] <- colMeans(exprs_mat[genes_present, , drop=FALSE])
  }
}

sig_results <- data.frame()
for(sig_name in rownames(sig_scores)) {
  fhl_val <- sig_scores[sig_name, group == "FHL"]
  hc_val <- sig_scores[sig_name, group == "HC"]
  if(length(fhl_val) > 1 && length(hc_val) > 1) {
    tt <- t.test(fhl_val, hc_val)
    sig_results <- rbind(sig_results, data.frame(
      Signature = sig_name, Mean_FHL = mean(fhl_val), Mean_HC = mean(hc_val),
      Delta = mean(fhl_val) - mean(hc_val),
      P_value = tt$p.value, stringsAsFactors = FALSE
    ))
  }
}
sig_results$FDR <- p.adjust(sig_results$P_value, method="BH")

cat("  Published signature validation in GSE26050:\n")
for(i in 1:nrow(sig_results)) {
  r <- sig_results[i,]
  cat(sprintf("    %-30s Δ=%+7.3f FDR=%.2e\n", r$Signature, r$Delta, r$FDR))
}

write.csv(sig_results, file.path(OUT_DIR, "Table_External_Signature_Validation.csv"), row.names=FALSE)

# ============================================================
# Step 7: Ferroptosis gene - immune cell type correlation
# ============================================================
cat("\n=== Step 7: Ferroptosis-immune cell correlations ===\n")

# Correlate key ferroptosis genes with immune cell type scores
key_ferr_genes <- c("SLC7A11", "SLC25A37", "GCLM", "GPX4", "FTH1", "FTL",
                     "TFRC", "HMOX1", "SLC40A1", "NCOA4", "ACSL4", "NFE2L2",
                     "STAT3", "JAK2", "IL1B", "TNF", "IL6", "CXCL8", "IFNG")

immune_cell_types <- rownames(immune_scores)

ferr_immune_cor <- matrix(NA, nrow=length(key_ferr_genes), ncol=length(immune_cell_types))
rownames(ferr_immune_cor) <- key_ferr_genes
colnames(ferr_immune_cor) <- immune_cell_types

for(g in key_ferr_genes) {
  if(g %in% rownames(exprs_mat)) {
    for(ct in immune_cell_types) {
      ct_complete <- complete.cases(exprs_mat[g, ], immune_scores[ct, ])
      if(sum(ct_complete) >= 10) {
        ferr_immune_cor[g, ct] <- cor(exprs_mat[g, ct_complete],
                                       immune_scores[ct, ct_complete],
                                       method="spearman")
      }
    }
  }
}

cat("  Top ferroptosis gene correlations with immune cell types:\n")
for(g in c("SLC7A11", "GPX4", "FTH1", "STAT3", "IL1B")) {
  if(g %in% rownames(ferr_immune_cor)) {
    cors <- sort(abs(ferr_immune_cor[g, ]), decreasing=TRUE)
    top3 <- names(cors)[1:3]
    cat(sprintf("    %s:", g))
    for(ct in top3) {
      cat(sprintf(" %s(ρ=%.2f)", ct, ferr_immune_cor[g, ct]))
    }
    cat("\n")
  }
}

write.csv(ferr_immune_cor, file.path(OUT_DIR, "Table_Ferroptosis_Immune_Correlations.csv"))

# ============================================================
# Step 8: Generate publication-quality figures
# ============================================================
cat("\n=== Step 8: Generating figures ===\n")

# Set up plotting parameters
plt_theme <- theme_bw(base_size=9) + theme(
  plot.title = element_text(size=11, face="bold"),
  axis.title = element_text(size=9),
  axis.text = element_text(size=7),
  legend.text = element_text(size=7),
  legend.title = element_text(size=8)
)

# --- Figure A: Ferroptosis GSVA score boxplot ---
fig_a_data <- data.frame(
  Sample = colnames(exprs_mat),
  Group = group,
  Ferroptosis_Core = as.numeric(gsva_ferroptosis["Ferroptosis_Core", ]),
  Ferroptosis_Defense = as.numeric(gsva_ferroptosis["Ferroptosis_Defense", ]),
  Iron_Homeostasis = as.numeric(gsva_ferroptosis["Iron_Homeostasis", ])
)

p1 <- ggplot(reshape2::melt(fig_a_data, id.vars=c("Sample", "Group"),
                              variable.name="Pathway", value.name="GSVA_Score"),
             aes(x=Group, y=GSVA_Score, fill=Group)) +
  geom_boxplot(outlier.size=0.5, alpha=0.8) +
  geom_jitter(width=0.15, size=0.8, alpha=0.5) +
  facet_wrap(~Pathway, ncol=3) +
  scale_fill_manual(values=c("HC"="#3498DB", "FHL"="#E74C3C")) +
  plt_theme + labs(title="Ferroptosis Pathway Scores (GSVA)", y="GSVA Enrichment Score")

ggsave(file.path(OUT_DIR, "Figure_Ferroptosis_GSVA_Scores.png"), p1,
       width=10, height=4, dpi=300)
ggsave(file.path(OUT_DIR, "Figure_Ferroptosis_GSVA_Scores.tif"), p1,
       width=10, height=4, dpi=300)
cat("  Saved ferroptosis GSVA figure\n")

# --- Figure B: Immune cell type landscape ---
immune_long <- data.frame(
  Sample = colnames(exprs_mat),
  Group = group
)
for(ct in rownames(immune_scores)) {
  immune_long[[ct]] <- immune_scores[ct, ]
}
immune_long_m <- reshape2::melt(immune_long, id.vars=c("Sample", "Group"),
                                 variable.name="CellType", value.name="Score")

p2 <- ggplot(immune_long_m, aes(x=CellType, y=Score, fill=Group)) +
  geom_boxplot(outlier.size=0.3, alpha=0.8, lwd=0.3) +
  scale_fill_manual(values=c("HC"="#3498DB", "FHL"="#E74C3C")) +
  plt_theme + coord_flip() +
  labs(title="Immune Cell Type Signature Scores", x="", y="Mean Marker Expression (log2)")

ggsave(file.path(OUT_DIR, "Figure_Immune_CellType_Boxplot.png"), p2,
       width=8, height=6, dpi=300)
ggsave(file.path(OUT_DIR, "Figure_Immune_CellType_Boxplot.tif"), p2,
       width=8, height=6, dpi=300)
cat("  Saved immune cell type figure\n")

# --- Figure C: TF Activity Heatmap ---
tf_scores_clean <- tf_scores[complete.cases(tf_scores), ]
tf_scaled <- t(scale(t(tf_scores_clean)))

# Prepare annotation
anno_col <- data.frame(Group=group, row.names=colnames(tf_scaled))
anno_colors <- list(Group=c("FHL"="#E74C3C", "HC"="#3498DB"))

library(pheatmap)
p3 <- pheatmap(tf_scaled, annotation_col=anno_col, annotation_colors=anno_colors,
               cluster_rows=TRUE, cluster_cols=TRUE,
               color=colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
               main="TF Activity Inference (Target Gene Expression)",
               fontsize=8, fontsize_row=9, fontsize_col=6,
               border_color=NA, silent=TRUE)

# Save pheatmap
png(file.path(OUT_DIR, "Figure_TF_Activity_Heatmap.png"), width=10, height=6,
    units="in", res=300)
grid::grid.draw(p3$gtable)
dev.off()
tiff(file.path(OUT_DIR, "Figure_TF_Activity_Heatmap.tif"), width=10, height=6,
     units="in", res=300)
grid::grid.draw(p3$gtable)
dev.off()
cat("  Saved TF activity heatmap\n")

# --- Figure D: WGCNA Module-Trait Correlation ---
if(!is.null(module_stats) && nrow(module_stats) > 0) {
  sig_mods <- module_stats[module_stats$FDR < 0.05, ]
} else {
  sig_mods <- data.frame()
}
if(nrow(sig_mods) > 0) {
  p4_data <- sig_mods[order(sig_mods$Cor_FHL), ]
  p4_data$Module <- factor(p4_data$Module, levels=p4_data$Module)

  p4 <- ggplot(p4_data, aes(x=Module, y=Cor_FHL, fill=Cor_FHL)) +
    geom_bar(stat="identity", color="black", lwd=0.3) +
    geom_text(aes(label=sprintf("P=%.1e\nn=%d", FDR, Size)),
              hjust=ifelse(p4_data$Cor_FHL > 0, -0.1, 1.1), size=2.5) +
    scale_fill_gradient2(low="#3498DB", mid="white", high="#E74C3C") +
    coord_flip() +
    plt_theme + labs(title="WGCNA Modules Associated with FHL Status",
                     y="Correlation with FHL (r)", x="Module")

  ggsave(file.path(OUT_DIR, "Figure_WGCNA_FHL_Modules.png"), p4,
         width=8, height=5, dpi=300)
  ggsave(file.path(OUT_DIR, "Figure_WGCNA_FHL_Modules.tif"), p4,
         width=8, height=5, dpi=300)
  cat("  Saved WGCNA module figure\n")
} else {
  cat("  Skipped WGCNA module figure because no WGCNA module statistics were available.\n")
}

# --- Figure E: Ferroptosis-Immune Correlation Heatmap ---
# Focus on significant relationships
ferr_immune_clean <- ferr_immune_cor
ferr_immune_clean[is.na(ferr_immune_clean)] <- 0

p5 <- pheatmap(ferr_immune_clean,
               color=colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
               main="Ferroptosis Gene vs Immune Cell Type Correlation (Spearman ρ)",
               fontsize=8, fontsize_row=8, fontsize_col=8,
               border_color=NA, silent=TRUE)

png(file.path(OUT_DIR, "Figure_Ferroptosis_Immune_CorHeatmap.png"), width=10, height=8,
    units="in", res=300)
grid::grid.draw(p5$gtable)
dev.off()
tiff(file.path(OUT_DIR, "Figure_Ferroptosis_Immune_CorHeatmap.tif"), width=10, height=8,
     units="in", res=300)
grid::grid.draw(p5$gtable)
dev.off()
cat("  Saved ferroptosis-immune correlation heatmap\n")

# ============================================================
# Step 9: Summary and integration
# ============================================================
cat("\n=== Step 9: Summary ===\n")

# Create comprehensive summary
summary_list <- list(
  Ferroptosis_GSVA = gsva_results,
  Immune_CellTypes = immune_results,
  TF_Activity = tf_results,
  External_Signatures = sig_results
)
if(!is.null(module_stats)) {
  summary_list$WGCNA_Modules <- module_stats
}

cat("\n  Key Findings:\n")
cat("  -------------\n")

# Top findings
if(!is.null(module_stats) && nrow(module_stats) > 0) {
  top_modules <- head(module_stats[order(module_stats$FDR), ], 3)
  cat(sprintf("  1. WGCNA identified %d co-expression modules; %d significantly associated with FHL\n",
              nrow(module_stats), sum(module_stats$FDR < 0.05)))
  if(nrow(top_modules) > 0) {
    cat(sprintf("     Top module: %s (r=%.3f, FDR=%.2e)\n",
                top_modules$Module[1], top_modules$Cor_FHL[1], top_modules$FDR[1]))
  }
} else {
  cat("  1. WGCNA module summary unavailable in this run; run wgcna_standalone.R before deeper_analysis_v3.R for module integration.\n")
}

tf_activated <- tf_results[tf_results$Delta > 0, ]
tf_repressed <- tf_results[tf_results$Delta < 0, ]
most_activated <- if(nrow(tf_activated) > 0) {
  tf_activated$TF[which.max(tf_activated$Delta)]
} else {
  "none detected"
}
most_repressed <- if(nrow(tf_repressed) > 0) {
  tf_repressed$TF[which.min(tf_repressed$Delta)]
} else {
  "none detected"
}
cat(sprintf("  2. TF activity inference: most activated: %s; most repressed: %s\n",
            most_activated, most_repressed))

top_cell <- head(immune_results[order(immune_results$P_value), ], 3)
cat(sprintf("  3. Top enriched immune cell types: %s (Δ=%.3f)\n",
            paste(sprintf("%s", top_cell$CellType), collapse=", "),
            top_cell$Delta[1]))

cat(sprintf("  4. Ferroptosis Core pathway: GSVA Δ=%.3f (FDR=%.2e)\n",
            gsva_results$Delta[gsva_results$GeneSet == "Ferroptosis_Core"],
            gsva_results$FDR[gsva_results$GeneSet == "Ferroptosis_Core"]))

cat(sprintf("\n  Output directory: %s\n", OUT_DIR))
cat("  Generated files:\n")
for(f in sort(list.files(OUT_DIR))) {
  if(grepl("Figure|Table", f)) {
    size <- file.info(file.path(OUT_DIR, f))$size
    cat(sprintf("    %s (%s bytes)\n", f, format(size, big.mark=",")))
  }
}

cat("\n=== Deeper Analysis v3 Complete ===\n")
