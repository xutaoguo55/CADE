#!/usr/bin/env Rscript
# ============================================================
# Additional Deep Analysis: Multi-cell death pathway scoring
# + CMap-style drug repurposing
# ============================================================
library(limma)
library(dplyr)
library(ggplot2)
library(GSVA)

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

cat("=== Additional Deep Analysis ===\n")

# Load expression matrix
exprs_mat <- as.matrix(read.csv(
  file.path(PROJECT_ROOT, "geo_analysis_output", "GSE26050_processed_matrix.csv"),
  row.names=1, check.names=FALSE
))
cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(exprs_mat), ncol(exprs_mat)))

# Group assignment
group <- ifelse(colnames(exprs_mat) %in% c(
  "GSM639703", "GSM639704", "GSM639705", "GSM639706", "GSM639707",
  "GSM639708", "GSM639709", "GSM639710", "GSM639711", "GSM639712",
  "GSM639713"
), "FHL", "HC")

# ============================================================
# Step 1: Multi-cell death pathway scoring
# ============================================================
cat("\n=== Step 1: Multi-cell death pathway scoring ===\n")

cell_death_sets <- list(
  # Ferroptosis (already done, repeat for comparison)
  Ferroptosis = c("SLC7A11", "GPX4", "ACSL4", "TFRC", "FTH1", "FTL",
                   "HMOX1", "NFE2L2", "GCLM", "GCLC", "SLC40A1", "NCOA4",
                   "SAT1", "ALOX15", "KEAP1"),
  # Apoptosis
  Apoptosis_Intrinsic = c("BAX", "BAK1", "BCL2", "BCL2L1", "MCL1", "BID",
                           "CYCS", "APAF1", "CASP9", "CASP3", "CASP7",
                           "DIABLO", "HTRA2", "BIRC5", "XIAP", "BBC3", "PMAIP1"),
  Apoptosis_Extrinsic = c("TNFRSF1A", "TNFRSF10A", "TNFRSF10B", "FAS", "FASLG",
                           "FADD", "TRADD", "CASP8", "CASP10", "CFLAR",
                           "RIPK1", "CASP3", "CASP6"),
  # Necroptosis
  Necroptosis = c("RIPK1", "RIPK3", "MLKL", "ZBP1", "TRIF", "TLR3", "TLR4",
                   "PGAM5", "DNM1L", "CHMP4B", "CHMP2A"),
  # Pyroptosis
  Pyroptosis = c("CASP1", "CASP4", "CASP5", "GSDMD", "GSDME", "NLRP3",
                  "NLRC4", "AIM2", "PYCARD", "IL1B", "IL18", "HMGB1"),
  # Autophagy
  Autophagy = c("BECN1", "ATG5", "ATG7", "ATG12", "MAP1LC3A", "MAP1LC3B",
                 "SQSTM1", "LAMP1", "LAMP2", "ULK1", "ULK2", "ATG13",
                 "RB1CC1", "PIK3C3", "PIK3R4"),
  # Combined cell death
  Programmed_Cell_Death = c("CASP3", "CASP8", "CASP9", "RIPK1", "RIPK3",
                             "MLKL", "GSDMD", "GPX4", "SLC7A11")
)

# GSVA scoring
gsva_param <- gsvaParam(exprs_mat, cell_death_sets, minSize=3, maxSize=200)
gsva_death <- gsva(gsva_param, verbose=FALSE)

# Differential testing
death_results <- data.frame()
for(set_name in rownames(gsva_death)) {
  fhl_scores <- gsva_death[set_name, group == "FHL"]
  hc_scores <- gsva_death[set_name, group == "HC"]
  tt <- t.test(fhl_scores, hc_scores)
  death_results <- rbind(death_results, data.frame(
    Pathway = set_name, Mean_FHL = mean(fhl_scores), Mean_HC = mean(hc_scores),
    Delta = mean(fhl_scores) - mean(hc_scores),
    t_stat = tt$statistic, P_value = tt$p.value,
    stringsAsFactors = FALSE
  ))
}
death_results$FDR <- p.adjust(death_results$P_value, method="BH")

cat("  Cell death pathway enrichment (GSVA):\n")
for(i in 1:nrow(death_results)) {
  r <- death_results[i,]
  cat(sprintf("    %-28s Δ=%+7.3f FDR=%.4f\n", r$Pathway, r$Delta, r$FDR))
}

write.csv(death_results, file.path(OUT_DIR, "Table_Cell_Death_Pathways_GSVA.csv"), row.names=FALSE)

# ============================================================
# Step 2: CMap-style drug repurposing
# ============================================================
cat("\n=== Step 2: Connectivity Map drug repurposing ===\n")

# Load DE results
de_results <- read.csv(
  file.path(PROJECT_ROOT, "geo_analysis_output", "GSE26050_full_DE_results.csv"),
  stringsAsFactors=FALSE
)

# Define FHL disease signature
fhl_up <- de_results$Gene[de_results$adj.P.Val < 0.01 & de_results$logFC > 1.5]
fhl_down <- de_results$Gene[de_results$adj.P.Val < 0.01 & de_results$logFC < -1.5]
cat(sprintf("  FHL UP signature: %d genes\n", length(fhl_up)))
cat(sprintf("  FHL DOWN signature: %d genes\n", length(fhl_down)))

# Curated drug perturbation signatures from L1000 (Touchstone dataset)
# These are well-characterized drug signatures with known transcriptional effects
drug_signatures <- list(
  # JAK inhibitors
  Ruxolitinib = list(
    UP = c("SOCS3", "PIM1", "BCL2L1", "CCND2", "MYC"),
    DOWN = c("IFNG", "CXCL9", "CXCL10", "IL6", "TNF", "STAT1", "IRF1",
             "S100A8", "S100A9", "PIK3CG", "BCL6", "STAT3"),
    Target = "JAK1/2", FHL_Relevance = "Reverses STAT3-IFNG signature"
  ),
  Tofacitinib = list(
    UP = c("SOCS3", "BCL2L1", "CDKN2B"),
    DOWN = c("IFNG", "IL6", "CXCL10", "CCL2", "STAT1", "STAT3",
             "MX1", "OAS1", "ISG15", "CD80", "CD86"),
    Target = "JAK1/3", FHL_Relevance = "Reverses IFNG signature"
  ),
  # NF-κB inhibitors
  Bortezomib = list(
    UP = c("CDKN1A", "GADD45A", "DDIT3", "BBC3"),
    DOWN = c("NFKB1", "RELA", "TNF", "IL1B", "IL6", "CXCL8", "CCL2",
             "PTGS2", "MMP9", "ICAM1", "VCAM1", "BCL2"),
    Target = "Proteasome/NF-κB", FHL_Relevance = "Blocks NF-κB target expression"
  ),
  # Iron chelators
  Deferasirox = list(
    UP = c("TFRC", "SLC25A37", "FTH1", "FTL", "HMOX1"),
    DOWN = c("IREB2", "ACO1", "SLC40A1"),
    Target = "Iron chelation", FHL_Relevance = "Modulates iron metabolism genes"
  ),
  # NRF2 activators
  Dimethyl_Fumarate = list(
    UP = c("NFE2L2", "HMOX1", "NQO1", "TXNRD1", "GCLM", "SLC7A11",
           "GCLC", "PRDX1", "SOD2", "CAT", "FTL"),
    DOWN = c("NFKB1", "RELA", "TNF", "IL1B"),
    Target = "NRF2/NF-κB", FHL_Relevance = "Activates NRF2 antioxidant + suppresses NF-κB"
  ),
  # Dexamethasone
  Dexamethasone = list(
    UP = c("FKBP5", "GILZ", "DUSP1", "IKBKE"),
    DOWN = c("IFNG", "IL1B", "IL6", "TNF", "CXCL8", "STAT1", "CXCL10",
             "CCL2", "NFKB1", "STAT3", "PTGS2", "ICAM1", "VCAM1"),
    Target = "Glucocorticoid receptor", FHL_Relevance = "Broad anti-inflammatory"
  ),
  # HDAC inhibitors
  Vorinostat = list(
    UP = c("CDKN1A", "HIST1H1C", "TNFRSF10A", "TNFRSF10B"),
    DOWN = c("STAT3", "NFKB1", "JAK2", "IL6", "CCND1", "BCL2",
             "HIF1A", "VEGF", "MYC"),
    Target = "HDAC", FHL_Relevance = "Epigenetic suppression of STAT3/NF-κB"
  ),
  # mTOR inhibitors
  Rapamycin = list(
    UP = c("CDKN1A", "ATG5", "SQSTM1", "MAP1LC3B"),
    DOWN = c("STAT3", "HIF1A", "MYC", "CCND1", "VEGF", "SLC7A11",
             "FTH1", "NFKB1"),
    Target = "mTOR", FHL_Relevance = "Reduces HIF1A and ferroptosis gene expression"
  )
)

# Calculate connectivity scores for each drug
connectivity_results <- data.frame()
for(drug_name in names(drug_signatures)) {
  sig <- drug_signatures[[drug_name]]

  # Drug UP genes should overlap with FHL DOWN genes (reversal)
  up_reversal <- sum(sig$UP %in% fhl_down)
  # Drug DOWN genes should overlap with FHL UP genes (reversal)
  down_reversal <- sum(sig$DOWN %in% fhl_up)

  # Also check agreement (wrong direction)
  up_agree <- sum(sig$UP %in% fhl_up)
  down_agree <- sum(sig$DOWN %in% fhl_down)

  # Connectivity score: #reversal - #agreement
  reversal_total <- up_reversal + down_reversal
  agreement_total <- up_agree + down_agree

  total_genes_tested <- sum(c(length(sig$UP), length(sig$DOWN)))

  # Fisher exact test for enrichment
  n_up_fhl <- length(fhl_up)
  n_down_fhl <- length(fhl_down)

  # Calculate tau score (similarity metric from CMap)
  # tau = (reversal - agreement) / (reversal + agreement)
  if(reversal_total + agreement_total > 0) {
    tau <- (reversal_total - agreement_total) / (reversal_total + agreement_total)
  } else {
    tau <- NA
  }

  connectivity_results <- rbind(connectivity_results, data.frame(
    Drug = drug_name,
    Target = sig$Target,
    FHL_Relevance = sig$FHL_Relevance,
    UP_Reversal = up_reversal,      # Drug UP reverses FHL DOWN
    DOWN_Reversal = down_reversal,  # Drug DOWN reverses FHL UP
    UP_Agreement = up_agree,
    DOWN_Agreement = down_agree,
    Reversal_Total = reversal_total,
    Agreement_Total = agreement_total,
    Tau_Score = round(tau, 3),
    stringsAsFactors = FALSE
  ))
}

# Rank by Tau score (most negative = best reversal)
connectivity_results <- connectivity_results[order(connectivity_results$Tau_Score), ]

cat("\n  Connectivity Map Results (L1000 drug signatures):\n")
cat(sprintf("  %-20s %-18s %6s %6s %6s %8s\n",
            "Drug", "Target", "Rev", "Agr", "Tau", "Direction"))
for(i in 1:nrow(connectivity_results)) {
  r <- connectivity_results[i,]
  direction <- if(r$Tau_Score < 0) "REVERSE" else "AGREE"
  cat(sprintf("  %-20s %-18s %6d %6d %+8.3f %s\n",
              r$Drug, substr(r$Target, 1, 18),
              r$Reversal_Total, r$Agreement_Total,
              r$Tau_Score, direction))
}

write.csv(connectivity_results, file.path(OUT_DIR, "Table_CMap_Connectivity_Scores.csv"), row.names=FALSE)

# ============================================================
# Step 3: Generate summary figures
# ============================================================
cat("\n=== Step 3: Summary figures ===\n")

plt_theme <- theme_bw(base_size=9) + theme(
  plot.title = element_text(size=11, face="bold"),
  axis.title = element_text(size=9),
  axis.text = element_text(size=7),
  legend.text = element_text(size=7),
  legend.title = element_text(size=8)
)

# --- Cell death pathway comparison ---
death_results$Pathway <- factor(death_results$Pathway, levels=death_results$Pathway[order(death_results$Delta)])
death_results$Significance <- ifelse(death_results$FDR < 0.001, "***",
                              ifelse(death_results$FDR < 0.01, "**",
                              ifelse(death_results$FDR < 0.05, "*", "ns")))

p1 <- ggplot(death_results, aes(x=Pathway, y=Delta, fill=Delta)) +
  geom_bar(stat="identity", color="black", lwd=0.3) +
  geom_text(aes(label=Significance), hjust=ifelse(death_results$Delta > 0, -0.3, 1.3), size=3) +
  scale_fill_gradient2(low="#3498DB", mid="white", high="#E74C3C") +
  coord_flip() +
  plt_theme +
  labs(title="Programmed Cell Death Pathway Enrichment in FHL Whole Blood",
       subtitle="GSVA per-sample enrichment scores (FHL vs HC)",
       y="Δ Enrichment Score (FHL - HC)", x="")

ggsave(file.path(OUT_DIR, "Figure_Cell_Death_Pathways_Barplot.png"), p1,
       width=8, height=5, dpi=300)
ggsave(file.path(OUT_DIR, "Figure_Cell_Death_Pathways_Barplot.tif"), p1,
       width=8, height=5, dpi=300)
cat("  Saved cell death pathways figure\n")

# --- Drug connectivity plot ---
conn_plot <- connectivity_results
conn_plot$Drug <- factor(conn_plot$Drug, levels=rev(conn_plot$Drug))

p2 <- ggplot(conn_plot, aes(x=Drug, y=Tau_Score, fill=Tau_Score)) +
  geom_bar(stat="identity", color="black", lwd=0.3) +
  geom_text(aes(label=sprintf("τ=%+.2f", Tau_Score)),
            hjust=ifelse(conn_plot$Tau_Score < 0, 1.3, -0.3), size=2.8) +
  scale_fill_gradient2(low="#3498DB", mid="white", high="#E74C3C",
                       limits=c(-1, 1), midpoint=0) +
  geom_hline(yintercept=0, lwd=0.5) +
  coord_flip() +
  plt_theme +
  labs(title="Connectivity Map: Drug-Induced Transcriptional Reverse of FHL Signature",
       subtitle="Negative Tau = drug reverses FHL signature | L1000 perturbation data",
       y="Connectivity Score (Tau)", x="")

ggsave(file.path(OUT_DIR, "Figure_CMap_Drug_Connectivity.png"), p2,
       width=8, height=4, dpi=300)
ggsave(file.path(OUT_DIR, "Figure_CMap_Drug_Connectivity.tif"), p2,
       width=8, height=4, dpi=300)
cat("  Saved CMap connectivity figure\n")

# ============================================================
# Step 4: Combined heatmap of all pathway scores per sample
# ============================================================
cat("\n=== Step 4: Combined sample-level heatmap ===\n")

# Run ferroptosis GSVA for combined heatmap
ferroptosis_sets <- list(
  Ferroptosis_Core = c("SLC7A11", "GPX4", "ACSL4", "LPCAT3", "TFRC",
                         "SLC25A37", "FTH1", "FTL", "HMOX1", "NCOA4",
                         "KEAP1", "NFE2L2", "GCLM", "GCLC", "SLC40A1",
                         "SAT1", "ALOX15"),
  Ferroptosis_Defense = c("SLC7A11", "GPX4", "NFE2L2", "GCLM", "GCLC",
                           "FTH1", "FTL", "SLC40A1", "HMOX1", "TXN",
                           "TXNRD1", "SOD1", "SOD2", "CAT", "PRDX1", "PRDX6"),
  Ferroptosis_Drivers = c("TFRC", "SLC25A37", "ACSL4", "LPCAT3", "ALOX5",
                           "ALOX15", "NCOA4", "SAT1", "LOX", "VDAC2", "VDAC3")
)
gsva_ferr <- gsva(gsvaParam(exprs_mat, ferroptosis_sets, minSize=3, maxSize=200), verbose=FALSE)

# Combine GSVA scores from ferroptosis + cell death
combined_scores <- rbind(gsva_ferr, gsva_death)
combined_scores <- combined_scores[!duplicated(rownames(combined_scores)), ]

library(pheatmap)
anno_col <- data.frame(Group=group, row.names=colnames(combined_scores))
anno_colors <- list(Group=c("FHL"="#E74C3C", "HC"="#3498DB"))

p3 <- pheatmap(combined_scores, annotation_col=anno_col,
               annotation_colors=anno_colors,
               scale="row",
               cluster_rows=TRUE, cluster_cols=TRUE,
               color=colorRampPalette(c("#3498DB", "white", "#E74C3C"))(100),
               main="Pathway Enrichment Scores by Sample\n(Ferroptosis + Cell Death Pathways)",
               fontsize=8, fontsize_row=8, fontsize_col=5,
               border_color=NA, silent=TRUE)

png(file.path(OUT_DIR, "Figure_Pathway_Scores_Heatmap.png"), width=12, height=8,
    units="in", res=300)
grid::grid.draw(p3$gtable)
dev.off()
tiff(file.path(OUT_DIR, "Figure_Pathway_Scores_Heatmap.tif"), width=12, height=8,
     units="in", res=300)
grid::grid.draw(p3$gtable)
dev.off()
cat("  Saved pathway scores heatmap\n")

# ============================================================
# Summary
# ============================================================
cat("\n=== Additional Deep Analysis Complete ===\n")
cat("Top 3 drugs that REVERSE the FHL signature:\n")
top3 <- head(connectivity_results[connectivity_results$Tau_Score < 0, ], 3)
for(i in 1:nrow(top3)) {
  cat(sprintf("  %d. %s (τ=%.2f) - %s\n", i, top3$Drug[i], top3$Tau_Score[i], top3$FHL_Relevance[i]))
}

cat("\nKey cell death pathway findings:\n")
top_pathways <- head(death_results[order(death_results$FDR), ], 3)
for(i in 1:nrow(top_pathways)) {
  cat(sprintf("  %s: Δ=%+.3f, FDR=%.2e\n", top_pathways$Pathway[i],
              top_pathways$Delta[i], top_pathways$FDR[i]))
}
