# Install packages automatically if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler", update = FALSE, ask = FALSE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# All unique Ligand & Receptor genes from our top metastatic interactions
lr_all_genes <- unique(c(
  "PKM", "HGF", "FN1", "FGL1", "CD44", "EGFR",
  "CALM1", "FAS", "B2M", "KIR3DL2", "KIR2DL4",
  "ANXA1", "CD48", "CD2"
))

# Run GO Biological Process (BP) Enrichment
ego_bp <- enrichGO(
  gene          = lr_all_genes,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.1
)

# Save enrichment table
write.csv(as.data.frame(ego_bp), "Metastatic_HCC_GO_Enrichment.csv", row.names = FALSE)

# Plot top 12 enriched biological pathways
p_go <- dotplot(ego_bp, showCategory = 12) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "bold", size = 10),
    plot.title  = element_text(face = "bold", size = 12)
  ) +
  labs(
    title = "Enriched Biological Processes in Metastatic HCC Crosstalk",
    x     = "Gene Ratio"
  )

print(p_go)
ggsave("Metastatic_HCC_GO_Enrichment.png", plot = p_go, width = 9.5, height = 6.5, dpi = 300)