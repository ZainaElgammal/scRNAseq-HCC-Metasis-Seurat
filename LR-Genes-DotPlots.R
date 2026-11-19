library(Seurat)
library(ggplot2)

# Top Ligand & Receptor genes identified in our interaction analysis
lr_genes <- c("PKM", "HGF", "FN1", "FGL1", "CD44", "EGFR", "CALM1", "FAS", "B2M", "CD48", "CD2")

# Create a combined label (CellType + Site) for direct Primary vs Metastasis comparison
hcc_seurat$celltype_site <- paste(hcc_seurat$celltype, hcc_seurat$site, sep = " | ")
Idents(hcc_seurat) <- "celltype_site"

# Generate comparative DotPlot
p_dot <- DotPlot(
  hcc_seurat,
  features  = lr_genes,
  cols      = c("#f0f0f0", "#a63603"),
  dot.scale = 7
) +
  RotatedAxis() +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", size = 13)
  ) +
  labs(
    title = "Expression of Key Ligand-Receptor Genes: Primary vs Metastasis",
    x     = "Ligand / Receptor Genes",
    y     = "Cell Type | Site"
  )

print(p_dot)
ggsave("HCC_LR_Genes_DotPlot.png", plot = p_dot, width = 10, height = 8, dpi = 300)