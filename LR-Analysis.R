library(Seurat)
library(BulkSignalR)
library(SingleCellSignalR)
library(ggplot2)

run_scsr <- function(seu_obj, prefix) {
  expr_mat   <- as.matrix(GetAssayData(seu_obj, assay = "RNA", layer = "data"))
  cell_types <- as.character(seu_obj$celltype)
  
  # 1. Build network object
  scsr_obj <- SCSRNoNet(
    expr_mat,
    normalize       = FALSE,
    method          = "log-only",
    min.count       = 1,
    prop            = 0.01,
    log.transformed = TRUE,
    populations     = cell_types
  )
  
  # 2. Infer paracrine interactions with balanced sensitivity
  scsr_obj <- performInferences(
    scsr_obj,
    verbose   = TRUE,
    min.logFC = 0.01,
    max.pval  = 0.1
  )
  
  # 3. Extract all cell-pair tables
  para_list <- paracrines(scsr_obj)
  res_table <- do.call(rbind, lapply(names(para_list), function(pair) {
    df <- para_list[[pair]]
    if (!is.null(df) && nrow(df) > 0) {
      df$cell_pair <- pair
      df$condition <- prefix
    }
    return(df)
  }))
  
  write.csv(res_table, paste0(prefix, "_LR_interactions.csv"), row.names = FALSE)
  return(list(obj = scsr_obj, table = res_table))
}

# Run on both Primary Tumor and Metastasis
pt_res  <- run_scsr(pt_obj, "Primary_Tumor")
met_res <- run_scsr(met_obj, "Metastasis")

# Check interaction counts
cat("Primary interactions:", nrow(pt_res$table), "\n")
cat("Metastasis interactions:", nrow(met_res$table), "\n")

# Combine top interactions and plot comparative Bubble Plot
combined_lr <- rbind(pt_res$table, met_res$table)
combined_lr$LR_pair <- paste(combined_lr$L, combined_lr$R, sep = " -> ")

top_pairs <- unique(combined_lr[order(-combined_lr$LR.score), "LR_pair"])[1:min(20, length(unique(combined_lr$LR_pair)))]
plot_df   <- combined_lr[combined_lr$LR_pair %in% top_pairs, ]

p <- ggplot(plot_df, aes(x = cell_pair, y = LR_pair, size = LR.score, color = condition)) +
  geom_point(alpha = 0.85, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Primary_Tumor" = "#2b5c8f", "Metastasis" = "#d95f02")) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Ligand-Receptor Signaling: Primary HCC vs Metastasis",
    x     = "Cell-Cell Interaction",
    y     = "Ligand -> Receptor",
    size  = "LR Score",
    color = "Site"
  )

print(p)
ggsave("HCC_LR_Comparison_BubblePlot.png", plot = p, width = 10, height = 7, dpi = 300)