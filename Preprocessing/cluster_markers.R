library(Seurat)
library(presto)
library(SingleR)
library(SingleCellExperiment)
library(tidyverse)
library(patchwork)

seurat_obj <- readRDS("data/Seurat Objects/GSE217494_processed.rds")

# Finding cluster markers
seurat_obj.markers <- FindAllMarkers(seurat_obj, min.pct = 0.5, only.pos = TRUE)

head(arrange(filter(seurat_obj.markers, avg_log2FC > 1), desc(avg_log2FC)), 20)

all_top10 <- seurat_obj.markers %>%
  group_by(cluster) %>%
  filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup()

p6 <- DoHeatmap(
  seurat_obj,
  features = all_top10$gene
) +
  scale_fill_gradientn(
    colors = c("magenta", "black", "yellow"),
    breaks = c(-2, 2),
    labels = c("Low", "High")
  )
ggsave(
  filename = "data/plots/markers_heatmap.png",
  plot = p6,
  device = "png",
  width = 15,
  height = 10
)

# Using SingleR as reference
heart_ref <- celldex::HumanPrimaryCellAtlasData()

cluster_results <- SingleR(
  test = as.SingleCellExperiment(seurat_obj),
  ref = heart_ref,
  labels = heart_ref$label.main
)

seurat_obj$singler_labels <- cluster_results$labels

p7 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "singler_labels",
  pt.size = 0.1,
  label = FALSE,
  repel = TRUE
)
ggsave(
  filename = "data/plots/singleR_umap.png",
  plot = p7,
  device = "png",
  width = 20,
  height = 15
)

# Labeling clusters
cluster_labels <- c(
  "0" = "Endothelial Cells",
  "1" = "Endothelial Cells",
  "2" = "Macrophages",
  "3" = "Fibroblasts",
  "4" = "Pericytes",
  "5" = "Vascular SMC",
  "6" = "Fibroblasts",
  "7" = "T/NK Cells",
  "8" = "Dendritic Cells",
  "9" = "Monocytes",
  "10" = "Glial Cells",
  "11" = "B Cells",
  "12" = "Cardiomyocytes"
)

seurat_obj <- RenameIdents(seurat_obj, cluster_labels)

table(Idents(seurat_obj))

# saveRDS(seurat_obj, file = "data/Seurat Objects/GSE217494_labeled.rds")

p8 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  pt.size = 0.1,
  label = FALSE,
  repel = TRUE
)
ggsave(
  filename = "data/plots/labeled_umap.png",
  plot = p8,
  device = "png",
  width = 20,
  height = 15
)
