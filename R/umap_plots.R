plot_umap <- function(seurat_obj) {
  library(Seurat)

  umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE)

  return(umap)
}