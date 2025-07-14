plot_umap <- function(seurat_obj) {
  library(Seurat)

  source("R/validation.R")

  if (is_seurat_obj(seurat_obj)) {
    umap <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE)

    return(umap)
  } else {
    validate(need(FALSE, "UMAP plotting function did not receive a seurat object"))
  }
}