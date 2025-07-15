plot_umap <- function(seurat_obj) {
  library(Seurat)

  source("R/validation.R")

  if (is_seurat_obj(seurat_obj)) {
    umap <- DimPlot(
      seurat_obj,
      reduction = "umap",
      pt.size = 0.1,
      label = TRUE,
      repel = TRUE
    )

    return(umap)
  } else {
    validate(need(FALSE, "UMAP plotting function did not receive a seurat object"))
  }
}

plot_feature <- function(seurat_obj, marker) {
  library(Seurat)

  source("R/validation.R")

  if (is_seurat_obj(seurat_obj)) {
    if (marker %in% rownames(seurat_obj)) {
      feature_plot <- FeaturePlot(
        seurat_obj,
        features = marker,
        pt.size = 0.1,
        cols = c("lightgrey", "blue")
      )

      return(feature_plot)
    }
  } else {
    validate(need(FALSE, "Feature plotting function did not receive a seurat object"))
  }
}