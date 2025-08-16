library(Seurat)

source("R/validation.R")

plot_umap <- function(seurat_obj) {
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
  if (is_seurat_obj(seurat_obj)) {
    if (marker %in% rownames(seurat_obj)) {
      vals <- FetchData(seurat_obj, vars = marker)[, 1]

      colors <- if (all(vals == 0)) {
        c("lightgrey", "lightgrey")
      } else {
        c("lightgrey", "blue")
      }

      feature_plot <- FeaturePlot(
        seurat_obj,
        features = marker,
        pt.size = 0.1,
        cols = colors
      )

      return(feature_plot)
    }
  } else {
    validate(need(FALSE, "Feature plotting function did not receive a seurat object"))
  }
}