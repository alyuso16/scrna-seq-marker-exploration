run_preprocessing <- function(file_df) {
  library(Seurat)
  library(sctransform)
  library(DoubletFinder)
  library(tidyverse)

  source("R/validation.R")

  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = NULL, value = 0)

  if (!is.null(progress)) progress$inc(0, message = "Loading data...")

  seurat_obj.data <- read_file(file_df)

  if (is_multimodal(seurat_obj.data)) {
    counts <- seurat_obj.data$"Gene Expression"
  } else {
    counts <- seurat_obj.data
  }

  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = "seurat_obj",
    min.cells = 3,
    min.features = 200
  )

  if (!is.null(progress)) progress$inc(0.05, message = "Filtering data...")

  seurat_obj <- PercentageFeatureSet(
    seurat_obj,
    pattern = "^MT-",
    col.name = "percent.mt"
  )
  upper_lim_nfeature <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.98)
  lower_lim_nfeature <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.02)
  upper_lim_ncount <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.98)

  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > lower_lim_nfeature & nFeature_RNA < upper_lim_nfeature & nCount_RNA < upper_lim_ncount & percent.mt < 15)

  if (!is.null(progress)) progress$inc(0.05, message = "Running SCTransform...")

  seurat_obj <- SCTransform(seurat_obj)

  sct_data <- GetAssayData(seurat_obj, assay = "SCT", slot = "data")
  nonzero_features <- rowSums(sct_data != 0) > 0
  seurat_obj <- subset(seurat_obj, features = rownames(sct_data)[nonzero_features])

  if (!is.null(progress)) progress$inc(0.1, message = "Running PCA...")

  seurat_obj <- RunPCA(seurat_obj)

  if (!is.null(progress)) progress$inc(0.1, message = "Detecting doublets...")

  # Doublet Detection
  sweep.res <- paramSweep(seurat_obj, PCs = 1:10, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  best.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  nExp_poi <- round(0.075 * ncol(seurat_obj))

  seurat_obj <- doubletFinder(
    seurat_obj,
    PCs = 1:10,
    pN = 0.25,
    pK = best.pK,
    nExp = nExp_poi,
    reuse.pANN = NULL,
    sct = TRUE
  )

  if (!is.null(progress)) progress$inc(0.3, message = "Filtering doublets...")

  doublet_col <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
  seurat_obj <- subset(seurat_obj, subset = !!as.name(doublet_col) == "Singlet")

  if (!is.null(progress)) progress$inc(0.1, message = "Selecting PCs...")

  # Choose PCs
  all_pc_std <- seurat_obj[["pca"]]@stdev
  pc_pct_var <- ((all_pc_std ^ 2) / sum(all_pc_std ^ 2)) * 100
  cum_pct_var <- cumsum(pc_pct_var)

  included_pcs <- which(cum_pct_var <= 90 & pc_pct_var >= 5)
  cutoff_pc <- max(included_pcs) + 1

  if (cutoff_pc < 10) {
    cutoff_pc <- 10
  }

  if (!is.null(progress)) progress$inc(0.1, message = "Finding neighbors, clustering, and performing UMAP...")

  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:cutoff_pc)
  seurat_obj <- FindClusters(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:cutoff_pc)

  if (!is.null(progress)) progress$set(message = "Finished", value = 1)

  return(seurat_obj)
}

get_features <- function(seurat_obj) {
  source("R/validation.R")

  if (is_seurat_obj(seurat_obj)) {
    markers <- rownames(seurat_obj)
    return(markers)
  } else {
    validate(need(FALSE, "get_features needs a Seurat object"))
  }
}