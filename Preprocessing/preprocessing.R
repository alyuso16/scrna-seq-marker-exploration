library(BPCells)
library(Seurat)
library(sctransform)
library(SingleCellExperiment)
library(scDblFinder)
library(tidyverse)
library(patchwork)

# Reading and merging samples
data_dir <- "data/GSE217494_RAW"
sample_names <- c("sample13", "sample28", "sample32", "sample33", "sample34", "sample39", "sample41", "sample4", "sample5", "sample8", "sample42")

combined_samples.data <- sapply(sample_names, function(i) {
  sample.data <- Read10X(file.path(data_dir, i))$"Gene Expression"
  colnames(sample.data) <- paste(sapply(strsplit(colnames(sample.data), split = "-"), '[[', 1L), i, sep = "-")
  sample.data
})

seurat_obj.data <- do.call("cbind", combined_samples.data)
counts <- seurat_obj.data

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "seurat_obj",
  min.cells = 20
)

saveRDS(seurat_obj, file = "data/Seurat Objects/GSE217494_raw.rds")

# Using BPCells
write_matrix_dir(
  mat = seurat_obj@assays$RNA@layers$counts,
  dir = "data/bpcells/GSE217494_bp_counts"
)

counts.mat <- open_matrix_dir(dir = "data/bpcells/GSE217494_bp_counts")
seurat_obj@assays$RNA@layers$counts <- counts.mat

saveRDS(seurat_obj, file = "data/Seurat Objects/GSE217494_raw_bp.rds")

# seurat_obj <- readRDS("data/Seurat Objects/GSE217494_raw_bp.rds")

# QC
seurat_obj <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-",
  col.name = "percent.mt"
)

p1 <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts"
)
ggsave(filename = "data/plots/raw_metrics_vln.png", plot = p1, device = "png")

upper_lim_nfeature <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.98)
lower_lim_nfeature <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.02)
upper_lim_ncount <- quantile(seurat_obj@meta.data$nCount_RNA, probs = 0.98)

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > lower_lim_nfeature & nFeature_RNA < upper_lim_nfeature & nCount_RNA < upper_lim_ncount & percent.mt < 15
)

p2 <- VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  layer = "counts"
)
ggsave(filename = "data/plots/post_QC_metrics_vln.png", plot = p2, device = "png")

# Normalization
seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2")

p3 <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)
ggsave(filename = "data/plots/counts_features_scatter.png", plot = p3, device = "png")

# Doublet detection
sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")
assay(sce, "counts") <- as(GetAssayData(seurat_obj, assay = "RNA", slot = "counts"), "dgCMatrix")
sce <- scDblFinder(sce)

seurat_obj$scDblFinder_score <- sce$scDblFinder.score
seurat_obj$scDblFinder_class <- sce$scDblFinder.class

seurat_obj <- subset(seurat_obj, subset = scDblFinder_class == "singlet")

# Choosing PCs
seurat_obj <- RunPCA(seurat_obj)

p4 <- ElbowPlot(
  seurat_obj,
  ndims = 50,
  reduction = "pca"
)
ggsave(filename = "data/plots/elbow_plot.png", plot = p4, device = "png")

cutoff_pc <- 11

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:cutoff_pc)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.3)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:cutoff_pc)

# saveRDS(seurat_obj, file = "data/Seurat Objects/GSE217494_processed.rds")
# seurat_obj <- readRDS("data/Seurat Objects/GSE217494_processed.rds")

p5 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  pt.size = 0.1,
  label = TRUE,
  repel = TRUE
)
ggsave(
  filename = "data/plots/unlabeled_umap.png",
  plot = p5,
  device = "png",
  width = 20,
  height = 15
)
