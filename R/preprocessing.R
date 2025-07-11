library(Seurat)
library(sctransform)
library(DoubletFinder)
library(tidyverse)
library(patchwork)

seurat_obj.data <- Read10X_h5(filename = "data/GSM8352048_MA5_filtered_feature_bc_matrix.h5")
seurat_obj <- CreateSeuratObject(counts = seurat_obj.data$"Gene Expression", project = "seurat_obj", min.cells = 3, min.features = 200)

str(seurat_obj)

seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
upper_lim <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.98)
lower_lim <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.02)

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > lower_lim & nFeature_RNA < upper_lim & percent.mt < 10)

seurat_obj <- SCTransform(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

