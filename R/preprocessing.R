library(Seurat)
library(sctransform)
library(DoubletFinder)
library(tidyverse)
library(patchwork)

seurat_obj.data <- Read10X_h5(filename = "data/GSM8352048_MA5_filtered_feature_bc_matrix.h5")
seurat_obj <- CreateSeuratObject(counts = seurat_obj.data$"Gene Expression", project = "seurat_obj", min.cells = 3, min.features = 200)

str(seurat_obj)

seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
upper_lim_nfeature <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.98)
lower_lim_nfeature <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.02)
upper_lim_ncount <- quantile(seurat_obj@meta.data$nFeature_RNA, probs = 0.98)

seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > lower_lim_nfeature & nFeature_RNA < upper_lim_nfeature & nCount_RNA < upper_lim_ncount & percent.mt < 10)

seurat_obj <- SCTransform(seurat_obj)
saveRDS(seurat_obj, file = "data/GSM8352048_processed.rds")
seurat_obj <- RunPCA(seurat_obj)

# Doublet Detection
sweep.res <- paramSweep(seurat_obj, PCs = 1:10, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
best.pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

nExp_poi <- round(0.075 * ncol(seurat_obj))

seurat_obj <- doubletFinder(seurat_obj, PCs = 1:10, pN = 0.25, pK = best.pK, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)

doublet_col <- grep("DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)
seurat_obj <- subset(seurat_obj, subset = !!as.name(doublet_col) == "Singlet")

all_pc_std <- Stdev(seurat_obj, reduction = "pca")
pct_change_std <- diff(all_pc_std) / all_pc_std[-length(all_pc_std)]
elbow_point <- which.min(pct_change_std) + 1
if(elbow_point < 10) {
    elbow_point <- 10
}

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:elbow_point)
seurat_obj <- FindClusters(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:elbow_point)

DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE)
