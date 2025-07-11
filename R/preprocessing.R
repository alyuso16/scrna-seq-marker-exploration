library(Seurat)
library(sctransform)
library(tidyverse)
library(patchwork)

seurat_obj.data <- Read10X_h5(filename = "data/GSM8352048_MA5_filtered_feature_bc_matrix.h5")

seurat_obj <- CreateSeuratObject(counts = seurat_obj.data$"Gene Expression", project = "seurat_obj", min.cells = 3, min.features = 200)

