library("BPCells")

is_multimodal <- function(seurat_obj.data) {
  if (is.list(seurat_obj.data)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is_seurat_obj <- function(obj) {
  if (!inherits(obj, "Seurat")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

get_filetype <- function(file_df) {
  filetype <- tolower(tools::file_ext(file_df$name))
  validate(need(filetype %in% c("h5", "h5ad", "rds"), "Invalid file type"))

  return(filetype)
}

read_file <- function(file_df) {
  filepath <- file_df$datapath
  filetype <- get_filetype(file_df)

  tryCatch(
    {
      if (filetype == "h5") {
        data <- Read10X_h5(filepath)
        return(data)
      } else if (filetype == "h5ad") {
        data <- ReadH5AD(filepath)
      }
    },
    error = function(e) {
      validate(need(FALSE, "Failed to read file"))
    }
  )
}

get_features <- function(seurat_obj) {
  if (is_seurat_obj(seurat_obj)) {
    markers <- FindAllMarkers(seurat_obj, min.pct = 0.1, only.pos = TRUE)
    return(markers)
  } else {
    validate(need(FALSE, "get_features needs a Seurat object"))
  }
}

get_cluster_labels <- function(seurat_obj) {
  if (is_seurat_obj(seurat_obj)) {
    cluster_labels <- levels(Idents(seurat_obj))
    return(cluster_labels)
  } else {
    validate(need(FALSE, "get_cluster_labels needs a Seurat object"))
  }
}

filter_markers <- function(markers, cluster, sort_order) {
  if (cluster != "All Clusters") {
    markers <- filter(markers, cluster == cluster)
  }

  if (sort_order == 1) {
    markers <- arrange(markers, desc(avg_log2FC))
  } else if (sort_order == 2) {
    markers <- arrange(markers, desc(pct.1))
  }

  return(rownames(markers))
}