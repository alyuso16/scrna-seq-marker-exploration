is_multimodal <- function(seurat_obj.data) {
  if (is.list(seurat_obj.data) && all(sapply(seurat_obj.data, is.matrix))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}