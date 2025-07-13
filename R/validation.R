is_multimodal <- function(seurat_obj.data) {
  if (is.list(seurat_obj.data) && all(sapply(seurat_obj.data, is.matrix))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

get_filetype <- function(file_df) {
  filetype <- tolower(tools::file_ext(file_df$name))
  validate(need(filetype %in% c("h5"), "Upload .h5 file"))

  return(filetype)
}