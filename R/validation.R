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
  validate(need(filetype %in% c("h5", "rds"), "Invalid file type"))

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
      }
    },
    error = function(e) {
      validate(need(FALSE, "Failed to read file"))
    }
  )
}