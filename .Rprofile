source("renv/activate.R")

if (requireNamespace("reticulate", quietly = TRUE)) {
  if (Sys.getenv("USE_CONDA_FOR_R") == "TRUE") {
    reticulate::use_condaenv("scrna_py", required = TRUE)
  }
}