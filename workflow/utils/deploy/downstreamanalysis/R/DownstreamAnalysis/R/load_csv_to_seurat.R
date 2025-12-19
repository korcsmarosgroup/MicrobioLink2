usethis::use_package("Seurat", type = "Imports")
#' load_cvs_to_seurat

load_cvs_to_seurat <- function(path, sample, log_file_2) {
  write_log_2(paste0("Loading CSV file ", path), log_file_2)
  df <- read.table(path, sep = "\t")

  # Create Seurat object
  data <- CreateSeuratObject(df, min.cells = 3, min.genes = 200)

  write_log_2(paste0("Successfully loaded CSV for sample ", sample), log_file_2)
}
