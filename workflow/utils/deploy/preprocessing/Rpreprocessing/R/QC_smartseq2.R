usethis::use_package("Seurat", type = "Imports")
usethis::use_package("ggplot2", type = "Imports")
#' QC_smartseq2
#'
#' Perform QC for Smart-seq2 or Smart-seq3 data:
#'
#' @export

QC_smartseq2 <- function(output_folder, log_file, platform = "smartseq2") {

  write_log(paste0("Running QC for ", platform, " samples..."), log_file)

  dir_list <- list.dirs(output_folder, recursive = FALSE, full.names = FALSE)
  for (sample in dir_list) {
    sample_path <- file.path(output_folder, sample)
    count_file <- file.path(sample_path, "counts.tsv")
    matrix_path <- file.path(sample_path, "matrix.mtx")

    # --- Detect matrix type ---
    if (file.exists(matrix_path)) {
      smartseq2_mtx <- Seurat::ReadMtx(
        mtx = matrix_path, features = "features.tsv",
        cells = "barcodes.tsv"
      )
      smartseq2_data <- Seurat::CreateSeuratObject(counts = expression_mtx)
    } else if (file.exists(count_file)) {
      counts <- read.table(count_file)
      counts <- t(counts)
      smartseq2_data <- Seurat::CreateSeuratObject(counts)
    } else {
      write_log(paste0("No count matrix found for ", sample, ". Skipping."), log_file)
      next
    }
  }

  # --- Calculate QC metrics ---
  smartseq2_data[["percent.mt"]] <- Seurat::PercentageFeatureSet(smartseq2_data, pattern = "MT")

  qc_dir <- file.path(sample_path, "QC")
  dir.create(qc_dir)

  # --- Violin plots ---
  smartseq2_vlnplot <- VlnPlot(smartseq2_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  smartseq2_VlnPlot_name <- paste0(sample, "_smartseq2_VlnPlot.png")
  ggplot2::ggsave(filename = smartseq2_VlnPlot_name, plot = smartseq2_vlnplot, path = qc_dir, width = 12, height = 8, dpi = 300)

  # --- Bar plot: total counts ---



  # --- Filtering ---
  smartseq2_filtered <- subset(smartseq2_data, subset = nFeature_RNA < 10000 & percent.mt < 15)

  write_log(paste0("Filtered ", sample, ": retained ", nrow(smartseq2_filtered), " cells."), log_file)

  # --- Save results ---
  saveRDS(smartseq2_filtered, file = paste0(sample, "_filtered_qc.rds"))
  write_log(paste0("QC results saved for ", sample), log_file)
}
