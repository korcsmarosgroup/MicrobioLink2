usethis::use_package("readr", type = "Imports")

#'  normalize_all_samples
#'
#'  Normalize all QC-filtered .h5ad files from each sample directory.
#'  For each sample in the output_dir:
#'    - Loads its 'QC_filtered.h5ad' file
#'    - Performs total-count normalization (counts per 10,000)
#'    - Applies log1p transformation
#'    - Saves the normalized file as 'normalized.h5ad' in the same folder
#'
#'  @details
#'  ERROR CODE 13: Normalization failed.
#'
#'  @param output_dir (str) Path to the main output directory containing per-sample subfolders.
#'  @param log_file (str): Path to the central log file.
#'
#'  @export

normalize_all_samples <- function(output_dir, HVG_selection, number_top_genes, log_file) {
  write_log(paste0("Starting normalisation of all samples in ", output_dir), log_file)

  tryCatch(
    {
      dir_list <- list.dirs(output_dir)
      for (sample in dir_list) {
        sample_path <- file.path(output_dir, sample)

        input_rds <- file.path(sample_path, "QC_filtered.rds")
        output_rds <- file.path(sample_path, "normalized.rds")
        output_tsv <- file.path(sample_path, "normalized.tsv")
        output_hvg_rds <- file.path(sample_path, "normalized_hvg.rds")
        output_hvg_tsv <- file.path(sample_path, "normalized_hvg.tsv")

        if (!file.exists(input_rds)) {
          write_log(paste0("Skipping ", sample, " : No QC_filtered.rds found"), log_file)
          next
        }

        write_log(paste0("Normalizing ", sample), log_file)
        sc10x <- readRDS(input_rds)

        if (isTRUE(HGV_selection)) {
          sc10x_hvg <- FindVariableFeatures(sc10x, selection.method = "vst", nfeatures = number_top_genes)
          saveRDS(sc10x_hvg, file = output_hvg_rds)

          active_assay <- DefaultAssay(sc10x_hvg)
          counts <- LayerData(sc10x, assay = active_assay, layer = "counts")
          sc10x_hvg_counts <- counts[sc10x_hvg, drop = FALSE]
          sc10x_hvg_dense <- as.matrix(t(sc10x_hvg_counts))

          write.table(
            sc10x_hvg_dense,
            file = output_hvg_tsv,
            sep = "\t",
            quote = FALSE,
            col.names = NA
          )
          write_log(paste0("Saved normalized HVG data for ", sample, " to ", output_hvg_rds), log_file)
        }
        sc10x_normalized <- NormalizeData(sc10x, normalization.method = "LogNormalize", scale.factor = 10000)
        saveRDS(sc10x_normalized, file = output_rds)
        active_assay <- DefaultAssay(sc10x)
        normalized_counts <- LayerData(sc10x_normalized, assay = active_assay, layer = "data")
        normalized_dense <- as.matrix(t(normalized_counts))
        write.table(
          normalized_dense,
          file = output_tsv,
          sep = "\t",
          quote = FALSE,
          col.names = NA
        )
      }
      write_log("All sample normalizations completed successfully", log_file)
    },
    error = function(e) {
      msg <- paste0("ERROR CODE 13: Normalization failed due to", conditionMessage(e))
      write_log(msg, log_file)
    }
  )
}
