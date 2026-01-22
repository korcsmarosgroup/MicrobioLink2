#' normalize_all_samples
#'
#' Normalize all QC-filtered .h5ad files from each sample directory.
#' For each sample in the output_dir:
#'    - Loads its 'QC_filtered.rds' file
#'    - Performs total-count normalization (counts per 10,000)
#'    - Applies log1p transformation
#'    - Saves the normalized file as 'normalized.rds' in the same folder
#'
#' @details
#' ERROR CODE 13: Normalization failed.
#'
#' @param output_dir (str) Path to the main output directory containing per-sample subfolders.
#' @param log_file (str): Path to the central log file.
#'
#' @export

normalize_all_samples <- function(output_dir, HVG_selection, number_top_genes, log_file) {

  output_dir <- file.path(output_dir, "STARsolo_out")
  write_log(paste0("Starting normalisation of all samples in ", output_dir), log_file)

  tryCatch(
    {
      # Loop through all sample folders
      dir_list <- list.dirs(output_dir)
      for (sample in dir_list) {
        sample_path <- file.path(output_dir, sample)
        if (!dir.exists(sample_path)) next

        input_rds <- file.path(sample_path, "QC_filtered.rds")
        output_rds <- file.path(sample_path, "normalized.rds")
        output_tsv <- file.path(sample_path, "normalized.tsv")
        output_hvg_rds <- file.path(sample_path, "normalized_hvg.rds")
        output_hvg_tsv <- file.path(sample_path, "normalized_hvg.tsv")

        if (!file.exists(input_rds)) {
          write_log(paste0("Skipping ", sample, ": No QC_filtered.rds found"), log_file)
          next
        }

        write_log(paste0("Normalizing sample: ", sample), log_file)
        sc <- readRDS(input_rds)

        # Normalization steps

        sc_normalized <- NormalizeData(sc)
        all_genes <- rownames(sc_normalized)

        # HVG selection

        if (isTRUE(HVG_selection)) {
          sc_normalized_hvg <- FindVariableFeatures(sc_normalized, selection.method = "vst", nfeatures = number_top_genes)

          # HVG volcano plots
          # top10 <- head(VariableFeatures(sc_normalized_hvg), 10)
          #
          # sc_normalized_hvg_volcano_plot <- VariableFeaturePlot(sc_normalized_hvg)
          # sc_normalized_hvg_volcano_plot_name <- paste0(sample, "_sc_normalized_hvg_volcano_plot.png")
          # ggplot2::ggsave(filename = sc_normalized_hvg_volcano_plot_name, plot = sc_normalized_hvg_volcano_plot, path = qc_dir, width = 12, height = 8, dpi = 300)
          #
          # sc_normalized_hvg_volcano_plot_top10_labelled_name <- paste0(sample, "_sc_normalized_hvg_volcano_plot_top10_labelled.png")
          # sc_normalized_hvg_volcano_plot_top10_labelled <- LabelPoints(plot = sc_normalized_hvg_volcano_plot_top10_labelled, points = top10, repel = TRUE)
          # ggplot2::ggsave(filename = sc_normalized_hvg_volcano_plot_top10_labelled_name, plot = sc_normalized_hvg_volcano_plot_top10_labelled, path = qc_dir, width = 12, height = 8, dpi = 300)

          # Scale data
          sc_normalized_scaled_hvg <- ScaleData(sc_normalized_hvg, features = all_genes)

          # Save as RDS
          saveRDS(sc_normalized_scaled_hvg, file = output_hvg_rds)
          write_log(paste0("Saved normalized and scaled HVG data for ", sample, " to ", output_hvg_rds), log_file)

          # Save as tsv
          active_assay <- DefaultAssay(sc_normalized_scaled_hvg)
          counts <- LayerData(sc, assay = active_assay, layer = "counts")
          sc_normalized_scaled_hvg_counts <- counts[sc_normalized_scaled_hvg, drop = FALSE]
          sc_normalized_scaled_hvg_counts_dense <- as.matrix(t(sc_normalized_scaled_hvg_counts))

          write.table(
            sc_normalized_scaled_hvg_counts_dense,
            file = sc_normalized_scaled_hvg_counts_tsv,
            sep = "\t",
            quote = FALSE,
            col.names = NA
          )
          write_log(paste0("Saved normalized and scaled HVG data for ", sample, " to ", output_hvg_tsv), log_file)

        } else {

          # Scale data
          sc_normalized_scaled <- ScaleData(sc_normalized, features = all_genes)

          # Save as RDS
          saveRDS(sc_normalized_scaled, file = output_rds)

          # Save as tsv
          active_assay <- DefaultAssay(sc)
          normalized_scaled_counts <- LayerData(sc_normalized_scaled, assay = active_assay, layer = "data")
          normalized_scaled_dense <- as.matrix(t(normalized_scaled_counts))
          write.table(
            normalized_scaled_dense,
            file = output_tsv,
            sep = "\t",
            quote = FALSE,
            col.names = NA
          )
        }
      }
      write_log("All sample normalizations completed successfully", log_file)
    },
    error = function(e) {
      msg <- paste0("ERROR CODE 13: Normalization failed due to", conditionMessage(e))
      write_log(msg, log_file)
    }
  )
}
