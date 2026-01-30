usethis::use_package("Matrix", type = "Imports")
usethis::use_package("Seurat", type = "Imports")
usethis::use_package("ggplot2", type = "Imports")
usethis::use_package("R.utils", type = "Imports")
usethis::use_package("DoubletFinder", type = "Imports")
#' QC_10x
#'
#' Perform comprehensive QC on 10x Genomics data:
#'    - Load raw & filtered matrices
#'    - Compute QC metrics
#'    - Visualize QC (violin & bar plots)
#'    - Filter poor-quality cells
#'    - Detect potential doublets using Scrublet
#'
#' @export

# output folder = sample1_Solo.out

QC_10x <- function(output_folder, log_file) {
  files_to_gzip <- list(
    c("matrix.mtx", "matrix.mtx.gz"),
    c("barcodes.tsv", "barcodes.tsv.gz"),
    c("features.tsv", "features.tsv.gz")
  )

  if (platform == "10x") {
    checking_directory <- file.path(output_folder, "STARsolo_10x_out")
  } else if (platform == "DropSeq") {
    checking_directory <- file.path(output_folder, "STARsolo_DropSeq_out")
  }

  dir_list <- list.dirs(checking_directory, recursive = FALSE, full.names = FALSE)
  for (sample in dir_list) {
    sample_path <- file.path(checking_directory, sample)
    raw_dir <- file.path(sample_path, "Solo.out", "Gene", "raw")
    filtered_dir <- file.path(sample_path, "Solo.out", "Gene", "filtered")

    if (!dir.exists(filtered_dir)) {
      write_log(paste0("Skipping ", sample, " : no filtered directory found."), log_file)
      next
    }

    for (i in seq_along(files_to_gzip)) {
      in_path <- file.path(filtered_dir, files_to_gzip[[i]][[1]])
      out_path <- file.path(filtered_dir, files_to_gzip[[i]][[2]])

      # If gz file exists already, skip
      if (file.exists(out_path)) {
        next
      }

      # If ungzipped file exists, compress it
      if (file.exists(in_path)) {
        write_log(paste0("Gzipping ", in_name, " for ", sample), log_file)
        tryCatch(
          {
            R.utils::gzip(
              filename  = in_path,
              destname  = out_path,
              overwrite = FALSE,
              remove    = FALSE
            )
          },
          error = function(e) {
            write_log(paste0("ERROR gzipping ", name, " for ", sample, " : ", conditionMessage(e)), log_file)
            next
          }
        )
      } else {
        write_log(paste0("ERROR: ", in_name, " and ", out_name, " both missing for ", sample), log_file)
        next
      }
    }

    # --- Load data ---
    load_10x <- function(filtered_dir) {
      tryCatch(
        {
          write_log(paste0("Loading filtered 10x data for ", sample), log_file)
          sc_droplet_data <- Seurat::Read10X(filtered_dir, gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
          sc_droplet <- Seurat::CreateSeuratObject(counts = sc_droplet_data)
          return(sc_droplet)
        },
        error = function(e) {
          write_log(paste0("ERROR loading data for ", sample, " :", conditionMessage(e)), log_file)
          return(paste0("ERROR loading data for ", sample, " :", (conditionMessage(e))))
        }
      )
    }

    result <- load_10x(filtered_dir)
    if (grepl("error", result, ignore.case = TRUE)) next

    # --- Basic QC metrics & selecting cells for further analysis ---

    write_log(paste0("Calculating QC metrics for ", sample), log_file)
    sc_droplet[["percent.mt"]] <- Seurat::PercentageFeatureSet(sc_droplet, patter = "^MT-")

    # --- Violin Plots ---

    qc_plot_dir <- file.path(sample_path, "QC_plots")
    dir.create(qc_plot_dir)
    write_log(paste0("Generating QC violin plots for ", sample), log_file)

    sc_droplet_VlnPlot <- VlnPlot(sc_droplet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    sc_droplet_VlnPlot_name <- paste0(sample, "_sc_droplet_VlnPlot.png")
    ggplot2::ggsave(filename = sc_droplet_VlnPlot_name, plot = sc_droplet_VlnPlot, path = qc_plot_dir, width = 12, height = 8, dpi = 300)

    # --- FeatureScatter: mitochondiral proportion ---

    sc_droplet_FeatureScatter <- FeatureScatter(sc_droplet, feature1 = "nCount_RNA", feature2 = "percent.mt")
    sc_droplet_FeatureScatter_filename <- paste0(sample, "_mt_FeatureScatter.png")
    ggplot2::ggsave(filename = sc_droplet_FeatureScatter_filename, plot = sc_droplet_FeatureScatter, path = qc_plot_dir, width = 6, height = 6, dpi = 300)
    write_log(paste0("Saved FeatureScatter to: ", qc_plot_dir), log_file)

    # --- Filtering poor-quality cells ---

    write_log(paste0("Filtering low-quality cells for ", sample), log_file)
    ncells_before <- ncol(sc_droplet)
    sc_droplet_filtered <- subset(sc_droplet, subset = nFeature_RNA < 2500 & nCount_RNA < 1000 & percent.mt < 5)
    ncells_after <- ncol(sc_droplet_filtered)
    nfiltered_cells <- ncells_before - ncells_after
    write_log(paste0(sample, ": filtered ", nfiltered_cells, " cells; ", ncells_after, " remain."), log_file)

    # Handle empty Seurat object after filtering
    if (isTRUE(ncol(sc_droplet_filtered) == 0)) {
      write_log(paste0(sample, ": all cells filtered out, skipping save."), log_file)
      next
    }

    if (isTRUE(ncol(sc_droplet_filtered) <= 50)) {
      write_log(paste0("Skipping doublet detection for ", sample, " too few cells after filtering (n=", ncol(sc_droplet_filtered)), write_log)
    } else {
      # --- Doublet detection with DoubletFinder ---

      sc_droplet_filtered <- tryCatch(
        {
          # Clustering

          sc_droplet_filtered <- FindNeighbors(sc_droplet_filtered, reduction = "pca", dims = 1:10)
          sc_droplet_filtered <- FindClusters(sc_droplet_filtered, resolution = 0.3)

          viz_dimloadings_plot <- VizDimLoadings(sc_droplet_filtered, dims = 1:2, reduction = "pca")
          dimplot_plot <- DimPlot(sc_droplet_filtered, group.by = "seurat_clusters") + NoLegend()

          viz_dimploadings_filename <- paste0(sample, "viz_dimloadings_plot.png")
          dimplot_plot_filename <- paste0(sample, "dimplot_plot.png")
          ggplot2::ggsave(viz_dimploadings_filename, viz_dimloadings_plot, path = qc_plot_dir, width = 6, height = 6, dpi = 300)
          ggplot2::ggsave(dimplot_plot_filename, dimplot_plot, path = qc_plot_dir, width = 6, height = 6, dpi = 300)

          # Predicted doublet rate

          multiplet_rates <- data.frame(
            "cells_recovered" = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000),
            "multiplet_rate" = c(0.4, 0.8, 1.6, 2.4, 3.2, 4.0, 4.8, 5.6, 6.4, 7.2, 8.0)
          )
          predicted_rates <- lm(multiplet_rate ~ cells_recovered, data = multiplet_rates)
          actual_recovered_cells <- ncol(pbmc10k_filt)
          cell_number <- data.frame(cells_recovered = c(actual_recovered_cells))
          predicted_rate <- predict(predicted_rates, newdata = cell_number)
          predicted_rate_doubletfinder <- predicted_rate / 100
          predicted_rate_doubletfinder <- round(predicted_rate_doubletfinder, digits = 2)
          # Doublet detection

          sweep.res.list_sc_droplet_filtered <- paramSweep(sc_droplet_filtered, PCs = 1:10, sct = FALSE)
          sweep.stats_sc_droplet_filtered <- summarizeSweep(sweep.res.list_sc_droplet_filtered, GT = FALSE)
          bcmvn_sc_droplet_filtered <- find.pK(sweep.stats_sc_droplet_filtered)
          bcmvn_sc_droplet_filtered$pK <- as.numeric(as.character(bcmvn_sc_droplet_filtered$pK))
          optimal_pK <- bcmvn_sc_droplet_filtered$pK[which.max(bcmvn_sc_droplet_filtered$BCmetric)]

          homotypic.prop <- modelHomotypic(sc_droplet_filtered@meta.data$seurat_clusters)
          nExp_poi <- round(predicted_rate_doubletfinder * nrow(sc_droplet_filtered@meta.data))
          nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

          sc_droplet_filtered <- doubletFinder(sc_droplet_filtered, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)

          pANN_col <- grep("pANN", colnames(sc_droplet_filtered@meta.data), value = TRUE)
          sc_droplet_filtered <- doubletFinder(sc_droplet_filtered, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj, reuse.pANN = pANN_col, sct = FALSE)

          DF_classifications <- grep("DF.classifications", colnames(sc_droplet_filtered@meta.data), value = TRUE)
          DF_classification <- DF_classifications[-1]

          doublet_calls <- sc_droplet_filtered@meta.data[[DF_classification]]
          call_frequencies <- table(doublet_calls)
          write_log(paste0(sample, ": estimated doublet rate = ", call_frequencies[[1]]), log_file)

          sc_droplet_filtered <- sc_droplet_filtered[
            , sc_droplet_filtered@meta.data[[DF_classification]] == "Singlet"
          ]
        },
        warning = function(w) {
          write_log(paste0("WARNING: DoubletFinder failed for ", sample, conditionMessage(w)), write_log)
        }
      )

      # --- Save QC results ---

      qc_out <- file.path(sample_path, "QC_filtered.rds")
      saveRDS(sc_droplet_filtered, file = qc_out)
      write_log(paste0("QC complete for ", sample, ". Saved to ", qc_out, "."), log_file)
    }
  }
}
