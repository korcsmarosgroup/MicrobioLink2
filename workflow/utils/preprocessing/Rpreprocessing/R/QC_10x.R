usethis::use_package("Matrix", type = "Imports")
usethis::use_package("Seurat", type = "Imports")
#'  QC_10x
#'
#'  Perform comprehensive QC on 10x Genomics data:
#'    - Load raw & filtered matrices
#'    - Compute QC metrics
#'    - Visualize QC (violin & bar plots)
#'    - Filter poor-quality cells
#'    - Detect potential doublets using Scrublet
#'
#'  @export

# output folder = sample1_Solo.out

QC_10x <- function(output_folder, log_file) {
  # sample1, sample2, sample3, ...
  dir_list <- list.dirs(output_folder, recursive = FALSE, full.names = FALSE)
  # sample1
  for (sample in dir_list) {
    # output_folder/sample1
    sample_path <- file.path(ouput_folder, sample)
    # output_folder/sampl1/raw
    raw_dir <- file.path(sample_path, "raw")
    # output_folder/sample_1/filtered
    filtered_dir <- file.path(sample_path, "filtered")

    if (!dir.exists(filtered_dir)) {
      write_log(paste0("Skipping ", sample, " : no filtered directory found."), log_file)
    }

    # --- Load data ---

    load_10x <- function(filtered_dir) {
      tryCatch(
        {
          sc10x_data <- Seurat::Read10X(filtered_dir, gene.column = 2, cell.column = 1, unique.features = TRUE, strip.suffix = FALSE)
          sc10x <- Seurat::CreateSeuratObject(counts = sc10x_data)
          return(sc10x)
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
    sc10x[["percent.mt"]] <- Seurat::PercentageFeatureSet(sc10x, patter = "^MT-")

    # --- Violin Plots ---

    qc_plot_dir <- file.path(sample_path, "QC_plots")
    dir.create(qc_plot_dir)
    write_log(paste0("Generating QC violin plots for ", sample), log_file)

    sc10x_VlnPlot <- VlnPlot(sc10x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(sc10x_VlnPlot)
    ggplot2::ggsave(filename = "sc10x_VlnPlot", plot = sc10x_VlnPlot, path = qc_plot_dir, width = 12, height = 8, dpi = 300)

    # --- Filtering poor quality cells ---

    write_log(paste0("Filtering low quality cells for ", sample), log_file)
    ncells_before <- ncol(sc10x)
    sc10x <- subset(sc10x, subset = nFeature_RNA < 2500 & nCount_RNA < 1000 & percent.mt < 5)
    ncells_after <- ncol(sc10x)
    nfiltered_cells <- ncells_before - ncells_after
    write_log(paste0(sample, " : filtered ", nfiltered_cells, " cells; ", ncells_after, " remain."), log_file)

    # Handle empty Seurat object after filtering
    if (ncol(sc10x) == 0) {
      write_log(paste0(sample, " : all cells filtered out, skipping save."), log_file)
      next
    }

    # --- Doublet detection with DoubletFinder ---

    # Save QC results

    qc_out <- file.path(sample_path, "QC_filtered.rds")
    saveRDS(sc10x, file = qc_out)
    write_log(paste0("QC complete for ", sample, ". Saved to ", qc_out, "."), log_file)
  }
}
