#' run_star_solo_unified
#'
#' Run STAR or STARsolo depending on the experimental platform and configuration settings.
#' This function executes the STAR aligner for microwell-based experiments or STARsolo
#' for droplet-based single-cell RNA-seq data. It automatically detects input FASTQ pairs,
#' constructs appropriate STAR command-line arguments, and logs progress and errors.
#' The function supports both gzipped and uncompressed FASTQ files, and dynamically
#' includes custom STAR/STARsolo parameters defined in the configuration.
#'
#' @param configuration (dict): Dictionary containing run parameters and paths. Expected keys:
#'  - "platform" (str): Sequencing platform ("microwell" or "droplet").
#'  - "input_dir" (str): Directory containing FASTQ files.
#'  - "genome_dir" (str): Directory of the STAR genome index.
#'  - "threads" (int, optional): Number of CPU threads for STAR. Defaults to 8.
#'  - "Fastq_file_format" (str, optional): Layout of FASTQ files ("merged" or "subdir").
#'  - "STAR_outdir" (str, optional): Output directory for STAR alignments.
#'  - "STAR_params" (dict, optional): Additional STAR command-line parameters.
#'  - "STARsolo_outdir" (str, optional): Output directory for STARsolo.
#'  - "STARsolo_params" (dict, optional): Additional STARsolo command-line parameters.
#'
#'  @param log_file (str): Path to a log file for recording messages and errors.
#'
#'  @param solo (bool, optional):
#'  - If True: Runs STARsolo (droplet-based single-cell data).
#'  - If False: Runs standard STAR alignment (microwell-based data).
#'  Defaults to False.
#'
#'  @return Returns:
#'
#'  Notes:
#'  - Requires STAR to be installed and accessible in the system PATH.
#'  - FASTQ files must contain '_R1' and '_R2' in their filenames.
#'  - Skips samples missing either R1 or R2.
#'  - Use
#'
#'  @export

run_STAR_unified <- function(configuration, log_file, solo = FALSE) {
  platform <- tolower(file$platform)
  input_dir <- file$input_dir
  genome_dir <- file$genome_dir
  threads <- file$threads
  layout <- file$Fastq_file_format

  if (solo == TRUE) {
    if (platform != "droplet") {
      write_log(paste0("Platform is not droplet: ", platform, ". Skipping STARsolo"), log_file)
      output_dir <- paste0(file$STARsolo_outdir, input_dir, file$STARsolo_out, sep = "/")
      params <- file$STARsolo_params
    } else {
      write_log(paste0("Platform is not droplet: ", platform, ". Skipping STAR"), log_file)
      output_dir <- paste0(file$STAR_outdir, input_dir, file$STAR_out, sep = "/")
      params <- file$STAR_params
    }
  }

  samples <- pair_fastqs(input_dir, layout)

  for (sample in names(samples)) {
    files <- sample[[samples]]
    found <- sapply(tags, function(tag) any(grepl(tag, files)))
    if (!(all(found))) {
      missing_tags <- tags[!found]
      message <- message(sample, ": missing -> ", paste0(missing_tags, collapse = ", "))
      write_log(message, log_file)
    }

    pattern_1 <- paste0(sample, "R1")
    pattern_2 <- paste0(sample, "R2")

    cmd <- "STAR"
    args <- c("--runThreadN", as.character(threads),
              "--genomeDir", genome_dir,
              "--readFilesIn", all(grep(pattern_1, files, values = TRUE)), all(grep(pattern_2, files, values = TRUE)),
              "--outFileNamePrefix", file.path(sample_out, ""))

    # gzipped FASTQs

    # Add STAR/STARsolo parameters

    result <- tryCatch({
      ifelse(solo, mode = "STARsolo", mode = "STAR")
    } error = function(e) cat("ERROR: ", mode, " failed for ", sample))

  }
}
