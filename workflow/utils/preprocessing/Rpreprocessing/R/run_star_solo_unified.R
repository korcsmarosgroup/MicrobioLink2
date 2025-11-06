usethis::use_package("processx", type = "Import")

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
  stop("function not implemented yet")

  platform <- tolower(file$platform)
  input_dir <- file$input_dir
  genome_dir <- file$genome_dir
  threads <- as.character(file$threads)
  layout <- file$Fastq_file_format

  if (isTRUE(solo == TRUE)) {
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
    files <- samples[[sample]]
    found <- sapply(tags, function(tag) any(grepl(tag, files)))
    if (!(all(found))) {
      missing_tags <- tags[!found]
      write_log(paste0("Skipping ", sample, " missing R1 or R2"), log_file)
    }

    sample_out <- file.path(output_dir, sample)
    dir.create(sample_out)

    file_R1 <- grep("R1", files, value = TRUE)
    file_R2 <- grep("R2", files, value = TRUE)

    cmd_base <- "STAR"
    cmd_args <- list(
      paste("--runThreadN", shQuote(threads)),
      paste("--genomeDir", shQuote(genome_dir)),
      paste("--readFilesIn", shQuote(file_R1), shQuote(file_R2)),
      paste("--outFileNamePrefix", shQuote(sample_out))
    )

    cmd <- c(cmd_base, cmd_args)
    pipeline <- paste(cmd, collapse = " ")

    # gzipped FASTQs

    gzipped <- FALSE
    zcat <- NULL
    cmd_unzip <- NULL

    if (isTRUE(all(any(grepl(".gz$", files))) && any(grepl("R[12]", files)))) {
      gzipped <- TRUE
      zcat <- "zcat"
      cmd_unzip <- list(
        paste("--readFilesCommand", shQuote(zcat))
      )
      cmd <- c(cmd, cmd_unzipped)
      pipeline <- paste(cmd, collapse = " ")
    }

    for (param in names(params)) {
      value <- params[[param]]
      if (!is.null(value)) {
        cmd_new <- paste("--", param, shQuote(value))
        cmd <- c(cmd, cmd_new)
        pipeline <- paste(cmd, collapse = " ")
      }
    }

    STARUnified_result <- tryCatch(
      {
        if (isTRUE(all(any(grepl(".gz$", files))) && any(grepl("R[12]", files)))) {
          if (isTRUE(mode == "STARsolo")) {
            write_log(paste0("Running ", mode, "for ", sample, " : ", pipeline), log_file)
            result <- processx::run(cmd_base, c(unlist(cmd_args), cmd_unzip), error_on_status = TRUE)
            write_log(paste0(mode, " completed for ", sample, " . Output: ", sample_out), log_file)
            return(result)
          } else if (isTRUE(mode == "STAR")) {
            write_log(paste0("Running ", mode, "for ", sample, " : ", pipeline), log_file)
            result <- processx::run(cmd_base, c(unlist(cmd_args), cmd_unzip), error_on_status = TRUE)
            write_log(paste0(mode, " completed for ", sample, " . Output: ", sample_out), log_file)
            return(result)
          }
        } else {
          if (isTRUE(mode == "STARsolo")) {
            write_log(paste0("Running ", mode, "for ", sample, " : ", pipeline), log_file)
            result <- processx::run(cmd_base, unlist(cmd_args), error_on_status = TRUE)
            write_log(paste0(mode, " completed for ", sample, " . Output: ", sample_out), log_file)
            return(result)
          } else if (isTRUE(mode == "STAR")) {
            write_log(paste0("Running ", mode, "for ", sample, " : ", pipeline), log_file)
            result <- processx::run(cmd_base, unlist(cmd_args), error_on_status = TRUE)
            write_log(paste0(mode, " completed for ", sample, " . Output: ", sample_out), log_file)
            return(result)
          }
        }
      },
      error = function(e) {
        write_log(paste0("ERROR: ", mode, " failed for: ", sample, result$stderr))
        cat("ERROR MESSAGE", mode, "failed for:", sample, e$message)
        return(NULL)
      }
    )
  }
}
