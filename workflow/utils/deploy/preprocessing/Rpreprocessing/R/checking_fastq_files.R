#' checking_fastq_files
#'
#' Validate FASTQ inputs and return only samples with complete paired-end data.
#' The function scans the inpu directory that contains either:
#' -  A flat collection of FASTQ files ("merged" layout), or
#' -  Multiple per-sample subdirectories ("subdir" layout).
#'
#' The function:
#' -  Identifies FASTQ files for each sample (R1, R2, and optionally I1).
#' -  Logs warnings when samples are incomplete (missing R1 or R2)
#' -  Logs errors and exits when the folder structure is invalid.
#' -  Produces a named list containing **only samples that have both R1 and R2 FASTQs**
#'
#' The returns named list is suitable for downstream tools that require complete paired FASTQ inputs (e.g. STAR, trimming tools, quantification).
#'
#' @param input_dir The system-based absolute path to the input directory containing the FASTQ files or sample folders.
#' @param Fastq_file_format Expected layout of the FASTQ input. Must be:
#'        - 'merged': all FASTQ files are in a single directory.
#'        - 'subdir': each sample has its own subdirectory containing FASTQ files.
#' @param log_file Path to the log file where warnings and information messages are written.
#'
#' @details
#' Error codes:
#'    ERROR CODE 5: Invalid Fastq_file_format argument.
#'    ERROR CODE 6: No FASTQ files found in merged layout.
#'    ERROR CODE 7: No sample subdirectories found in subdir layout.
#'
#' @return named list: A named list mapping sample prefixes to nested lists of file paths, e.g.:
#' $sample1
#' $R1
#' $"/path/to/SampleA_R1.fastq.gz"
#' $R2
#' $"/path/to/SampleA_R2.fastq.gz"
#'
#' $sample2
#' $R1
#' $...
#' $R2
#' $...
#'
#' Only samples that contain **both R1 and R2** files are included.
#'
#' @export
check_fastq_files <- function(input_dir, Fastq_file_format, log_file) {
  # ------------------------------------------------------------
  # Validate Fastq_file_format
  # ------------------------------------------------------------

  if (!(Fastq_file_format %in% c("merged", "subdir"))) stop("ERROR CODE 5: Invalid Fastq_file_format. Must be merged or subdir.")

  #
  # --------------------------------------------------------
  #   Case 1: MERGED layout - all FASTQs in one directory
  # --------------------------------------------------------
  #

  if (Fastq_file_format == "merged") {
    fastq.files <- list.files(input_dir, pattern = "fastq(\\.gz)?$", recursive = FALSE)
    if (length(fastq.files) == 0) {
      stop(paste0("ERROR CODE 6: No FASTQ files found in ", input_dir))
    }

    write_log("Using merged FASTQ layout", log_file)

    group_names <- sub("[._-]?R[12]\\.fastq(\\.gz)?$", "", fastq.files, perl = TRUE)
    group_names <- sub("[._-]$", "", group_names)
    group_names <- unique(group_names)

    samples <- list()

    for (g in group_names) {
      for (file in fastq.files) {
        if (grepl("(^|[._-])R1([._-]|$)", file) && (grepl("(^|[._-])R1([._-]|$)", file))) {
          samples[[g]][["R1"]] <- file
        } else if ((grepl("(^|[._-])R2([._-]|$)", file)) && (grepl("(^|[._-])R2([._-]|$)", file))) {
          samples[[g]][["R2"]] <- file
        } else {
          warning(paste0("Skipping unrecognised file name: ", f))
          write_log(paste0("WARNING: Skipping unrecognised file name", f), log_file)
          next
        }
      }
    }

    for (sample in names(samples)) {
      missing_tags <- c()
      if (!("R1" %in% names(samples[[sample]]))) {
        missing_tags <- paste0(sample, "_R1")
        warning(sample, " R1 file is missing")
      }
      if (!"R2" %in% names(samples[[sample]])) {
        missing_tags <- paste0(sample, "_R2")
        warning(sample, " R2 file is missing")
      }
      if (length(missing_tags) > 0) {
        warning(paste0(
          sample, " should have 2 FASTQ files. But only ", 2 - (length(missing_tags)), " files found."
        ))
      }
    }

    return(samples)
  }

  # ------------------------------------------------------------
  # Case 2: SUBDIR layout — one folder per sample
  # ------------------------------------------------------------

  else if (Fastq_file_format == "subdir") {
    subdirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)

    if (length(subdirs) == 0) {
      stop(paste0("ERROR CODE 7: No sample subdirectories found in ", input_dir))
    }

    write_log("Using per-sample subdirectory layout.", log_file)

    samples <- list()
    for (folder in subdirs) {
      base_name <- basename(folder)
      fastqs <- list.files(folder, pattern = "fastq(\\.gz)?$", recursive = FALSE)
      if (length(fastqs) == 0) {
        write_log(paste0("WARNING: ", base_name, " has no FASTQ files"), log_file)
        next
      } else {
        write_log(paste0("WARNING: ", base_name, " has FASTQ files"), log_file)
        group_name <- sub("[._-]?R[12]\\.fastq(\\.gz)?$", "", fastqs, perl = TRUE)
        group_name <- sub("[._-]$", "", group_name)
        group_name <- unique(group_name)

        samples_entry <- list()
        for (f in fastqs) {
          if (grepl("(^|[._-])R1([._-]|$)", f)) {
            if (!("R1" %in% names(samples_entry))) {
              samples_entry[["R1"]] <- f
            } else {
              stop(paste0("ERROR: Multiple R1 files found in folder ", folder))
            }
          } else if (grepl("(^|[._-])R2([._-]|$)", f)) {
            if (!("R2" %in% names(samples_entry))) {
              samples_entry[["R2"]] <- f
            } else {
              stop(paste0("ERROR: Multiple R2 files found in folder ", folder))
            }
          } else {
            warning(paste0("Skipping unrecognised file name: ", f))
            write_log(paste0("WARNING: Skipping unrecognised file name", f), log_file)
            next
          }
        }
        samples[[group_name]] <- samples_entry
      }
    }

    for (sample in names(samples)) {
      missing_tags <- c()
      if (!("R1" %in% names(samples[[sample]]))) {
        missing_tags <- paste0(sample, "_R1")
        warning(sample, " R1 file is missing")
      }
      if (!"R2" %in% names(samples[[sample]])) {
        missing_tags <- paste0(sample, "_R2")
        warning(sample, " R2 file is missing")
      }
      if (length(missing_tags) > 0) {
        warning(paste0(
          sample, " should have 2 FASTQ files. But only ", 2 - (length(missing_tags)), " files found."
        ))
      } else {
        write_log(paste0(sample, " has the necessary FASTQ files (R1, R2)"))
      }
    }
    return(samples)
  }
  # ------------------------------------------------------------
  # FINAL STEP - RETURN ONLY SAMPLES WITH BOTH R1 AND R2
  # ------------------------------------------------------------

  required <- c("R1", "R2")
  missing <- c()

  for (s in names(samples)) {
    sample <- samples[[s]]
    if (!(required %in% names(sample))) {
      missing <- c(missing, s)
    } else {
      next
    }
  }

  complete_samples <- samples[!names(samples) %in% missing]
  print(complete_samples)

  return(complete_samples)
}

