#' checking_fastq_files_2.R
#'
#' Checking the FASTQ input folder: verifies that all expected sequencing files are present
#' and correctly structured according to the specified Fastq_file_format ('merged' or 'subdir').
#'
#' If any expected file (R1, R2, or I1) is missing, or the folder structure does not match
#' the specified layout, an error is raised, and execution stops.
#'
#' @param input_dir The system-based absolute path to the input directory containing the FASTQ files or pe
#' -sample subdirectories.
#' @param Fastq_file_format The organization style of the FASTQ data to validate.
#'    Must be either:
#'        - 'merged': all FASTQ files are in a single directory.
#'        - 'subdir': each sample has its own subdirectory containing FASTQ files.
#' @param log_file The path of the log file.
#'
#' @details
#' Error codes:
#'    ERROR CODE 5: Invalid Fastq_file_format
#'    ERROR CODE 6: No FASTQ files found (merged layout)
#'    ERROR CODE 7: Subdirectory missing R1/R2/I1 (subdir layout)
#'    ERROR CODE 8: Subdirectory has wrong number of files (subdir layout)
#'
#' @return NULL
#'    The function does not return any value. If all checks pass, execution continues normally.
#'
#' @export


check_fastq_files <- function(input_dir, Fastq_file_format, log_file) {

  # ------------------------------------------------------------
  # Validate Fastq_file_format
  # ------------------------------------------------------------

  if (!(Fastq_file_format %in% c("merged", "subdir"))) stop("ERROR CODE 5: Invalid Fastq_file_format. Must be      merged or subdir.")

  #
  # --------------------------------------------------------
  #   Case 1: Flat layout - all FASTQs in one directory
  # --------------------------------------------------------
  #

  if (Fastq_file_format == "merged") {
    fastq.files <- list.files(input_dir, pattern = "fastq|fastq.gz", recursive = FALSE)
    if (length(fastq.files) == 0) {
      stop(paste0("ERROR CODE 6: No FASTQ files found in ", input_dir))
    }

    tags <- c("R1", "R2", "I1")

    group_names <- sub("^([^_]+)_.*$", "\\1", fastq.files)
    group_names <- unique(group_names)

    samples <- list()
    for (g in group_names) {
      samples[[g]] <- list(
        R1 = NULL,
        R2 = NULL,
        I1 = NULL
      )
    }

    for (sample in names(samples)) {
      pattern <- sample
      for (tag in (names(samples[[sample]]))) {
        for (f in fastq.files) {
          if (grepl(sample, f) && grepl(tag, f)) {
            samples[[sample]][[tag]] <- f
          }
        }
      }
    }

    for (tag in names(samples[[sample]])) {
      value <- samples[[sample]][[tag]]
      if (is.null(value) || is.na(value)) {
        message(sample, " ", tag, " file is missing")
      }
    }

    for (sample in names(samples)) {
      missing_tags <- names(samples[[sample]])[sapply(samples[[sample]], is.null)]
      if (length(missing_tags) > 0) {
        message(paste0(sample, " should have 3 FASTQ files. But only ", 3 - (length(missing_tags)), " files found."))
      }
    }
  }

  # ------------------------------------------------------------
  # Case 2: Subdirectory layout — one folder per sample
  # ------------------------------------------------------------

  else if (Fastq_file_format == "subdir") {
    subdirs <- list.dirs(input_dir, recursive = FALSE)

    if (length(subdirs) == 0) {
      stop(paste0("ERROR CODE 7: No sample subdirectories found in ", subdirs))
    }

    write_log("Using per-sample subdirectory layout", log_file)

    for (folder in subdirs) {
      sample_name <- basename(folder)
      fastqs <- list.files(folder, pattern = "fastq|fastq.gz", recursive = FALSE)
      if (length(fastqs) == 0) {
        write_log(paste0("WARNING: ", sample_name, " has no FASTQ files"), log_file)
        # TODO: Check how can we implement a python continue-equivalent method here
        # Aim: Not stop the whole script, just jump to the next folder in the for loop
      } else {
        write_log(paste0(sample_name, " has FASTQ files"), log_file)
      }

      group_names <- sub("^([^_]+)_.*$", "\\1", fastqs)
      group_names <- unique(group_names)

      samples <- list()
      for (g in group_names) {
        samples[[g]] <- list(
          R1 = NULL,
          R2 = NULL,
          I1 = NULL
        )
      }

      for (sample in names(samples)) {
        pattern <- sample
        for (tag in (names(samples[[sample]]))) {
          for (f in fastqs) {
            if (grepl(sample, f) && grepl(tag, f)) {
              samples[[sample]][[tag]] <- f
            }
          }
        }
      }

      for (tag in names(samples[[sample]])) {
        value <- samples[[sample]][[tag]]
        if (is.null(value) || is.na(value)) {
          warning(sample, " ", tag, " file is missing")
          write_log(paste0("WARNING: ", sample, " is missing", tag, " file"), log_file)
        }
      }

      for (sample in names(samples)) {
        missing_tags <- names(samples[[sample]])[sapply(samples[[sample]], is.null)]
        if (length(missing_tags) > 0) {
          message(paste0(
            sample, " should have 3 FASTQ files. But only ", 3 - (length(missing_tags)),
            " files found."
          ))
        }
      }
      write_log(paste0(sample_name, " has the necessary FASTQ files (R1, R2, R3"), log_file)
    }
  }
}
