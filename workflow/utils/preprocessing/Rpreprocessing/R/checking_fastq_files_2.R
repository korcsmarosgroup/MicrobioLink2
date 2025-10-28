check.fastq.files <- function (input_folder, Fastq_file_format) {

    stop("function not implemented yet")

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
      fastq.files <- list.files(input_folder, pattern="*.fastq$|*.fastq.gz$", recursive=FALSE)
      if (length(fastq.files) == 0) {
          stop(paste0("ERROR CODE 6: No FASTQ files found in ", input_folder))
      }

      write_log("Using flat FASTQ layout", log_file)

      # Validate file names

      tags <- c("R1", "R2", "I1")
      validate_file_names <- function(fastq.files, tags) {
        matches <- sapply(tags, function(t) grepl(t, fastq.files))
        if (any(rowSums(matches)) > 0) {
          warning("Unrecognised file name")
        }
      }

      # Group files by sample prefix

      fastq.files <- list.files(input_folder, pattern="fastq|fastq.gz", recursive=FALSE)
      group_names <- sub("^(.*)_(R1|R2|I1).*", "\\1", fastq.files)
      keys <- (unique(group_names))

      group_files <- function(keys, tags, fastq.files) {
        samples <- list()

        for (key in keys) {
          samples[[key]] <- list()
          for (tag in tags) {
            matched <- grep(paste0(key, ".*", tag), fastq.files, value=TRUE)
            samples[[key]][[tag]] <- matched
          }
        }

        return(samples)
      }

      # lapply(groups, function(x)
      #     if (length(x) != 3) write_log("Sample should have 3 FASTQ files")
      #   )

    }

  # ------------------------------------------------------------
  # Case 2: Subdirectory layout — one folder per sample
  # ------------------------------------------------------------

  else if (Fastq_file_format == "subdir") {
      subdirectories <- list.dirs(input_folder, full.names = TRUE, recursive = FALSE)
      fastq.files <- list.files(subdirectories, pattern="*.fastq$|*.fastq.gz$", recursive = FALSE)
#      if (length(fastq.files) == 0) {
#          stop(paste0("ERROR CODE 7: No FASTQ files found in ", input_folder))
#      }

      write_log("Using per-sample subdirectory layout", log_file )

      # Validate file names

      tags <- c("R1", "R2", "I1")
      validate_file_names <- function(fastq.files, tags) {
        matches <- sapply(tags, function(t) grepl(t, fastq.files))
        if (any(rowSums(matches)) > 0) {
          warning("Unrecognised file name")
        }
      }

      # Group files by sample prefix

      fastq.files <- list.files(input_folder, pattern="fastq|fastq.gz", recursive=FALSE)
      group_names <- sub("^(.*)_(R1|R2|I1).*", "\\1", fastq.files)
      keys <- (unique(group_names))

      group_files <- function(keys, tags, fastq.files) {
        samples <- list()

        for (key in keys) {
          samples[[key]] <- list()
          for (tag in tags) {
            matched <- grep(paste0(key, ".*", tag), fastq.files, value=TRUE)
            samples[[key]][[tag]] <- matched
          }
        }

        return(samples)
      }
      lapply(groups, function(x)
          if (length(x) != 3) write_log("Sample should have 3 FASTQ files")
      )

    }

}
