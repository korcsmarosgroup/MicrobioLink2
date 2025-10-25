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

      group_names <- sub("^([^_]+)_.*$", "\\1", fastq.files)
      groups <- split(fastq.files, group_names)
      sample_order <- order(as.numeric(gsub("\\D", "", names(groups))))
      groups <- groups[sample_order]
      read_order <- c("R1", "R2", "I1")
      groups <- lapply(groups, function(files) {
        files[order(match(sub(".*_(R1|R2|I1).*", "\\1", files), read_order))]
      })

      pattern <- c("_R1", "_R2", "_I1")
      check_files <- function(groups, pattern) {
          for (i in seq_along(groups)) {
              if (all(sapply(pattern, function(x) any(grepl(x, groups[[i]]))))) {
                  message(sprintf("Sample %d contains all required FASTQ files", i))
                } else {
                  warning(sprintf("Sample %d is missing one or more FASTQ files", i))
                }
          }
      }

      print(check_files(groups, pattern))

      lapply(groups, function(x)
          if (length(x) != 3) write_log("Sample should have 3 FASTQ files")
        )

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

      group_names <- sub("^([^_]+)_.*$", "\\1", fastq.files)
      groups <- split(fastq.files, group_names)
      sample_order <- order(as.numeric(gsub("\\D", "", names(groups))))
      groups <- groups[sample_order]
      read_order <- c("R1", "R2", "I1")
      groups <- lapply(groups, function(files) {
          files[order(match(sub(".*_(R1|R2|I1).*", "\\1", files), read_order))]
      })

      pattern <- c("_R1", "_R2", "_I1")
      check_files <- function(groups, pattern) {
          for (i in seq_along(groups)) {
              if (all(sapply(pattern, function(x) any(grepl(x, groups[[i]]))))) {
                  message(sprintf("Sample %d contains all required FASTQ files", i))
                } else {
                  warning(sprintf("Sample %d is missing one or more FASTQ files", i))
                }
          }
      }

      lapply(groups, function(x)
          if (length(x) != 3) write_log("Sample should have 3 FASTQ files")
      )

    }

}
