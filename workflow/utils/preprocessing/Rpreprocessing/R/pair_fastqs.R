#' pair_fastqs
#'
#' Match R1, R2, and optional I1 FASTQ files per sample.
#'
#' This helper function scans the input directory (or per-sample subdirectories)
#' and pairs sequencing files for each sample based on filename prefixes. It works
#' for both 'merged' (all files in one directory) and 'subdir' (one folder per sample) layouts.
#'
#' Notes:
#' - Expects FASTQ filenames to contain '_R1', '_R2', and optionally '_I1'.
#' - Does not validate that all paired files exist; just groups them.
#' - Returns a dictionary suitable for feeding into STAR/STARsolo commands.
#'
#' @param input_dir (is.character) Absolute or relative path to the input FASTQ folder.
#' @param layout (is.character) Organization style of FASTQ files:
#' - 'merged': all FASTQs in a single directory.
#' - 'subdir': one subdirectory per sample containing FASTQs.
#'
#' @return dict: Nested list of samples and their corresponding FASTQ paths.
#'
#' @keywords internal

pair_fastqs <- function(input_dir, layout) {
  if (isTRUE(layout == "merged")) {
    fastq.files <- list.files(input_folder, pattern = "*.fastq$|*.fastq.gz$", full.names = TRUE, recursive = FALSE)
  } else if (isTRUE(layout == "subdir")) {
    subdirs <- list.dirs(input_folder, recursive = FALSE)
    fastq.files <- lapply(subdirs, function(sd) list.files(sd, pattern = ".R"))
    fastq.files <- unlist(fastq.files)
  }
  print(fastq.files)

  sample_names <- sub("^([^_]+)_.*$", "\\1", fastq.files)
  samples <- split(fastq.files, sample_names)
  sample_order <- order(as.numeric(gsub("\\D", "", names(samples))))
  samples <- samples[sample_order]
  read_order <- c("R1", "R2", "I1")
  samples <- lapply(samples, function(files) {
    files[order(match(sub(".*_(R1|R2|I1).*", "\\1", files), read_order))]
  })
  print(samples)
}
