check.fastq.files <- function (input_folder) {

  # list out all the fastq and fastq.gz files
  fastq.files <- list.files(input_folder, pattern="fastq|fastq.gz", recursive=TRUE)

  # group output by 3
  # R doesn't have dictionaries -> nested lists? named lists?
  # could do split(x = data, f = groups)? groups = unique(grep("sample"))

  group_names <- sub("^([^_]+)_.*$", "\\1", fastq.files)
  groups <- split(fastq.files, group_names)

  # check if groups have 3 items

  for (l in split_lengths) {
    if (l != 3) {
      warning(paste0("WARNING: Sample ", l, " should have 3 FASTQ files, but found ", ))
    }
  }
}
