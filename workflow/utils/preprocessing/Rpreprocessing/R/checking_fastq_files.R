check.fastq.files <- function (input_folder) {

  # list out all the fastq and fastq.gz files
  fastq.files <- list.files(".", pattern="fastq|fastq.gz", recursive=TRUE)

  # group output by 3

  group_names <- sub("^([^_]+)_.*$", "\\1", fastq.files)
  groups <- split(fastq.files, group_names)

  # check if groups have 3 items

  for (l in split_lengths) {
    if (l != 3) {
      warning(paste0("WARNING: Sample ", l, " should have 3 FASTQ files, but found ", ))
    }
  }
}

group_names <- sub("^([^_]+)_.*$", "\\1", fastq.files)

groups <- split(fastq.files, group_names)

sample_order <- order(as.numeric(gsub("\\D", "", names(groups))))

groups <- groups[sample_order]

read_order <- c("R1", "R2", "I1")

groups <- lapply(groups, function(files) {

  files[order(match(sub(".*_(R1|R2|I1).*", "\\1", files), read_order))]

})
