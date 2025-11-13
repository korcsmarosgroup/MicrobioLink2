source("./R/pair_fastqs.R")

test_that("pair_fastqs pairs and orders FASTQs in merged layout", {
  input_dir <- file.path(tempdir(), "pair_fastqs_merged")
  dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(input_dir, recursive = TRUE), add = TRUE)

  files_to_create <- file.path(
    input_dir,
    c(
      "sample1_R2.fastq",
      "sample1_R1.fastq",
      "sample1_I1.fastq",
      "sample2_R2.fastq.gz",
      "sample2_R1.fastq.gz"
    )
  )
  file.create(files_to_create)

  result <- pair_fastqs(input_dir, "merged")
  result_names <- basename(names(result))
  result <- stats::setNames(result, result_names)

  expect_equal(sort(result_names), c("sample1", "sample2"))
  expect_equal(
    basename(result[["sample1"]]),
    c("sample1_R1.fastq", "sample1_R2.fastq", "sample1_I1.fastq")
  )
  expect_equal(
    basename(result[["sample2"]]),
    c("sample2_R1.fastq.gz", "sample2_R2.fastq.gz")
  )
})

test_that("pair_fastqs pairs FASTQs when samples are stored in subdirectories", {
  input_dir <- file.path(tempdir(), "pair_fastqs_subdir")
  dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(input_dir, recursive = TRUE), add = TRUE)

  sample_dirs <- file.path(input_dir, c("sample1", "sample2"))
  lapply(sample_dirs, dir.create, showWarnings = FALSE)

  file.create(file.path(sample_dirs[[1]], c("sample1_R2.fastq", "sample1_R1.fastq")))
  file.create(file.path(sample_dirs[[2]], c("sample2_R2.fastq", "sample2_R1.fastq")))

  result <- pair_fastqs(input_dir, "subdir")

  expect_equal(names(result), c("sample1", "sample2"))
  expect_equal(
    result[["sample1"]],
    c("sample1_R1.fastq", "sample1_R2.fastq")
  )
  expect_equal(
    result[["sample2"]],
    c("sample2_R1.fastq", "sample2_R2.fastq")
  )
})
