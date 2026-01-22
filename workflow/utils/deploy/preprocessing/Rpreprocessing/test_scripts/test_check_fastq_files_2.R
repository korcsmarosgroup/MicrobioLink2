repo_root <- normalizePath(
  file.path(testthat::test_path(), "..", "..", "..", "..", "..", ".."),
  mustWork = TRUE
)

script_dir <- file.path(repo_root, "workflow", "utils", "preprocessing", "Rpreprocessing", "R")

source(file.path(script_dir, "add_log_message.R"))
source(file.path(script_dir, "checking_fastq_files_2.R"))

create_temp_log <- function() {
  log_file <- tempfile("preprocessing_", fileext = ".log")
  file.create(log_file)
  log_file
}

test_that("invalid Fastq_file_format raises error code 5", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  expect_error(
    check_fastq_files(tmp_dir, "invalid", log_file),
    "ERROR CODE 5: Invalid Fastq_file_format"
  )
})

test_that("merged layout without FASTQ files throws error code 6", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  expect_error(
    check_fastq_files(tmp_dir, "merged", log_file),
    "ERROR CODE 6: No FASTQ files found"
  )
})

test_that("merged layout returns paired samples for R1/R2 files", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  file.create(file.path(tmp_dir, "SampleA_R1.fastq"))
  file.create(file.path(tmp_dir, "SampleA_R2.fastq"))

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  sample_map <- check_fastq_files(tmp_dir, "merged", log_file)

  expect_true("SampleA" %in% names(sample_map))
  expect_identical(sample_map[["SampleA"]][["R1"]], "SampleA_R1.fastq")
  expect_identical(sample_map[["SampleA"]][["R2"]], "SampleA_R2.fastq")
})

test_that("subdir layout without sample subfolders triggers error code 7", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  expect_error(
    check_fastq_files(tmp_dir, "subdir", log_file),
    "ERROR CODE 7: No sample subdirectories found"
  )
})

test_that("subdir layout warns when mates are missing", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  sample_dir <- file.path(tmp_dir, "SampleB")
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  file.create(file.path(sample_dir, "SampleB_R1.fastq"))

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  warning_messages <- character()

  withCallingHandlers(
    check_fastq_files(tmp_dir, "subdir", log_file),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_true(any(grepl("R2 file is missing", warning_messages)))

  log_entries <- readLines(log_file)
  expect_true(any(grepl("SampleB has FASTQ files", log_entries)))
})
