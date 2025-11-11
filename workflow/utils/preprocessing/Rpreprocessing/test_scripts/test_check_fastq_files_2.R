library(testthat)

repo_root <- getwd()
script_path <- file.path(repo_root, "workflow", "utils", "preprocessing", "Rpreprocessing", "R", "checking_fastq_files_2.R")
if (!file.exists(script_path)) {
  repo_root <- normalizePath(file.path(repo_root, "..", "..", "..", "..", ".."), mustWork = TRUE)
  script_path <- file.path(repo_root, "workflow", "utils", "preprocessing", "Rpreprocessing", "R", "checking_fastq_files_2.R")
}

source(file.path(repo_root, "workflow", "utils", "preprocessing", "Rpreprocessing", "R", "add_log_message.R"))
source(script_path)

original_write_log <- write_log

create_temp_log <- function() {
  log_file <- tempfile("preprocessing_", fileext = ".log")
  file.create(log_file)
  log_file
}

test_that("check_fastq_files rejects invalid Fastq_file_format", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  expect_error(
    check_fastq_files(tmp_dir, "unsupported", log_file),
    "ERROR CODE 5: Invalid Fastq_file_format"
  )
})

test_that("merged layout without FASTQ files triggers error code 6", {
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

test_that("merged layout reports missing tag via messages", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  file.create(file.path(tmp_dir, "SampleA_R1.fastq"))
  file.create(file.path(tmp_dir, "SampleA_R2.fastq"))

  expect_message(
    check_fastq_files(tmp_dir, "merged", log_file),
    "SampleA I1 file is missing"
  )
})

test_that("subdir layout with no sample folders errors early", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  expect_error(
    check_fastq_files(tmp_dir, "subdir", log_file),
    "ERROR CODE 7: No sample subdirectories found in"
  )
})

test_that("subdir layout logs empty sample warning before failure", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  sample_dir <- file.path(tmp_dir, "SampleB")
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  expect_error(
    check_fastq_files(tmp_dir, "subdir", log_file),
    "attempt to select less than one element in get1index"
  )

  log_entries <- readLines(log_file)
  expect_true(any(grepl("WARNING: SampleB has no FASTQ files", log_entries)))
})

test_that("subdir layout emits warnings for missing mates", {
  tmp_dir <- tempfile("fastq_input_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  sample_dir <- file.path(tmp_dir, "SampleC")
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)

  file.create(file.path(sample_dir, "SampleC_R1.fastq"))

  log_file <- create_temp_log()
  on.exit(unlink(log_file, force = TRUE), add = TRUE)

  default_log <- log_file
  assign(
    "write_log",
    function(message, log_file = NULL) {
      target <- if (is.null(log_file)) default_log else log_file
      original_write_log(message, target)
    },
    envir = globalenv()
  )
  on.exit(assign("write_log", original_write_log, envir = globalenv()), add = TRUE)

  warnings <- character()
  messages <- character()

  withCallingHandlers(
    check_fastq_files(tmp_dir, "subdir", log_file),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_true(any(grepl("R2 file is missing", warnings)))
  expect_true(any(grepl("I1 file is missing", warnings)))
  expect_true(any(grepl("should have 3 FASTQ files", messages)))

  log_entries <- readLines(log_file)
  expect_true(any(grepl("SampleC has FASTQ files", log_entries)))
  expect_true(any(grepl("WARNING: SampleC is missingR2 file", log_entries)))
  expect_true(any(grepl("WARNING: SampleC is missingI1 file", log_entries)))
})
