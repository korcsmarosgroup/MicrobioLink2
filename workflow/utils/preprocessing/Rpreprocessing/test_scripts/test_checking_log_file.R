library(testthat)

source("C:/Users/bogla/PycharmProjects/MicrobioLink2/workflow/utils/preprocessing/Rpreprocessing/R/checking_log_file.R")

test_that(
  "checking_log_file creates a log file",
  {
    tmp_dir <- tempdir()
    config_folder <- file.path(tmp_dir, "config")
    dir.create(config_folder, showWarnings = FALSE)
    log_file <- file.path(tmp_dir, "preprocessing.log")
    result <- checking_log_file(config_folder)
    expect_equal(result, log_file)
    expect_true(file.exist(result))
  }
)
