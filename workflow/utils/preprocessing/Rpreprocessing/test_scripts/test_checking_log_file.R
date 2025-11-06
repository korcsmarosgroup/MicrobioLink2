library(testthat)

source("C:/Users/bogla/PycharmProjects/MicrobioLink2/workflow/utils/preprocessing/Rpreprocessing/R/checking_log_file.R")

test_that(
  "checking_log_file creates a log file",
  {
    # create a temporary config folder
    tmp_dir <- tempdir()
    config_folder <- file.path(tmp_dir)
    dir.create(config_folder, showWarnings = FALSE)

    # create file path where log_file is expected
    log_file <- file.path(config_folder, "preprocessing.log")

    # run function
    result <- checking_log_file(config_folder)

    # log_file is where it is expected
    expect_equal(result, log_file)

    # file exists
    expect_true(file.exists(result))

    # file is empty
    content <- readLines(result, warn = FALSE)
    expect_true(length(content) == 0)
  }
)
