source("C:/Users/bogla/PycharmProjects/MicrobioLink2/workflow/utils/preprocessing/Rpreprocessing/R/add_log_message.R")

test_that("checking if messages are added to log_file", {
  # create temp dir
  tmp_dir <- tempdir()
  config_folder <- file.path(tmp_dir, "config")
  dir.create(config_folder, showWarnings = FALSE)

  # create temp log_file
  log_file <- tempfile("preprocessing", tmpdir = config_folder, ".log")
  local_file(log_file)
  file.create(log_file)

  # check if content is added
  write_log("This is a test.", log_file)
  content <- readLines(log_file)
  expect_true(length(content) > 0)
})
