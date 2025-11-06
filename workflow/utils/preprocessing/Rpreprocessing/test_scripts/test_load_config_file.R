source("C:/Users/bogla/PycharmProjects/MicrobioLink2/workflow/utils/preprocessing/Rpreprocessing/R/loading_config_file.R")

test_that("check if error is raised if config file is empty", {
  # create temp dir
  tmp_dir <- tempdir()
  config_folder <- file.path(tmp_dir, "config")
  dir.create(config_folder, showWarnings = FALSE)

  # create temporary configuration.yaml file
  configuration_file_path <- file.path(config_folder, "configuration.yaml")
  withr::local_file(configuration_file_path)
  file.create(configuration_file_path)

  # create temporary log_file
  log_file <- tempfile("preprocessing", tmpdir = config_folder, ".log")
  local_file(log_file)
  file.create(log_file)

  # run function
  expect_error(loading_config_file(config_folder, log_file), "Configuration file is empty")
})
