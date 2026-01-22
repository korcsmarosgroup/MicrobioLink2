source("./R/loading_config_file.R")

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
  expect_error(loading_config_file(config_folder, log_file), "ERROR CODE 3: Configuration file is empty")
})

test_that("check if valid file content is returned", {
  # create temp dir
  tmp_dir <- tempdir()
  config_folder <- file.path(tmp_dir, "config")
  dir.create(config_folder, showWarnings = FALSE)

  # create dummy configuration.yaml file
  configuration_file_path <- file.path(config_folder, "configuration.yaml")
  withr::local_file(configuration_file_path)
  file.create(configuration_file_path)
  data <- list(input_dir = "/input_dir", output_dir = "./results")
  yaml::write_yaml(data, configuration_file_path)

  # create temporary log_file
  log_file <- tempfile("preprocessing", tmpdir = config_folder, ".log")
  local_file(log_file)
  file.create(log_file)

  # run function
  result <- loading_config_file(configuration_file_path, log_file)

  # check config file structure
  expect_true(is.list(result))
  expect_true(result$input_dir == "/input_dir")
  expect_true(result$output_dir == "./results")
})

test_that("check if error is raised if config file cannot be parsed correctly", {
  # create a temp dir
  tmp_dir <- tempdir()
  config_folder <- file.path(tmp_dir, "config")
  dir.create(config_folder, showWarnings = FALSE)

  # create temporary unreadable "bad_configuration.yaml" file
  configuration_file_path <- file.path(config_folder, "bad_configuration.yaml")
  withr::local_file(configuration_file_path)
  file.create(configuration_file_path)
  writeLines("::: bad yaml :::", configuration_file_path)

  # create temporary log_file
  log_file <- tempfile("preprocessing", tmpdir = config_folder, ".log")
  withr::local_file(log_file)
  file.create(log_file)

  # run function
  error_msg <- paste0("ERROR CODE 2: Unable to read configuration file:", sep = " ", configuration_file_path)
  expect_error(loading_config_file(configuration_file_path, log_file), error_msg, fixed = TRUE)
})
