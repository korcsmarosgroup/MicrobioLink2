source("./R/checking_configuration.R")

test_that("checking if config file path is returned by checking_configuration()", {
  # create temporary directory
  tmp_dir <- tempdir()
  config_folder <- file.path(tmp_dir, "config")
  dir.create(config_folder, showWarnings = FALSE)

  # create temporary configuration.yaml file
  configuration_file_path <- file.path(config_folder, "configuration.yaml")
  withr::local_file(configuration_file_path)
  file.create(configuration_file_path)

  # run function
  result <- checking_configuration(config_folder)

  # configuration file is where it is expected
  expect_equal(result, configuration_file_path)

  # configuration file exists
  expect_true(file.exists(configuration_file_path))
})

test_that("checking if error is raised if file doesn't exist", {
  # create temp dir
  tmp_dir <- tempdir()
  config_folder <- file.path(tmp_dir, "config")
  dir.create(config_folder, showWarnings = FALSE)

  # error is raised if config file doesn't exist
  expect_error(checking_configuration(config_folder))
})
