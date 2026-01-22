source("./R/preparing_start_genome.R")
source("./R/add_log_message.R")

test_that("run_STAR_unified builds STAR command for microwell data", {
  testthat::skip_if_not_installed("mockery")

  configuration <- define_configuration(
    mode = "STAR",
    platform = "microwell",
    params = list(outSAMtype = "BAM SortedByCoordinate")
  )

  samples_list <- list(Sample1 = c("Sample1_R1.fastq", "Sample1_R2.fastq"))
  tag_values <- c("R1", "R2")

  output_dir <- build_output_dir(configuration$STAR_outdir, configuration$input_dir, configuration$STAR_out)
  log_file <- create_temp_log()

  on.exit(unlink(c(log_file, output_dir), recursive = TRUE, force = TRUE), add = TRUE)

  # Function relies on globals: file, samples, tags, mode, AND (when solo=FALSE) output_dir + params
  setup_globals(configuration, samples_list, tag_values)
  assign("output_dir", output_dir, envir = globalenv())
  assign("params", configuration$STAR_params, envir = globalenv())

  # Clean up globals so tests don't leak into each other
  on.exit(
    rm(
      list = intersect(
        c("file", "samples", "tags", "mode", "output_dir", "params"),
        ls(envir = globalenv())
      ),
      envir = globalenv()
    ),
    add = TRUE
  )

  calls <- list()
  mocked_result <- list(status = 0, stdout = "", stderr = "")

  mock_run <- function(command, args, error_on_status) {
    calls <<- append(calls, list(list(command = command, args = args, error_on_status = error_on_status)))
    mocked_result
  }

  mockery::stub(run_STAR_solo_unified, "processx::run", mock_run)

  result <- run_STAR_solo_unified(configuration, log_file, solo = FALSE)

  expect_equal(length(calls), 1)
  expect_identical(calls[[1]]$command, "STAR")
  expect_true(any(grepl("--runThreadN", calls[[1]]$args)))
  expect_true(any(grepl(configuration$genome_dir, calls[[1]]$args, fixed = TRUE)))
  expect_true(any(grepl("Sample1_R1.fastq", calls[[1]]$args)))
  expect_true(any(grepl("Sample1_R2.fastq", calls[[1]]$args)))

  # Params are currently reflected in the *logged pipeline* (not in processx::run args)
  log_entries <- readLines(log_file, warn = FALSE)
  expect_true(any(grepl("--outSAMtype", log_entries)))

  expect_equal(result, mocked_result)
  expect_true(dir.exists(file.path(output_dir, "Sample1")))
})


test_that("run_STAR_unified logs and uses droplet settings when solo is requested", {
  testthat::skip_if_not_installed("mockery")

  configuration <- define_configuration(
    mode = "STARsolo",
    platform = "droplet",
    params = list(soloType = "CB_UMI_Simple")
  )

  samples_list <- list(CellA = c("CellA_R1.fastq", "CellA_R2.fastq"))
  tag_values <- c("R1", "R2")

  output_dir <- build_output_dir(configuration$STAR_outdir, configuration$input_dir, configuration$STAR_out)
  log_file <- create_temp_log()
  on.exit(unlink(c(log_file, output_dir), recursive = TRUE, force = TRUE), add = TRUE)

  setup_globals(configuration, samples_list, tag_values)

  on.exit(
    rm(
      list = intersect(
        c("file", "samples", "tags", "mode"),
        ls(envir = globalenv())
      ),
      envir = globalenv()
    ),
    add = TRUE
  )

  calls <- list()
  mocked_result <- list(status = 0, stdout = "", stderr = "")

  mock_run <- function(command, args, error_on_status) {
    calls <<- append(calls, list(list(command = command, args = args, error_on_status = error_on_status)))
    mocked_result
  }

  mockery::stub(run_STAR_solo_unified, "processx::run", mock_run)

  result <- run_STAR_solo_unified(configuration, log_file, solo = TRUE)

  log_entries <- readLines(log_file, warn = FALSE)
  expect_true(any(grepl("Skipping STAR", log_entries)))

  expect_equal(length(calls), 1)
  expect_identical(calls[[1]]$command, "STAR")

  # Params currently show up in the logged pipeline
  expect_true(any(grepl("--soloType", log_entries)))

  expect_equal(result, mocked_result)
})
