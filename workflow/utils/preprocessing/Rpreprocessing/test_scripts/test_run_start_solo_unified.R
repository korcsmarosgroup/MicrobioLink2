pkg_wd <- getwd()
script_path <- file.path(pkg_wd, "R", "run_star_solo_unified.R")
if (!file.exists(script_path)) {
  pkg_wd <- normalizePath(file.path(pkg_wd, "..", "..", "..", "..", ".."), mustWork = TRUE)
  script_path <- file.path(pkg_wd, "R", "run_star_solo_unified.R")
}

source(file.path(pkg_wd, "R", "add_log_message.R"))
source(script_path)

create_temp_log <- function() {
  log_file <- tempfile("preprocessing_", fileext = ".log")
  file.create(log_file)
  log_file
}

define_configuration <- function(mode, platform = "microwell", params = list()) {
  temp_root <- tempdir()
  list(
    platform = platform,
    input_dir = "inputs",
    genome_dir = file.path(temp_root, "genome_index"),
    threads = 2,
    Fastq_file_format = "merged",
    STAR_outdir = file.path(temp_root, "star"),
    STAR_out = "output",
    STAR_params = params,
    STARsolo_outdir = file.path(temp_root, "starsolo"),
    STARsolo_out = "solo_output",
    STARsolo_params = params,
    mode = mode
  )
}

setup_globals <- function(configuration, samples_list, tag_values) {
  assign("file", configuration, envir = globalenv())
  assign("samples", samples_list, envir = globalenv())
  assign("tags", tag_values, envir = globalenv())
  assign("mode", configuration$mode, envir = globalenv())
}

build_output_dir <- function(base_dir, input_dir, suffix) {
  output_dir <- paste0(base_dir, input_dir, suffix, sep = "/")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir
}

test_that("run_STAR_unified builds STAR command for microwell data", {
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

  setup_globals(configuration, samples_list, tag_values)

  calls <- list()

  mock_run <- function(command, args, error_on_status) {
    calls <<- append(calls, list(list(command = command, args = args)))
    mocked_result
  }

  stub(run_STAR_solo_unified, "processx::run", mock_run)

  result <- run_STAR_solo_unified(configuration, log_file, solo = FALSE)

  expect_equal(length(calls), 1)
  expect_identical(calls[[1]]$command, "STAR")
  expect_true(any(grepl("--runThreadN", calls[[1]]$args)))
  expect_true(any(grepl(configuration$genome_dir, calls[[1]]$args)))
  expect_true(any(grepl("Sample1_R1.fastq", calls[[1]]$args)))
  expect_true(any(grepl("--outSAMtype", calls[[1]]$args)))
  expect_equal(result, mocked_result)
  expect_true(dir.exists(file.path(output_dir, "Sample1")))
})

test_that("run_STAR_unified logs and uses droplet settings when solo is requested", {
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

  calls <- list()
  mocked_result <- list(status = 0, stdout = "", stderr = "")

  result <- with_mocked_bindings(
    run_STAR_solo_unified(configuration, log_file, solo = TRUE),
    `processx::run` = function(command, args, error_on_status) {
      calls <<- append(calls, list(list(command = command, args = args, error_on_status = error_on_status)))
      mocked_result
    }
  )

  log_entries <- readLines(log_file)
  expect_true(any(grepl("Skipping STAR", log_entries)))
  expect_equal(length(calls), 1)
  expect_true(any(grepl("soloType", calls[[1]]$args)))
  expect_equal(result, mocked_result)
})
