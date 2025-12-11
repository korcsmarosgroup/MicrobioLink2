source("./R/preparing_start_genome.R")
source("./R/add_log_message.R")

#' Helper: create minimal genome and splice junction parameter lists for testing
minimal_genome_index_params <- function() {
  list(
    genomeChrBinNbits = 18,
    genomeSAindexNbases = 12,
    genomeSAsparseD = 1,
    genomeSuffixLengthMax = 0,
    genomeTransformType = NULL,
    genomeTransformVCF = tempfile("vcf")
  )
}

minimal_splice_junction_params <- function(tmp_dir) {
  list(
    sjdbGTFfile = file.path(tmp_dir, "annotation.gtf"),
    sjdbOverhang = 100,
    sjdbFileChrStartEnd = file.path(tmp_dir, "sjdb.tsv"),
    sjdbGTFchrPrefix = "",
    sjdbGTFfeatureExon = "exon",
    sjdbGTFtagExonParentTranscript = "transcript_id",
    sjdbGTFtagExonParentGene = "gene_id",
    sjdbGTFtagExonParentGeneName = "gene_name",
    sjdbGTFtagExonParentGeneType = "gene_type",
    sjdbScore = 2,
    sjdbInsertSave = "Basic"
  )
}

test_that("prepare_STAR_genome exits early when index exists", {
  tmp_dir <- tempdir()
  genome_dir <- file.path(tmp_dir, "existing_genome")
  dir.create(genome_dir, showWarnings = FALSE)
  file.create(file.path(genome_dir, "Genome"))
  file.create(file.path(genome_dir, "SA"))
  file.create(file.path(genome_dir, "SAIndex"))

  log_file <- tempfile("preprocessing", tmpdir = tmp_dir, fileext = ".log")
  withr::local_file(log_file)

  result <- prepare_STAR_genome(
    genome_dir,
    fasta_file = file.path(tmp_dir, "genome.fa"),
    genome_index_params = minimal_genome_index_params(),
    splice_junction_params = minimal_splice_junction_params(tmp_dir),
    log_file = log_file
  )

  expect_null(result)
  log_lines <- readLines(log_file)
  expect_true(any(grepl("already exists", log_lines)))
})

test_that("prepare_STAR_genome builds and runs STAR command when index is missing", {
  tmp_dir <- tempdir()
  genome_dir <- file.path(tmp_dir, "new_genome")
  if (dir.exists(genome_dir)) unlink(genome_dir, recursive = TRUE)

  log_file <- tempfile("preprocessing", tmpdir = tmp_dir, fileext = ".log")
  withr::local_file(log_file)

  captured_args <- list()
  mock_run <- function(command, args, error_on_status) {
    captured_args$command <<- command
    captured_args$args <<- args
    captured_args$error_on_status <<- error_on_status
    list(status = 0, stdout = "mocked")
  }

  # Save original
  original_run <- get("run", envir = asNamespace("processx"))

  # Patch run() in processx namespace
  assignInNamespace("run", mock_run, ns = "processx")

  # Ensure restoration even if test fails
  withr::defer(assignInNamespace("run", original_run, ns = "processx"))

  params_genome <- minimal_genome_index_params()
  params_splice <- minimal_splice_junction_params(tmp_dir)

  result <- prepare_STAR_genome(
    genome_dir,
    fasta_file = file.path(tmp_dir, "genome.fa"),
    genome_index_params = params_genome,
    splice_junction_params = params_splice,
    log_file = log_file
  )

  expect_true(dir.exists(genome_dir))
  expect_equal(result$status, 0)
  expect_equal(captured_args$command, "STAR")
  expect_true(any(grepl("--runMode", captured_args$args)))
  expect_true(any(grepl("--genomeDir", captured_args$args)))
  expect_true(any(grepl("--sjdbGTFfile", captured_args$args)))

  log_lines <- readLines(log_file)
  expect_true(any(grepl("Running STAR command", log_lines)))
  expect_true(any(grepl("genome index generated", log_lines)))
})
