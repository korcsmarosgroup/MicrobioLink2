usethis::use_package("this.path", type = "Imports")
#' main
#'
#' Main function for the script.
#'
#' @details
#' ERROR CODE 4: Invalid or missing 'input_dir' in configuration file
#' ERROR CODE 9: Missing 'Fastq_file_format' in configuration file.
#' ERROR CODE 10: Missing required genome preparation parameters in configuration file
#' ERROR CODE 12: Platform is '{platform}'. No processing available for this platform. Exiting
#'
#' @export
main <- function() {
  config_folder <- dirname(this.path::this.dir())

  log_file <- checking_log_file(config_folder)

  log_file <- file.path(config_folder, "preprocessing.log")

  write_log("Log file created", log_file)

  configuration_file_path <- checking_configuration(config_folder)

  write_log("Configuration file is good, ready to start", log_file)

  configuration <- loading_config_file(configuration_file_path, log_file)

  write_log("Configuration file is successfully loaded", log_file)

  input_dir <- configuration$input_dir
  Fastq_file_format <- configuration$Fastq_file_format # must be "merged" or "subdir"

  if (isTRUE((is.null(input_dir)) || (!dir.exists(input_dir)))) {
    stop("ERROR CODE 4: Invalid or missing 'input_dir' in configuration file.")
    write_log("ERROR CODE 4: Invalid or missing /'input_dir/' in configuration file.", log_file)
  }

  if (is.null(Fastq_file_format)) {
    stop("ERROR CODE 9: Missing 'Fastq_file_format' in configuration file.")
    write_log("ERROR CODE 9: Missing 'Fastq_file_format' in configuration file.", log_file)
  }

  write_log(paste0("Checking FASTQ files in: ", input_dir, " (mode: ", Fastq_file_format, " )"), log_file)
  check_fastq_files(input_dir, Fastq_file_format, log_file)
  write_log("FASTQ file structure successfully validated.", log_file)

  genome_dir <- configuration$genome_dir
  fasta_file <- configuration$fasta_file
  genome_index_params <- configuration$genome_index_params
  splice_junction_params <- configuration$splice_junction_params
  HVG_selection <- configuration$HVG_selection
  number_top_genes <- configuration$number_top_genes

  # Validate config values before running
  if (any(sapply(list(genome_dir, fasta_file, gtf_file, read_length), is.null))) {
    stop("ERROR CODE 10: Missing required genome preparation parameters in configuration file.")
    write_log("ERROR CODE 10: Missing required genome preparation parameters in configuration file.", log_file)
  }

  # Prepare STAR genome index
  prepare_STAR_genome(genome_dir, fasta_file, gtf_file, read_length, log_file)
  # check platform and run star
  platform <- tolower(configuration$platform)

  if (isTRUE(platform == "droplet")) {
    run_STAR_solo_unified(configuration, log_file, solo = TRUE)
    QC_10x(output_folder, log_file)
  } else if (isTRUE(platform == "microwell")) {
    run_STAR_unified(configuration, log_file, solo = FALSE)
  } else {
    msg <- paste0("Platform is ", platform, " No processing available for this platform. Exiting.")
    write_log(msg, log_file)
    stop(paste0("ERROR CODE 12: ", msg))
  }

  normalize_all_samples(output_dir, HVG_selection, number_top_genes, log_file)
}
