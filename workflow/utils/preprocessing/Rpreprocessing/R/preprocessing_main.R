usethis::use_package("this.path", type="Imports")

#' preprocessing_main
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

  config_folder <- this.path::this.dir()

  log_file <- checking_log_file(config_folder)

  print(log_file)

  write_log("Log file created", log_file) # TO-DO: why message not displayed in log_file -> bad log_file path

  configuration_file_path <- checking_configuration(config_folder)

  write_log("Configuration file is good, ready to start", log_file)

  configuration <- loading_config_file(configuration_file_path, log_file)

  write_log("Configuration file is successfully loaded", log_file)

}
