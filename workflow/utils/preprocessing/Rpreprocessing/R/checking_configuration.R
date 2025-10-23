#' checking_configuration
#'
#' Checking the configuration file: if it exists, it does nothing.
#' If it doesn't, it throws an error.
#'
#' @param config_foler The system based absolute path of this script (because
#' the configuration file needs to be next to it).
#'
#' @details
#' Error codes:
#'
#'   ERROR CODE 1: The configuration file does not exist.
#'
#' @return String: the absolute path of the configuration file.
#'
#' @export

checking_configuration <- function (config_folder) {

  configuration_file_path <- paste0(dirname(config_folder), sep="/", "configuration.yaml")

  if (!file.exists(configuration_file_path)) {
    stop(paste0("ERROR CODE 1: The configuration file does not exist: ", sep="\n", configuration_file_path))
  }

  return(configuration_file_path)

}
