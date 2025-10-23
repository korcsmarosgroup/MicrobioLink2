#' checking_configuration
#'
#' Checking the configuration file: if it exists, it does nothing.
#' If it doesn't, it throws an error.
#'
#' @param config_foler The system based absolute path of this script (because
#' the configuration file needs to be next to it).
#'
#' Error codes:
#'
#'   ERROR CODE 1: The configuration file does not exist.
#'
#' @return String: the absolute path of the configuration file.
#'
#' @export

checking.configuration <- function (config_folder) {

  configuration.file.path <- paste0(config_folder, sep="/", "configuration.yaml")

  if (!file.exists(configuration.file.path)) {
    stop(paste0("ERROR CODE 1: The configuration file does not exist: ", sep="\n", configuration.file.path))
  }

  return(configuration.file.path)

}
