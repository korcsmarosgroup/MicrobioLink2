#' checking_log_file
#'
#' Checking the log file and if it exists, it will delete it.
#'
#' @param config_folder The system-based absolute path of this script (because the configuration file needs to be
#' next to it.)
#'
#' @return The full log file path
#'
#' @export

checking.log.file <- function (config_folder) {

  log.file <- paste0(config_folder, sep="/", "preprocessing.log")

  if (file.exists(log.file)) {
    file.remove(log.file)
  }

  return(log.file)

}
