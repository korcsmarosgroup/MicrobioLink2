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

checking_log_file <- function (config_folder) {

  log_file <- paste0(dirname(config_folder), sep="/", "preprocessing.log")

  if (file.exists(log_file)) {
    file.remove(log_file)
  }

  file.create(log_file)

  return(log_file)

}
