#' get_log_file
#'
#' Create a timestamped log filename.
#'
#' @param config_folder The system-based absolute path of this script (because the configuration file needs to be
#' next to it.)
#'
#' @return The full log file path
#'
#' @export

get_log_file <- function(config_folder) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file <- file.path(config_folder, paste0("preprocessing_", timestamp, ".log"))
  return(log_file)
}
