#' write_log_2
#' Generate a timestamped log file path inside the provided configuration folder.
#' @param config_folder str. Path to the folder containing the analysis configuration.
#' @return str Full file path to a new log file, named as analysis_<timestamp>.log.
#'
#' @export

get_log_file <- function(config_folder) {
  ts <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  log_file_2 <- file.path(config_folder, "analysis_", ts, ".log")
  return(log_file_2)
}
