#' add_log_message
#'
#' Log a message with a time stamp.
#'
#' @param message The actual commit message before push
#'
#' @param log_file The path of the log file
#' @keywords internal

write_log <- function (message, log_file) {

  time_stamp <- format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")

  write((paste0(time_stamp, sep=" ", message)), file = log_file, append = TRUE)

}

