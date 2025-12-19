#' add_log_message_2
#'
#' Log a message with a time stamp.
#'
#' @param message The actual commit message before push
#' @param log_file_2 The path of the downstream analysis log file
#'
#' @keywords internal

write_log_2 <- function(message, log_file_2) {
  time_stamp <- format(Sys.time(), format = "%Y-%m-%d %H:%M:%S")
  write((paste0(time_stamp, sep = " ", message)), file = log_file_2, append = TRUE)
}
