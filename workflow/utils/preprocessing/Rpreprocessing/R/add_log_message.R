#' add_log_message
#'
#' Log a message with a time stamp.
#'
#' @param message The actual commit message before push
#'
#' @param log.file The path of the log file
#' @keywords internal

log <- function (message, log.file) {

  write((paste0(format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"), sep=" ", message)), file = log.file)

}

