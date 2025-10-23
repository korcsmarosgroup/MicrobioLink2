# Log a message with a timestamp.
#
# Args:
#
#   message:
#     The actual commit message before push
#
#   log_file:
#     The path of the log file

log <- function (message, log_file) {

  write(paste0(Sys.time(), sep=" ", message), file = log_file)

}

