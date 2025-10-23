# Log a message with a timestamp.
#
# Args:
#
#   message:
#     The actual commit message before push
#
#   log_file:
#     The path of the log file

log <- function (message, log.file) {

  write((paste0(format(Sys.time(), format = "%Y-%m-%d %H:%M:%S"), sep=" ", message)), file = log.file)

}

