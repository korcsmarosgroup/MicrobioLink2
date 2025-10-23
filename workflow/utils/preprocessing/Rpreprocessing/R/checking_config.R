# Checking the configuration file: if it exists, it does nothing.
# If it doesn't, it throws an error.
#
# Args:
#
#   config_folder:
#     The system based absolute path of this script (because
#     the configuration file needs to be next to it).
#
# Error codes:
#
#   ERROR CODE 1: The configuration file does not exist.
#
# Returns:
#
#   String: the absolute path of the configuration file.


checking.configuration <- function (config_folder) {

  configuration.file.path <- paste0(config_folder, sep="/", "configuration.yaml")

  if (!file.exists(configuration.file.path)) {
    stop(paste0("ERROR CODE 1: The configuration file does not exist:", configuration.file.path))
  }

  return(configuration.file.path)

}
