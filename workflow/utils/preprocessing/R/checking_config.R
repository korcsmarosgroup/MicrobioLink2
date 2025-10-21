"""


"""

checking.configuration <- function (config_folder) {

  configuration.file.path <- paste0(config_folder, sep="/", "configuration.yaml")

  if (!file.exists(configuration.file.path)) {
    stop(paste0("ERROR CODE 1: The configuration file does not exist:", configuration.file.path))
  }

  return(configuration.file.path)

}
