checking_configuration <- function(config_folder) {
  config_path <- file.path(config_folder, "analysis_config.yaml")
  if (!(file.exists(config_path))) {
    stop(paste("ERROR CODE 1: Missing analysis_config.yaml at ", config_path))
  }
  return(config_path)
}
