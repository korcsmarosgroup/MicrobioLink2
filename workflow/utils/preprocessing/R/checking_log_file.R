install.packages("logr")

checking.log.file <- function (config_folder) {

  log.file <- paste0(config_folder, sep="/", "preprocessing.log")

  if (file.exists(log.file)) {
    file.remove(log.file)
  }

  return(log_file)

}
