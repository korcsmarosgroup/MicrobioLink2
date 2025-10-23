usethis::use_package("yaml", type="Imports")

#' loading_config_file
#'
#' Loading the content from the configuration file.
#'
#' @param configuration_file_path The absolute path of the configuration file.
#'
#' @param log_file The path of the log file.
#'
#' Error codes:
#'
#'   ERROR CODE 2: Unable to read configuration file.
#'   ERROR CODE 3: Configuration file is empty.
#'
#' @return Yaml object: contains the configuration file's content.
#'
#' @export

loading.config.file <- function(configuration.file.path, log.file) {

  if (file.size(configuration.file.path) == 0) {

    msg_1 <- "Configuration file is empty"

    cat(msg_1, log.file, append=TRUE)

    stop(paste0("The configuration file is empty"))

  }

  configuration <- tryCatch(

    {

      file <- yaml::read_yaml(configuration.file.path, stringsAsFactors=FALSE)
      return(file)

    },

    error=function(e) {

      msg_2 <- paste0("Error reading configuration file", configuration.file.path)

      cat(msg_2, file=log.file, append=TRUE)

      stop(paste0("Unable to read configuration file", configuration.file.path))

    }

  )

  return(configuration)

}
