install.packages("yaml")

loading.config.file <- function(configuration.file.path, log.file) {

  if (file.size(configuration.file.path) == 0) {

    msg_1 <- "Configuration file is empty"

    cat(msg_1, log.file, append=TRUE)

    stop(paste0("The configuration file is empty"))

  }

  configuration <- tryCatch(

    {

      file <- read_yaml(configuration.file.path, stringsAsFactors=FALSE)
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
