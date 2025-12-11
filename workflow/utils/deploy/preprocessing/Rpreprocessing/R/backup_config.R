usethis::use_package("tools", type = "Imports")
#' backup_config
#'
#' Copy YAML files (*.yaml, *.yml) from config_folder into output_folder.
#' Adds a timestamp to each file name instead of creating a backup folder.
#' Example: config.yaml → config_2025-02-10_15-32-18.yaml
#'
#' @export

backup_config <- function(config_folder, output_folder) {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  dir.create(output_folder)

  config_files <- list.files(config_folder, pattern = "\\.(yml|yaml)$")

  for (filename in config_files) {
    src_file <- file.path(config_folder, filename)

    # Split filename into name + extension
    name <- tools::file_path_sans_ext(filename)
    ext <- tools::file_ext(filename)

    # New filename with timestamp
    new_filename <- paste0(name, "_", timestamp, ext)
    dst_folder <- file.path(output_folder, "/")

    file.copy(src_file, dst_folder)
  }

  return(output_folder)
}
