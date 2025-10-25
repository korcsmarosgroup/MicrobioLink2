prepare_STAR_genome <- function (genome_dir, fasta_file, gtf_file, read_length, log_file) {

  stop("function not implemented yet")


#' Preparing the STAR genome index directory.
#'
#' Checks whether a valid STAR genome index already exists in the specified
#' genome directory. If all required index files are present, the function
#' logs this and skips index generation. If not, it generates a new STAR
#' genome index using the provided FASTA and GTF annotation files.
#'
#' @param genome_dir The system-based absolute path of the directory where the STAR
#' genome index will be stored.
#'
#' @param fasta_file The absolute path to the reference genome FASTA file.
#' @param gtf_file The absolute path to the GTF annotation file.
#' @param read_length Integer, the sequencing read length. STAR will use (read_length - 1)
#' as the sjdbOverhang value.
#' @param log_file The path of the log file to record events and errors.
#'
#' @details
#' Error codes:
#'    ERROR CODE 8: STAR genome generation failed.
#'
#' Behavior:
#'  - If genome_dir exists and contains valid STAR index files (Genome, SA, and SAindex), the generation step is
#'  skipped.
#'  - If the index is missing, a new one is created with STAR.
#'
#'  @return None The function logs progress and errors but returns no value.

    # ------------------------------------------------------------
    # Check for existing valid STAR index
    # ------------------------------------------------------------

    if (dir.exists(genome_dir)) {

        existing_files <- unlist(list.files(genome_dir))
        required_files <- c("Genome", "SA", "SAIndex")

        if ((length(existing_files %in% required_files)) == 3) {
            write_log(paste0("STAR genome index already exists at ", genome_dir, " Skipping generation."), log_file)
            return(NULL)
       }

     }

    # ------------------------------------------------------------
    # Create genome directory if missing
    # ------------------------------------------------------------

    else {
        write_log(paste0("STAR genome index not found. Generating new index at ", genome_dir), log_file)
        dir.create(genome_dir, recursive = TRUE)

     }

    # ------------------------------------------------------------
    # Build STAR command
    # ------------------------------------------------------------

    cmd <- "STAR"
    args <- c("--runThreadN", "8",
               "--runMode", "genomeGenerate",
               "--genomeDir", genome_dir,
               "--genomeFastaFiles", fasta_file,
               "--sjdbGTFfile", gtf_file,
               "--sjdbOverhang", is.character(read_length - 1))

    # ------------------------------------------------------------
    # Run STAR to generate the genome index
    # ------------------------------------------------------------


  tryCatch(

      write_log("Running STAR command", log_file),
      system2(cmd, args = args),
      write_log(paste0("STAR genome index generated at ", genome_dir), log_file),

      error = function(e) cat("ERROR CODE 11: START genome generation failed")

    )

}
