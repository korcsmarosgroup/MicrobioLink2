usethis::use_package("processx", type = "Imports")

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

prepare_STAR_genome <- function(genome_dir, fasta_file, genome_index_params, splice_junction_params, log_file) {
  if (dir.exists(genome_dir)) {
    existing_files <- unlist(list.files(genome_dir))
    required_files <- c("Genome", "SA", "SAIndex")

    if ((length(required_files %in% existing_files)) == 3) {
      write_log(paste0("STAR genome index already exists at ", genome_dir, ". Skipping generation."), log_file)
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
  # Build STAR genomeGenerate command
  # ------------------------------------------------------------

  sjdbGTFfile <- splice_junction_params$sjdbGTFfile
  sjdbOverhang <- as.character(splice_junction_params$sjdbOverhang)
  genomeChrBinNbits <- as.character(genome_index_params$genomeChrBinNbits)
  genomeSAindexNbases <- as.character(genome_index_params$genomeSAindexNbases)
  genomeSAsparseD <- as.character(genome_index_params$genomeSAsparseD)
  genomeSuffixLengthMax <- as.character(genome_index_params$genomeSuffixLengthMax)

  cmd_base <- "STAR"
  cmd_args <- list(
    paste("--runThreadN", "8"),
    paste("--runMode", "genomeGenerate"),
    paste("--genomeDir", shQuote(genome_dir)),
    paste("--genomeFastaFiles", shQuote(fasta_file)),
    paste("--sjdbGTFfile", shQuote(sjdbGTFfile)),
    paste("--sjdbOverhang", shQuote(sjdbOverhang)),
    paste("--genomeChrBinNbits", shQuote(genomeChrBinNbits)),
    paste("--genomeSAindexNbases", shQuote(genomeSAindexNbases)),
    paste("--genomeSAsparseD", shQuote(genomeSAsparseD)),
    paste("--genomeSuffixLengthMax", shQuote(genomeSuffixLengthMax))
  )

  # Optional genome transform
  genome_transform_type <- genome_index_params$genomeTransformType
  genome_transform_VCF <- genome_index_params$genomeTransformVCF
  if (!is.null(genome_transform_type)) {
    cmd_transform <- list(
      paste("--genomeTransformType", shQuote(genome_transform_type))
    )
    cmd_args <- c(cmd_args, cmd_transform)
    if (file.exists(genome_transform_VCF)) {
      cmd_VCF <- list(
        paste("--genomeTransformVCF", shQuote(genome_transform_VCF))
      )
      cmd_args <- c(cmd_args, cmd_VCF)
    }
  }

  # Optional splice junction file
  sjdb_file_chr_start_end <- splice_junction_params$sjdbFileChrStartEnd
  if (file.exists(sjdb_file_chr_start_end)) {
    cmd_sjdb_file <- paste("--sjdbFileChrStartEnd", shQuote(sjdb_file_chr_start_end))
    cmd_args <- c(cmd_args, cmd_sjdb_file)
  }

  # Splice junction tags
  sjdbGTFchrPrefix <- splice_junction_params$sjdbGTFchrPrefix
  sjdbGTFfeatureExon <- splice_junction_params$sjdbGTFfeatureExon
  sjdbGTFtagExonParentTranscript <- splice_junction_params$sjdbGTFtagExonParentTranscript
  sjdbGTFtagExonParentGene <- splice_junction_params$sjdbGTFtagExonParentGene
  sjdbGTFtagExonParentGeneName <- splice_junction_params$sjdbGTFtagExonParentGeneName
  sjdbGTFtagExonParentGeneType <- splice_junction_params$sjdbGTFtagExonParentGeneType
  sjdbScore <- as.character(splice_junction_params$sjdbScore)

  cmd_splice_junction_tags <- list(
    paste("--sjdbGTFchrPrefix", shQuote(sjdbGTFchrPrefix)),
    paste("--sjdbGTFfeatureExon", shQuote(sjdbGTFfeatureExon)),
    paste("--sjdbGTFtagExonParentTranscript", shQuote(sjdbGTFtagExonParentTranscript)),
    paste("--sjdbGTFtagExonParentGene", shQuote(sjdbGTFtagExonParentGene)),
    paste("--sjdbGTFtagExonParentGeneName", shQuote(sjdbGTFtagExonParentGeneName)),
    paste("--sjdbGTFtagExonParentGeneType", shQuote(sjdbGTFtagExonParentGeneType)),
    paste("--sjdbScore", shQuote(sjdbScore)),
    paste("--sjdbInsertSave", shQuote(splice_junction_params$sjdbInsertSave))
  )

  cmd <- c(cmd_base, cmd_args, cmd_splice_junction_tags)

  pipeline <- paste(cmd, collapse = " ")

  # ------------------------------------------------------------
  # Run STAR to generate the genome index
  # ------------------------------------------------------------

  STAR_result <- tryCatch(
    {
      write_log(paste0("Running STAR command: ", pipeline), log_file)
      result <- processx::run(cmd_base, c(unlist(cmd_args), unlist(cmd_splice_junction_tags)), error_on_status = TRUE)
      write_log(paste0("STAR genome index generated at ", genome_dir), log_file)
      return(result)
    },
    error = function(e) {
      cat("ERROR CODE 11: STAR genome generation failed", conditionMessage(e))
      return(NULL)
    },
    warning = function(w) {
      cat("WARNING:", conditionMessage(w))
      return(NULL)
    }
  )
  return(STAR_result)
}
