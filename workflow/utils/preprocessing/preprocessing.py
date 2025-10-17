# This module/script will take either FASTQ files from a scRNA-seq analysis
# or it will take a raw count matrix from a scRNA-seq analysis


from time import strftime
import scanpy as sc
import subprocess
import sys
import os
import yaml
import glob

def _log(message, log_file):
    """
    Log a message with a timestamp.
    
    Args:
        message:
            The actual commit message before push
        log_file:
            The path of the log file
    """

    with open(log_file, "a") as log:
        log.write(f"[{strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")


def CheckingLogFile(config_folder):
    """
    Checking the log file and if it exists if will delete it.

    Args:
        config_folder:
            The system-based absolute path of this script (because the 
            configuration file needs to be next to it).
    """

    log_file = f"{config_folder}/preprocessing.log"

    if os.path.isfile(log_file):
        os.remove(log_file)


def CheckingConfiguration(config_folder):
    """
    Checking the configuration file: if it exists, then do nothing. If it doesn't
    then throws an error
    if not.

    Args:
        config_folder:
            The system-based absolute path of this script (because the 
            configuration file needs to be next to it).

    Error codes:
        ERROR CODE 1: The configuration file does not exist.

    Returns:
        String: the absolute path of the configuration file
    """

    configuration_file_path = f"{config_folder}/configuration.yaml"

    if not os.path.isfile(configuration_file_path):
        sys.stderr.write(f"ERROR MESSAGE: The configuration file does not "
                         f"exists: {configuration_file_path}")
        sys.exit(1)

    return configuration_file_path

def LoadingConfiguration(configuration_file_path, log_file):
    """
    Loading the content from the configuration file.

    Args:
        configuration_file_path:
            The absolute path of the configuration file.
        log_file:
            The path of the log file.

    Error codes:
        ERROR CODE 2: Unable to read configuration file.
        ERROR CODE 3: Configuration file is empty.

    Returns:
        Yaml object: contains the configuration file's content
    """

    try:
        with open(configuration_file_path, "r") as file:
            configuration = yaml.safe_load(file)

    except Exception as e:
        sys.stderr.write(f"ERROR MESSAGE: Unable to read configuration file: {e}\n")
        _log(f"Error reading configuration file: {e}", log_file)
        sys.exit(2)

    if configuration is None:
        sys.stderr.write("ERROR MESSAGE: Configuration file is empty.\n")
        _log("Configuration file is empty.", log_file)
        sys.exit(3)

    return configuration

def CheckFASTQFiles(input_folder, Fastq_file_format, log_file):
    """
    Checking the FASTQ input folder: verifies that all expected sequencing files are present
    and correctly structured according to the specified Fastq_file_format ('merged' or 'subdir').

    If any expected file (R1, R2, or I1) is missing, or the folder structure does not match 
    the specified layout, an error is raised, and execution stops.

    Args:
        input_folder:
            The system-based absolute path to the input directory containing the FASTQ files
            or per-sample subdirectories.
        Fastq_file_format:
            The organization style of the FASTQ data to validate. Must be either:
                - 'merged': all FASTQ files are in a single directory.
                - 'subdir': each sample has its own subdirectory containing FASTQ files.
        log_file:
            The path of the log file.

    Error codes:
        ERROR CODE 4: Invalid Fastq_file_format
        ERROR CODE 5: No FASTQ files found (merged layout)
        ERROR CODE 6: Subdirectory missing R1/R2/I1 (subdir layout)
        ERROR CODE 7: Subdirectory has wrong number of files (subdir layout)

    Returns:
        None
            The function does not return any value. If all checks pass, execution continues normally.
    """
    # ------------------------------------------------------------
    # Validate Fastq_file_format
    # ------------------------------------------------------------
    if Fastq_file_format not in ("merged", "subdir"):
        sys.stderr.write(f"ERROR: Invalid Fastq_file_format '{Fastq_file_format}'. Must be 'merged' or 'subdir'.\n")
        sys.exit(4)

    # ------------------------------------------------------------
    # Case 1: Flat layout — all FASTQs in one directory
    # ------------------------------------------------------------
    if Fastq_file_format == "merged":
        fastq_files = glob.glob(os.path.join(input_folder, "*.fastq")) + \
                      glob.glob(os.path.join(input_folder, "*.fastq.gz"))

        if len(fastq_files) == 0:
            sys.stderr.write(f"ERROR: No FASTQ files found in {input_folder}\n")
            sys.exit(5)

        _log("Using flat FASTQ layout", log_file)

        # Group files by sample prefix
        samples = {}
        for f in fastq_files:
            base = os.path.basename(f)
            if "_R1" in base:
                sample = base.split("_R1")[0]
            elif "_R2" in base:
                sample = base.split("_R2")[0]
            elif "_I1" in base:
                sample = base.split("_I1")[0]
            else:
                sys.stderr.write(f"WARNING: Skipping unrecognized FASTQ file name: {base}\n")
                _log(f"WARNING: Skipping unrecognized FASTQ file name: {base}")
                continue

            if sample not in samples:
                samples[sample] = []
            samples[sample].append(f)

        # Validate each sample
        for sample, files in samples.items():
            found = {"R1": False, "R2": False, "I1": False}
            for f in files:
                fname = os.path.basename(f)
                for tag in found:
                    if tag in fname:
                        found[tag] = True

            missing = [k for k, v in found.items() if not v]
            if missing:
                _log(
                    f"WARNING: {sample} is missing {', '.join(missing)} file(s).", log_file
                )
                continue

            if len(files) != 3:
                _log(
                    f"WARNING: {sample} should have 3 FASTQ files, but {len(files)} found.", log_file
                )
                continue

    # ------------------------------------------------------------
    # Case 2: Subdirectory layout — one folder per sample
    # ------------------------------------------------------------
    elif Fastq_file_format == "subdir":
        subdirs = [d for d in glob.glob(os.path.join(input_folder, "*")) if os.path.isdir(d)]

        if len(subdirs) == 0:
            sys.stderr.write(f"ERROR: No sample subdirectories found in {input_folder}\n")
            sys.exit(6)

        _log("Using per-sample subdirectory layout.\n", log_file)

        for folder in subdirs:
            sample_name = os.path.basename(folder)
            fastqs = glob.glob(os.path.join(folder, "*.fastq")) + \
                     glob.glob(os.path.join(folder, "*.fastq.gz"))

            if len(fastqs) == 0:
                sys.stderr.write(f"ERROR: No FASTQ files found in {folder}")
                sys.exit(7)

            found = {"R1": False, "R2": False, "I1": False}
            for f in fastqs:
                fname = os.path.basename(f)
                for tag in found:
                    if tag in fname:
                        found[tag] = True

            missing = [k for k, v in found.items() if not v]
            if missing:
                _log(
                    f"WARNING: {sample_name} is missing the {', '.join(missing)} file(s).", log_file
                )
                continue

            if len(fastqs) != 3:
                _log(
                    f"WARNING: {sample_name} should have exactly 3 FASTQ files, but {len(fastqs)} found.", log_file
                )
                continue

def PrepareSTARGenome(genome_dir, fasta_file, gtf_file, read_length, log_file):
    """
    Preparing the STAR genome index directory.

    Checks whether a valid STAR genome index already exists in the specified
    genome directory. If all required index files are present, the function
    logs this and skips index generation. If not, it generates a new STAR
    genome index using the provided FASTA and GTF annotation files.

    Args:
        genome_dir:
            The system-based absolute path of the directory where the STAR
            genome index will be stored.
        fasta_file:
            The absolute path to the reference genome FASTA file.
        gtf_file:
            The absolute path to the GTF annotation file.
        read_length:
            Integer, the sequencing read length. STAR will use (read_length - 1)
            as the sjdbOverhang value.
        log_file:
            The path of the log file to record events and errors.

    Error codes:
        ERROR CODE 12: STAR genome generation failed.

    Behavior:
        - If genome_dir exists and contains valid STAR index files
          (Genome, SA, and SAindex), the generation step is skipped.
        - If the index is missing, a new one is created with STAR.

    Returns:
        None
            The function logs progress and errors but returns no value.
    """

    # ------------------------------------------------------------
    # Check for existing valid STAR index
    # ------------------------------------------------------------
    if os.path.isdir(genome_dir):
        existing_files = os.listdir(genome_dir)
        required_files = {"Genome", "SA", "SAindex"}
        if required_files.issubset(set(existing_files)):
            _log(f"STAR genome index already exists at {genome_dir}. Skipping generation.", log_file)
            return

    # ------------------------------------------------------------
    # Create genome directory if missing
    # ------------------------------------------------------------
    _log(f"STAR genome index not found. Generating new index at {genome_dir}", log_file)
    os.makedirs(genome_dir, exist_ok=True)

    # ------------------------------------------------------------
    # Build STAR command
    # ------------------------------------------------------------
    cmd = [
        "STAR",
        "--runThreadN", "8",
        "--runMode", "genomeGenerate",
        "--genomeDir", genome_dir,
        "--genomeFastaFiles", fasta_file,
        "--sjdbGTFfile", gtf_file,
        "--sjdbOverhang", str(read_length - 1)
    ]

    # ------------------------------------------------------------
    # Run STAR to generate the genome index
    # ------------------------------------------------------------
    try:
        _log(f"Running STAR command: {' '.join(cmd)}", log_file)
        subprocess.run(cmd, check=True)
        _log(f"STAR genome index successfully generated at {genome_dir}", log_file)

    except subprocess.CalledProcessError as e:
        error_msg = f"ERROR: STAR genome generation failed: {e}"
        _log(error_msg, log_file)
        sys.stderr.write(f"ERROR MESSAGE: STAR genome generation failed: {e}\n")
        sys.exit(12)

def main():
    """
    Main function for the script.
    """

    config_folder = os.path.dirname(os.path.abspath(__file__))
    logfile = CheckingLogFile(config_folder)
    logfile = f"{config_folder}/preprocessing.log"  

    _log("Log file created", logfile)

    configuration_file_path = CheckingConfiguration(config_folder)

    _log("Configuration file is good, ready to start", logfile)
    
    configuration = LoadingConfiguration(configuration_file_path, logfile)
    _log("Configuration file successfully loaded.", logfile)
    
    input_dir = configuration.get("input_dir")
    Fastq_file_format = configuration.get("Fastq_file_format")  # must be "merged" or "subdir"

    if not input_dir or not os.path.isdir(input_dir):
        sys.stderr.write("ERROR: Invalid or missing 'input_dir' in configuration file.\n")
        _log("ERROR: Invalid or missing 'input_dir' in configuration file.", logfile)
        sys.exit(6)

    if not Fastq_file_format:
        sys.stderr.write("ERROR: Missing 'Fastq_file_format' in configuration file.\n")
        _log("ERROR: Missing 'Fastq_file_format' in configuration file.", logfile)
        sys.exit(6)

    _log(f"Checking FASTQ files in: {input_dir} (mode: {Fastq_file_format})", logfile)
    CheckFASTQFiles(input_dir, Fastq_file_format)
    _log("FASTQ file structure successfully validated.", logfile)
    
    genome_dir = configuration.get("genome_dir")
    fasta_file = configuration.get("fasta_file")
    gtf_file = configuration.get("gtf_file")
    read_length = configuration.get("read_length")
    
    # Validate config values before running
    if not all([genome_dir, fasta_file, gtf_file, read_length]):
        sys.stderr.write("ERROR: Missing required genome preparation parameters in configuration file.\n")
        _log("ERROR: Missing required genome preparation parameters in configuration file.", logfile)
        sys.exit(6)

    # Prepare STAR genome index
    PrepareSTARGenome(genome_dir, fasta_file, gtf_file, read_length, logfile)
        





if __name__ == '__main__':
    main()
