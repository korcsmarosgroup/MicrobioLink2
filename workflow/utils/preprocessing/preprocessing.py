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

def CheckFASTQFiles(input_folder):
    """
    Checks that each sample in the input folder has exactly 3 FASTQ files.
    Assumes files are named like:
        sample1_R1.fastq.gz, sample1_R2.fastq.gz, sample1_index.fastq.gz

    #TODO: Check if we have more than 3 FASTQ file pairs, because of subsamples
    """
    
    fastq_files = glob.glob(os.path.join(input_folder, "*.fastq")) + \
              glob.glob(os.path.join(input_folder, "*.fastq.gz"))

    if len(fastq_files) == 0:
        sys.stderr.write(f"ERROR: No FASTQ files found in {input_folder}\n")
        sys.exit(4)

    # Group by sample
    samples = {}
    for f in fastq_files:
        base = os.path.basename(f)
        sample_R1_name = base.split("R1")[0]
        sample_R2_name = base.split("R2")[0]
        sample_I1_name = base.split("I1")[0]        

        if sample_name not in samples:
            samples[sample_name] = []

        samples[sample_name].append(f)

    # Check each sample has 3 files
    for sample, files in samples.items():
        if len(files) != 3:
            sys.stderr.write(f"ERROR: Sample '{sample}' should have 3 FASTQ files, but found {len(files)}.\n")
            sys.exit(5)


def main():
    """
    Main function for the script.
    """

    config_folder = os.path.dirname(os.path.abspath(__file__))
    logfile = CheckingLogFile(config_folder)

    _log("Log file created", logfile)

    configuration_file_path = CheckingConfiguration(config_folder)

    _log("Configuration file is good, ready to start", logfile)
    
    configuration = LoadingConfiguration(configuration_file_path, logfile)
    _log("Configuration file successfully loaded.", logfile)
    

    
        





if __name__ == '__main__':
    main()
