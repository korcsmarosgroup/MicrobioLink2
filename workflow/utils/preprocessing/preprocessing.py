# This module/script will take either FASTQ files from a scRNA-seq analysis
# or it will take a raw count matrix from a scRNA-seq analysis


from time import strftime
import subprocess
import scanpy as sc
import sys
import os
import yaml
import glob

def _log(message, log_file):
    """Log a message with a timestamp."""

    with open(log_file, "a") as log:
        log.write(f"[{strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")


def CheckingLogFile(config_folder):
    """Checking the log file: create new one at every start of the module."""

    log_file = f"{config_folder}/preprocessing.log"

    if os.path.isfile(log_file):
        os.remove(log_file)

    return log_file


def CheckingConfiguration(config_folder):
    """
    Checking the configuration file: if it is existing and throws an error
    if not.
    """

    configuration_file_path = f"{config_folder}/configuration.yaml"

    if not os.path.isfile(configuration_file_path):
        sys.stderr.write(f"ERROR MESSAGE: The configuration file does not "
                         f"exists: {configuration_file_path}")
        sys.exit(1)
        
    return configuration_file_path

def LoadingConfiguration(configuration_file_path, log_file):

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

def check_fastq_files(input_folder):
    """
    Checks that each sample in the input folder has exactly 3 FASTQ files.
    Assumes files are named like:
        sample1_R1.fastq.gz, sample1_R2.fastq.gz, sample1_index.fastq.gz
    """
    
    fastq_files = glob.glob(os.path.join(input_folder, "*.fastq")) + \
              glob.glob(os.path.join(input_folder, "*.fastq.gz"))

    if not fastq_files:
        sys.stderr.write(f"ERROR: No FASTQ files found in {input_folder}\n")
        sys.exit(4)

    # Group by sample
    samples = {}
    for f in fastq_files:
        base = os.path.basename(f)
        sample_name = base.split("_")[0]

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
    # changed to a universial for all operating systems
    config_folder = os.path.dirname(os.path.abspath(__file__))
     #config_folder = "/".join(os.path.abspath(__file__).split("/")[:-1])


    logfile = CheckingLogFile(config_folder)

    _log("Log file created", logfile)

    configuration_file_path = CheckingConfiguration(config_folder)

    _log("Configuration file is good, ready to start", logfile)
    
    configuration = LoadingConfiguration(configuration_file_path, logfile)
    _log("Configuration file successfully loaded.", logfile)
    

    
        





if __name__ == '__main__':
    main()
