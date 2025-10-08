# This module/script will take either FASTQ files from a scRNA-seq analysis
# or it will take a raw count matrix from a scRNA-seq analysis


from time import strftime
import subprocess
import scanpy
import sys
import os


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


def main():
    """
    Main function for the script.
    """

    config_folder = "/".join(os.path.abspath(__file__).split("/")[:-1])

    logfile = CheckingLogFile(config_folder)

    _log("Log file created", logfile)

    CheckingConfiguration(config_folder)

    _log("Configuration file is good, ready to start", logfile)


if __name__ == '__main__':
    main()
