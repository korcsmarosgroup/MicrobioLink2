import pytest
from preprocessing import LoadingConfiguration
from preprocessing import CheckingConfiguration
from preprocessing import CheckFASTQFiles
from preprocessing import PrepareSTARGenome
import subprocess

from pathlib import Path
import os
import yaml


@pytest.fixture
def config_folder():
    return Path(__file__).resolve().parent



def test_CheckingLogFile(config_folder):

    test_log_file = f"{config_folder}/test_preprocessing.log"

    if os.path.isfile(test_log_file):
        os.remove(test_log_file)

    assert not os.path.exists(test_log_file)


def test_CheckingConfiguration(tmp_path, capsys):

    config_folder = tmp_path

    with pytest.raises(SystemExit) as error:
        CheckingConfiguration(str(config_folder))

    assert error.value.code == 1

    error_message = capsys.readouterr()
    assert f"ERROR MESSAGE: The configuration file does not exists:" in error_message.err

def test_LoadingConfiguration_valid_yaml(tmp_path):
    """Should load a valid YAML file and return its contents."""
    config_file = tmp_path / "config.yaml"
    log_file = tmp_path / "log.txt"

    data = {"input_dir": "data/", "output_dir": "results/"}
    config_file.write_text(yaml.dump(data))

    result = LoadingConfiguration(str(config_file), str(log_file))
    assert isinstance(result, dict)
    assert result["input_dir"] == "data/"
    assert result["output_dir"] == "results/"


def test_LoadingConfiguration_unreadable_file(tmp_path, capsys):
    """Should exit with code 2 if YAML file cannot be read."""
    config_file = tmp_path / "bad_config.yaml"
    log_file = tmp_path / "log.txt"

    # create a malformed file (invalid YAML syntax)
    config_file.write_text("::: not valid yaml :::")

    with pytest.raises(SystemExit) as e:
        LoadingConfiguration(str(config_file), str(log_file))

    assert e.value.code == 2

    out, err = capsys.readouterr()
    assert "unable to read configuration file" in (out + err).lower()


def test_LoadingConfiguration_empty_file(tmp_path, capsys):
    """Should exit with code 3 if YAML file is empty."""
    config_file = tmp_path / "empty_config.yaml"
    log_file = tmp_path / "log.txt"

    config_file.write_text("")  # empty file

    with pytest.raises(SystemExit) as e:
        LoadingConfiguration(str(config_file), str(log_file))

    assert e.value.code == 3

    out, err = capsys.readouterr()
    assert "configuration file is empty" in (out + err).lower()
    
    
def create_fastq(folder: Path, name: str):
    """Helper to create dummy FASTQ files."""
    f = folder / name
    f.write_text("dummy FASTQ content")
    return f


# ---------------------------------------------------------
# TESTS
# ---------------------------------------------------------

def test_CheckFASTQFiles_invalid_format(tmp_path, capsys):
    """Should exit with code 5 if format is invalid."""
    log_file = tmp_path / "test.log"

    with pytest.raises(SystemExit) as e:
        CheckFASTQFiles(str(tmp_path), "invalid_format", str(log_file))

    assert e.value.code == 5

    out, err = capsys.readouterr()
    assert "Invalid Fastq_file_format" in (out + err)


def test_CheckFASTQFiles_merged_no_fastqs(tmp_path, capsys):
    """Should exit with code 6 if merged layout has no FASTQ files."""
    log_file = tmp_path / "test.log"

    with pytest.raises(SystemExit) as e:
        CheckFASTQFiles(str(tmp_path), "merged", str(log_file))

    assert e.value.code == 6
    out, err = capsys.readouterr()
    assert "No FASTQ files found" in (out + err)


def test_CheckFASTQFiles_merged_valid_files(tmp_path):
    """Should pass with correct R1, R2, and I1 files in merged layout."""
    log_file = tmp_path / "test.log"

    create_fastq(tmp_path, "sample_R1.fastq")
    create_fastq(tmp_path, "sample_R2.fastq")
    create_fastq(tmp_path, "sample_I1.fastq")

    # Should not exit
    CheckFASTQFiles(str(tmp_path), "merged", str(log_file))

    # Log file should exist and contain expected message
    assert log_file.exists()
    assert "Using flat FASTQ layout" in log_file.read_text()


def test_CheckFASTQFiles_subdir_no_subdirs(tmp_path, capsys):
    """Should exit with code 7 if no sample subdirectories exist."""
    log_file = tmp_path / "test.log"

    with pytest.raises(SystemExit) as e:
        CheckFASTQFiles(str(tmp_path), "subdir", str(log_file))

    assert e.value.code == 7
    out, err = capsys.readouterr()
    assert "No sample subdirectories found" in (out + err)


# def test_CheckFASTQFiles_subdir_empty_folder(tmp_path, capsys):
#     """Should exit with code 8 if a sample folder has no FASTQ files."""
#     log_file = tmp_path / "test.log"
#     sample_folder = tmp_path / "sampleA"
#     sample_folder.mkdir()

#     with pytest.raises(SystemExit) as e:
#         CheckFASTQFiles(str(tmp_path), "subdir", str(log_file))

#     assert e.value.code == 8
#     out, err = capsys.readouterr()
#     assert "No FASTQ files found" in (out + err)


def test_CheckFASTQFiles_subdir_valid_files(tmp_path):
    """Should pass with correct files inside sample subdirectory."""
    log_file = tmp_path / "test.log"
    sample_folder = tmp_path / "sampleB"
    sample_folder.mkdir()

    create_fastq(sample_folder, "sampleB_R1.fastq.gz")
    create_fastq(sample_folder, "sampleB_R2.fastq.gz")
    create_fastq(sample_folder, "sampleB_I1.fastq.gz")

    # Should not exit
    CheckFASTQFiles(str(tmp_path), "subdir", str(log_file))

    # Check that log file has correct message
    assert log_file.exists()
    assert "sampleB has the necessary FASTQ files (R1, R2, I1)" in log_file.read_text()
    
def test_PrepareSTARGenome_skips_existing_index(tmp_path):
    """Should skip STAR index generation if Genome, SA, SAindex exist."""
    # Create fake genome_dir with required files
    genome_dir = tmp_path / "genome"
    genome_dir.mkdir()
    for f in ["Genome", "SA", "SAindex"]:
        (genome_dir / f).write_text("dummy")

    fasta_file = tmp_path / "genome.fa"
    gtf_file = tmp_path / "annotation.gtf"
    log_file = tmp_path / "log.txt"

    # Run function — should just log and return, not call STAR
    PrepareSTARGenome(str(genome_dir), str(fasta_file), str(gtf_file), 100, str(log_file))

    # Confirm log file was written
    assert log_file.exists()
    log_content = log_file.read_text()
    assert "STAR genome index already exists" in log_content


def test_PrepareSTARGenome_runs_star_when_missing(monkeypatch, tmp_path):
    """Should call STAR command when genome index is missing."""
    genome_dir = tmp_path / "new_index"
    fasta_file = tmp_path / "genome.fa"
    gtf_file = tmp_path / "annotation.gtf"
    log_file = tmp_path / "log.txt"

    # Create dummy fasta and gtf
    fasta_file.write_text(">chr1\nACTG")
    gtf_file.write_text("chr1\tsource\tgene\t1\t10\t.\t+\t.\tgene_id \"g1\";")

    # Mock subprocess.run to avoid actually running STAR
    called = {}

    def fake_run(cmd, check):
        called["cmd"] = cmd
        return 0

    monkeypatch.setattr(subprocess, "run", fake_run)

    # Run the function
    PrepareSTARGenome(str(genome_dir), str(fasta_file), str(gtf_file), 100, str(log_file))

    # Check that STAR command was called with expected params
    assert "STAR" in called["cmd"][0]
    assert "--genomeDir" in called["cmd"]
    assert str(genome_dir) in called["cmd"]

    # Log file should confirm generation started and finished
    content = log_file.read_text()
    assert "Running STAR command" in content
    assert "successfully generated" in content


def test_PrepareSTARGenome_handles_star_failure(monkeypatch, tmp_path, capsys):
    """Should exit with code 11 and log error if STAR fails."""
    genome_dir = tmp_path / "bad_index"
    fasta_file = tmp_path / "genome.fa"
    gtf_file = tmp_path / "annotation.gtf"
    log_file = tmp_path / "log.txt"

    fasta_file.write_text(">chr1\nACTG")
    gtf_file.write_text("chr1\t.\t.\t1\t10\t.\t+\t.\tgene_id \"g1\";")

    def fake_run_fail(cmd, check):
        raise subprocess.CalledProcessError(1, cmd, "simulated failure")

    monkeypatch.setattr(subprocess, "run", fake_run_fail)

    with pytest.raises(SystemExit) as e:
        PrepareSTARGenome(str(genome_dir), str(fasta_file), str(gtf_file), 75, str(log_file))

    # Ensure exit code 11
    assert e.value.code == 11

    out, err = capsys.readouterr()
    assert "STAR genome generation failed" in (out + err)

    # Check the log also contains the error
    assert "STAR genome generation failed" in log_file.read_text()