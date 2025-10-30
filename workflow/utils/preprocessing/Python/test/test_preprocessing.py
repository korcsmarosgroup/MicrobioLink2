import pytest
from preprocessing import LoadingConfiguration
from preprocessing import CheckingConfiguration
from preprocessing import CheckFASTQFiles
from preprocessing import PrepareSTARGenome
from preprocessing import RunSTARUnified
from preprocessing import _pair_fastqs

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
    
def test_pair_fastqs_merged_basic(tmp_path):
    """Should correctly pair R1, R2, and I1 files in merged layout."""
    create_fastq(tmp_path, "sampleA_R1.fastq.gz")
    create_fastq(tmp_path, "sampleA_R2.fastq.gz")
    create_fastq(tmp_path, "sampleA_I1.fastq.gz")
    create_fastq(tmp_path, "sampleB_R1.fastq.gz")
    create_fastq(tmp_path, "sampleB_R2.fastq.gz")

    result = _pair_fastqs(str(tmp_path), "merged")

    assert set(result.keys()) == {"sampleA", "sampleB"}

    # Check all expected reads exist
    for sample in result:
        assert "R1" in result[sample]
        assert "R2" in result[sample]
        assert result[sample]["R1"].endswith(f"{sample}_R1.fastq.gz")
        assert result[sample]["R2"].endswith(f"{sample}_R2.fastq.gz")

    # Optional I1 file should be detected for sampleA only
    assert "I1" in result["sampleA"]
    assert "I1" not in result["sampleB"]


def test_pair_fastqs_subdir_layout(tmp_path):
    """Should correctly find and pair files in per-sample subdirectories."""
    sample1 = tmp_path / "sample1"
    sample2 = tmp_path / "sample2"
    sample1.mkdir()
    sample2.mkdir()

    create_fastq(sample1, "sample1_R1.fastq")
    create_fastq(sample1, "sample1_R2.fastq")
    create_fastq(sample2, "sample2_R1.fastq.gz")
    create_fastq(sample2, "sample2_R2.fastq.gz")

    result = _pair_fastqs(str(tmp_path), "subdir")

    assert set(result.keys()) == {"sample1", "sample2"}

    # Verify both R1 and R2 for each sample
    for sample in ["sample1", "sample2"]:
        assert "R1" in result[sample]
        assert "R2" in result[sample]
        assert result[sample]["R1"].endswith(f"{sample}_R1.fastq" if sample == "sample1" else f"{sample}_R1.fastq.gz")
        assert result[sample]["R2"].endswith(f"{sample}_R2.fastq" if sample == "sample1" else f"{sample}_R2.fastq.gz")


def test_pair_fastqs_handles_missing_pairs(tmp_path):
    """Should include samples even if only one read file is present."""
    create_fastq(tmp_path, "lonely_R1.fastq")

    result = _pair_fastqs(str(tmp_path), "merged")

    assert "lonely" in result
    assert "R1" in result["lonely"]
    assert "R2" not in result["lonely"]


def test_pair_fastqs_empty_directory(tmp_path):
    """Should return empty dict if no FASTQ files exist."""
    result = _pair_fastqs(str(tmp_path), "merged")
    assert result == {}


def test_pair_fastqs_subdir_with_empty_folder(tmp_path):
    """Should ignore subdirectories with no FASTQ files."""
    empty = tmp_path / "empty_sample"
    filled = tmp_path / "filled_sample"
    empty.mkdir()
    filled.mkdir()

    create_fastq(filled, "filled_sample_R1.fastq")
    create_fastq(filled, "filled_sample_R2.fastq")

    result = _pair_fastqs(str(tmp_path), "subdir")

    assert set(result.keys()) == {"filled_sample"}
    assert "R1" in result["filled_sample"]
    assert "R2" in result["filled_sample"]


def test_pair_fastqs_ignores_non_fastq(tmp_path):
    """Should ignore files that are not FASTQ format."""
    (tmp_path / "random.txt").write_text("not a fastq")
    create_fastq(tmp_path, "valid_R1.fastq")

    result = _pair_fastqs(str(tmp_path), "merged")

    assert set(result.keys()) == {"valid"}
    for sample, reads in result.items():
        for read_file in reads.values():
            assert read_file.endswith(".fastq")


def test_pair_fastqs_mixed_layout(tmp_path):
    """Should respect layout type and not mix merged and subdir results."""
    # merged-level file
    create_fastq(tmp_path, "merged_R1.fastq")

    # subdir files
    sub = tmp_path / "subsample"
    sub.mkdir()
    create_fastq(sub, "subsample_R1.fastq")
    create_fastq(sub, "subsample_R2.fastq")

    # merged layout -> only top-level FASTQs
    result_merged = _pair_fastqs(str(tmp_path), "merged")
    assert set(result_merged.keys()) == {"merged"}

    # subdir layout -> only files within subdirectories
    result_subdir = _pair_fastqs(str(tmp_path), "subdir")
    assert set(result_subdir.keys()) == {"subsample"}
    
def test_run_star_microwell_executes_star(monkeypatch, tmp_path):
    """Should build and execute STAR command for microwell samples."""
    config = {
        "platform": "microwell",
        "input_dir": str(tmp_path),
        "genome_dir": str(tmp_path / "genome"),
        "threads": 4,
        "Fastq_file_format": "merged",
        "STAR_outdir": str(tmp_path / "star_out"),
    }
    log_file = tmp_path / "log.txt"
    create_fastq(tmp_path, "sample_R1.fastq")
    create_fastq(tmp_path, "sample_R2.fastq")

    called = {}
    def fake_run(cmd, check):
        called["cmd"] = cmd
        return 0

    monkeypatch.setattr(subprocess, "run", fake_run)
    RunSTARUnified(config, str(log_file), solo=False)

    assert "STAR" in called["cmd"][0]
    assert "--genomeDir" in called["cmd"]
    assert os.path.exists(config["STAR_outdir"])


def test_run_star_solo_executes_starsolo(monkeypatch, tmp_path):
    """Should run STARsolo when solo=True and platform=droplet."""
    config = {
        "platform": "droplet",
        "input_dir": str(tmp_path),
        "genome_dir": str(tmp_path / "genome"),
        "threads": 2,
        "Fastq_file_format": "merged",
        "STARsolo_outdir": str(tmp_path / "solo_out"),
        "STARsolo_params": {"soloType": "CB_UMI_Simple", "soloCBwhitelist": "None"},
    }
    log_file = tmp_path / "log.txt"
    create_fastq(tmp_path, "cell_R1.fastq")
    create_fastq(tmp_path, "cell_R2.fastq")

    called = {}
    def fake_run(cmd, check):
        called["cmd"] = cmd

    monkeypatch.setattr(subprocess, "run", fake_run)
    RunSTARUnified(config, str(log_file), solo=True)

    assert "STAR" in called["cmd"][0]
    assert "--soloType" in called["cmd"]
    assert os.path.exists(config["STARsolo_outdir"])


def test_run_star_skips_wrong_platform(monkeypatch, tmp_path, capsys):
    """Should skip execution when platform doesn't match mode."""
    config = {
        "platform": "droplet",  # wrong for solo=False
        "input_dir": str(tmp_path),
        "genome_dir": str(tmp_path / "genome"),
        "Fastq_file_format": "merged",
    }
    log_file = tmp_path / "log.txt"

    RunSTARUnified(config, str(log_file), solo=False)

    assert log_file.exists(), "Expected log file to be created"
    log_contents = log_file.read_text()

    # Assert the exact skip message was logged
    assert "Platform 'droplet' is not microwell. Skipping STAR." in log_contents

def test_run_star_skips_incomplete_pairs(monkeypatch, tmp_path):
    """Should skip samples missing R2 files."""
    config = {
        "platform": "microwell",
        "input_dir": str(tmp_path),
        "genome_dir": str(tmp_path / "genome"),
        "Fastq_file_format": "merged",
    }
    log_file = tmp_path / "log.txt"
    create_fastq(tmp_path, "sample_R1.fastq")  # missing R2

    called = {"ran": False}
    def fake_run(cmd, check):
        called["ran"] = True

    monkeypatch.setattr(subprocess, "run", fake_run)
    RunSTARUnified(config, str(log_file), solo=False)

    # STAR shouldn't run
    assert not called["ran"]
    assert "missing R1 or R2" in log_file.read_text()


def test_run_star_handles_subprocess_error(monkeypatch, tmp_path, capsys):
    """Should log error and continue if subprocess.run fails."""
    config = {
        "platform": "microwell",
        "input_dir": str(tmp_path),
        "genome_dir": str(tmp_path / "genome"),
        "Fastq_file_format": "merged",
    }
    log_file = tmp_path / "log.txt"
    create_fastq(tmp_path, "sample_R1.fastq")
    create_fastq(tmp_path, "sample_R2.fastq")

    def fake_run_fail(cmd, check):
        raise subprocess.CalledProcessError(1, cmd, "simulated failure")

    monkeypatch.setattr(subprocess, "run", fake_run_fail)

    RunSTARUnified(config, str(log_file), solo=False)

    output = capsys.readouterr()
    log_contents = log_file.read_text()

    assert "failed for sample" in (output.err + log_contents)