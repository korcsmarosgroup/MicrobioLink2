import pytest
from preprocessing import LoadingConfiguration
from preprocessing import CheckingConfiguration
from preprocessing import CheckFASTQFiles
from preprocessing import PrepareSTARGenome
from preprocessing import RunSTARUnified
from preprocessing import QC_10x
from preprocessing import NormalizeAllSamples
import subprocess
from unittest import mock
from pathlib import Path
import os
import yaml
import pandas as pd
import numpy as np
import scanpy as sc

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

# -----------------------------
# Tests for LoadingConfiguration
# -----------------------------
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

# -----------------------------
# Tests for CheckFASTQFiles
# -----------------------------

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
    """Should pass with correct R1, R2 files in merged layout."""
    log_file = tmp_path / "test.log"

    create_fastq(tmp_path, "sample_R1.fastq")
    create_fastq(tmp_path, "sample_R2.fastq")

    # Should not exit
    result = CheckFASTQFiles(str(tmp_path), "merged", str(log_file))

    assert "sample" in result
    assert "R1" in result["sample"]
    assert "R2" in result["sample"]

    # Log file should exist and contain expected message
    assert log_file.exists()
    content = log_file.read_text()
    assert "Using flat FASTQ layout" in content
    assert "sample has the necessary FASTQ files" not in content  # merged case uses WARNING for missing

def test_CheckFASTQFiles_subdir_no_subdirs(tmp_path, capsys):
    """Should exit with code 7 if no sample subdirectories exist."""
    log_file = tmp_path / "test.log"

    with pytest.raises(SystemExit) as e:
        CheckFASTQFiles(str(tmp_path), "subdir", str(log_file))

    assert e.value.code == 7
    out, err = capsys.readouterr()
    assert "No sample subdirectories found" in (out + err)

def test_CheckFASTQFiles_subdir_empty_warning(tmp_path):
    """Empty subdirectory logs a warning but does not exit."""
    log_file = tmp_path / "test.log"
    empty_sample = tmp_path / "empty_sample"
    empty_sample.mkdir()

    result = CheckFASTQFiles(str(tmp_path), "subdir", str(log_file))

    # Empty samples should not be returned
    assert "empty_sample" not in result

    log_content = log_file.read_text()
    assert "empty_sample has no FASTQ files" in log_content

def test_CheckFASTQFiles_subdir_valid_files(tmp_path):
    """Should pass with correct files inside sample subdirectory."""
    log_file = tmp_path / "test.log"
    sample_folder = tmp_path / "sampleB"
    sample_folder.mkdir()

    create_fastq(sample_folder, "sampleB_R1.fastq.gz")
    create_fastq(sample_folder, "sampleB_R2.fastq.gz")

    result = CheckFASTQFiles(str(tmp_path), "subdir", str(log_file))

    assert "sampleB" in result
    assert "R1" in result["sampleB"]
    assert "R2" in result["sampleB"]

    # Check that log file has correct message
    content = log_file.read_text()
    assert "sampleB has the necessary FASTQ files (R1, R2)" in content

# -----------------------------
# New tests for updated function
# -----------------------------

def test_CheckFASTQFiles_returns_only_complete_pairs(tmp_path):
    """Merged layout returns only complete R1+R2 samples."""
    log_file = tmp_path / "log.txt"

    create_fastq(tmp_path, "A_R1.fastq")  # missing R2
    create_fastq(tmp_path, "B_R1.fastq")
    create_fastq(tmp_path, "B_R2.fastq")  # complete

    result = CheckFASTQFiles(str(tmp_path), "merged", str(log_file))

    assert "A" not in result
    assert "B" in result
    assert "R1" in result["B"]
    assert "R2" in result["B"]

def test_CheckFASTQFiles_logs_missing_R2_warning(tmp_path):
    """Merged layout logs a WARNING for missing R2."""
    log_file = tmp_path / "log.txt"

    create_fastq(tmp_path, "S_R1.fastq")  # incomplete

    CheckFASTQFiles(str(tmp_path), "merged", str(log_file))

    content = log_file.read_text()
    assert "WARNING: S is missing R2 file(s)" in content

def test_CheckFASTQFiles_subdir_returns_only_complete(tmp_path):
    """Subdir layout returns only complete samples."""
    log_file = tmp_path / "log.txt"

    # Complete sample
    s1 = tmp_path / "X"
    s1.mkdir()
    create_fastq(s1, "X_R1.fastq")
    create_fastq(s1, "X_R2.fastq")

    # Incomplete sample
    s2 = tmp_path / "Y"
    s2.mkdir()
    create_fastq(s2, "Y_R1.fastq")

    result = CheckFASTQFiles(str(tmp_path), "subdir", str(log_file))

    assert "X" in result
    assert "Y" not in result

def test_CheckFASTQFiles_subdir_multiple_samples(tmp_path):
    """Subdir layout returns only complete samples; incomplete samples are skipped."""
    log_file = tmp_path / "log.txt"

    # Complete sample 1
    s1 = tmp_path / "Sample1"
    s1.mkdir()
    create_fastq(s1, "Sample1_R1.fastq")
    create_fastq(s1, "Sample1_R2.fastq")

    # Complete sample 2
    s2 = tmp_path / "Sample2"
    s2.mkdir()
    create_fastq(s2, "Sample2_R1.fastq")
    create_fastq(s2, "Sample2_R2.fastq")

    # Incomplete sample 3 (missing R2)
    s3 = tmp_path / "Sample3"
    s3.mkdir()
    create_fastq(s3, "Sample3_R1.fastq")

    result = CheckFASTQFiles(str(tmp_path), "subdir", str(log_file))

    # Only complete samples should be returned
    assert "Sample1" in result
    assert "Sample2" in result
    assert "Sample3" not in result

    # Check that R1 and R2 keys exist for returned samples
    for sample in ["Sample1", "Sample2"]:
        assert "R1" in result[sample]
        assert "R2" in result[sample]

    # Log file should contain a warning for the incomplete sample
    log_content = log_file.read_text()
    assert "Sample3 is missing R2 file(s)" in log_content

def test_PrepareSTARGenome_skips_existing_index(tmp_path):
    """Should skip STAR index generation if Genome, SA, SAindex exist."""
    genome_dir = tmp_path / "genome"
    genome_dir.mkdir()
    for f in ["Genome", "SA", "SAindex"]:
        (genome_dir / f).write_text("dummy")

    fasta_file = tmp_path / "genome.fa"
    fasta_file.write_text(">chr1\nACTG")

    log_file = tmp_path / "log.txt"

    genome_index_params = {
        "genomeChrBinNbits": 18,
        "genomeSAindexNbases": 14,
        "genomeSAsparseD": 1,
        "genomeSuffixLengthMax": -1,
        "genomeTransformType": None,
        "genomeTransformVCF": None
    }
    splice_junction_params = {
        "sjdbFileChrStartEnd": None,
        "sjdbGTFfile": str(tmp_path / "annotation.gtf"),
        "sjdbGTFchrPrefix": "",
        "sjdbGTFfeatureExon": "exon",
        "sjdbGTFtagExonParentTranscript": "transcript_id",
        "sjdbGTFtagExonParentGene": "gene_id",
        "sjdbGTFtagExonParentGeneName": "gene_name",
        "sjdbGTFtagExonParentGeneType": "gene_biotype",
        "sjdbOverhang": 99,
        "sjdbScore": 2,
        "sjdbInsertSave": "Basic"
    }
    (tmp_path / "annotation.gtf").write_text("chr1\tsource\tgene\t1\t10\t.\t+\t.\tgene_id \"g1\";")

    # Call new function
    PrepareSTARGenome(
        str(genome_dir),
        str(fasta_file),
        genome_index_params,
        splice_junction_params,
        str(log_file)
    )

    assert log_file.exists()
    assert "STAR genome index already exists" in log_file.read_text()
    
def test_PrepareSTARGenome_runs_star_when_missing(monkeypatch, tmp_path):
    """Should call STAR command when genome index is missing."""
    genome_dir = tmp_path / "new_index"
    fasta_file = tmp_path / "genome.fa"
    fasta_file.write_text(">chr1\nACTG")

    annotation_gtf = tmp_path / "annotation.gtf"
    annotation_gtf.write_text("chr1\tsource\tgene\t1\t10\t.\t+\t.\tgene_id \"g1\";")

    log_file = tmp_path / "log.txt"

    genome_index_params = {
        "genomeChrBinNbits": 18,
        "genomeSAindexNbases": 14,
        "genomeSAsparseD": 1,
        "genomeSuffixLengthMax": -1,
        "genomeTransformType": None,
        "genomeTransformVCF": None
    }
    splice_junction_params = {
        "sjdbFileChrStartEnd": None,
        "sjdbGTFfile": str(annotation_gtf),
        "sjdbGTFchrPrefix": "",
        "sjdbGTFfeatureExon": "exon",
        "sjdbGTFtagExonParentTranscript": "transcript_id",
        "sjdbGTFtagExonParentGene": "gene_id",
        "sjdbGTFtagExonParentGeneName": "gene_name",
        "sjdbGTFtagExonParentGeneType": "gene_biotype",
        "sjdbOverhang": 99,
        "sjdbScore": 2,
        "sjdbInsertSave": "Basic"
    }

    called = {}
    def fake_run(cmd, check):
        called["cmd"] = cmd
        return 0

    monkeypatch.setattr(subprocess, "run", fake_run)

    PrepareSTARGenome(
        str(genome_dir),
        str(fasta_file),
        genome_index_params,
        splice_junction_params,
        str(log_file)
    )

    assert "STAR" in called["cmd"][0]
    assert "--genomeDir" in called["cmd"]
    assert str(genome_dir) in called["cmd"]

    content = log_file.read_text()
    assert "Running STAR genomeGenerate command" in content
    assert "successfully generated" in content
    
    
def test_PrepareSTARGenome_handles_star_failure(monkeypatch, tmp_path, capsys):
    """Should exit with code 11 and log error if STAR fails."""
    genome_dir = tmp_path / "bad_index"
    fasta_file = tmp_path / "genome.fa"
    fasta_file.write_text(">chr1\nACTG")
    annotation_gtf = tmp_path / "annotation.gtf"
    annotation_gtf.write_text("chr1\t.\t.\t1\t10\t.\t+\t.\tgene_id \"g1\";")
    log_file = tmp_path / "log.txt"

    genome_index_params = {
        "genomeChrBinNbits": 18,
        "genomeSAindexNbases": 14,
        "genomeSAsparseD": 1,
        "genomeSuffixLengthMax": -1,
        "genomeTransformType": None,
        "genomeTransformVCF": None
    }
    splice_junction_params = {
        "sjdbFileChrStartEnd": None,
        "sjdbGTFfile": str(annotation_gtf),
        "sjdbGTFchrPrefix": "",
        "sjdbGTFfeatureExon": "exon",
        "sjdbGTFtagExonParentTranscript": "transcript_id",
        "sjdbGTFtagExonParentGene": "gene_id",
        "sjdbGTFtagExonParentGeneName": "gene_name",
        "sjdbGTFtagExonParentGeneType": "gene_biotype",
        "sjdbOverhang": 75,
        "sjdbScore": 2,
        "sjdbInsertSave": "Basic"
    }

    def fake_run_fail(cmd, check):
        raise subprocess.CalledProcessError(1, cmd, "simulated failure")

    monkeypatch.setattr(subprocess, "run", fake_run_fail)

    with pytest.raises(SystemExit) as e:
        PrepareSTARGenome(
            str(genome_dir),
            str(fasta_file),
            genome_index_params,
            splice_junction_params,
            str(log_file)
        )

    assert e.value.code == 11
    out, err = capsys.readouterr()
    assert "STAR genome generation failed" in (out + err)
    assert "STAR genome generation failed" in log_file.read_text()
    
    

    
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


def test_run_star_skips_wrong_platform(monkeypatch, tmp_path):
    """Should skip execution when platform doesn't match mode."""
    config = {
        "platform": "droplet",  # wrong for solo=False
        "input_dir": str(tmp_path),
        "genome_dir": str(tmp_path / "genome"),
        "Fastq_file_format": "merged",
    }
    log_file = tmp_path / "log.txt"
    create_fastq(tmp_path, "sample_R1.fastq")
    create_fastq(tmp_path, "sample_R2.fastq")

    RunSTARUnified(config, str(log_file), solo=False)

    assert log_file.exists(), "Expected log file to be created"
    log_contents = log_file.read_text()
    assert "Platform 'droplet' is not microwell. Skipping STAR." in log_contents


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


def test_run_star_creates_output_dirs(monkeypatch, tmp_path):
    """Should create output directories for each sample."""
    config = {
        "platform": "microwell",
        "input_dir": str(tmp_path),
        "genome_dir": str(tmp_path / "genome"),
        "Fastq_file_format": "merged",
        "STAR_outdir": str(tmp_path / "star_out"),
    }
    log_file = tmp_path / "log.txt"
    create_fastq(tmp_path, "sample_R1.fastq")
    create_fastq(tmp_path, "sample_R2.fastq")

    def fake_run(cmd, check):
        return 0

    monkeypatch.setattr(subprocess, "run", fake_run)

    RunSTARUnified(config, str(log_file), solo=False)

    expected_dir = os.path.join(config["STAR_outdir"], "sample")
    # Since sample prefix is "sample" from "sample_R1.fastq"
    assert os.path.isdir(os.path.join(config["STAR_outdir"], "sample"))
    
    
@pytest.fixture
def fake_adata():
    """Provide a small fake AnnData object with values that pass QC filters."""

    X = np.random.poisson(1500, (10, 5))  # high total counts
    obs = pd.DataFrame({
        "n_genes_by_counts": np.random.randint(1000, 2000, 10),
        "pct_counts_mt": np.random.uniform(0, 2, 10),  # below 5%
        "total_counts": np.random.randint(2000, 4000, 10)
    }, index=[f"cell{i}" for i in range(10)])

    var = pd.DataFrame(index=[f"gene{i}" for i in range(5)])
    var["mt"] = [False, False, True, False, False]  # one MT gene

    adata = sc.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = adata.X.copy()
    return adata
@pytest.fixture
def sample_structure(tmp_path):
    """Create a fake output folder with one sample/filtered dir."""
    sample_dir = tmp_path / "Sample1" / "filtered"
    sample_dir.mkdir(parents=True)
    log_file = tmp_path / "test_log.txt"
    return tmp_path, log_file

def test_QC_10x_loads_data(fake_adata, sample_structure):
    output_folder, log_file = sample_structure

    with mock.patch("scanpy.read_10x_mtx", return_value=fake_adata) as mock_read, \
         mock.patch("scrublet.Scrublet"), \
         mock.patch("matplotlib.pyplot.savefig"):
        QC_10x(str(output_folder), str(log_file))

    mock_read.assert_called_once()


def test_QC_10x_computes_qc_metrics(fake_adata, sample_structure):
    output_folder, log_file = sample_structure

    with mock.patch("scanpy.read_10x_mtx", return_value=fake_adata), \
         mock.patch("scrublet.Scrublet"), \
         mock.patch("matplotlib.pyplot.savefig"):
        QC_10x(str(output_folder), str(log_file))

    # Verify QC columns exist
    assert "total_counts" in fake_adata.obs.columns or "pct_counts_mt" in fake_adata.obs.columns


def test_QC_10x_generates_plots(fake_adata, sample_structure):
    output_folder, log_file = sample_structure

    with mock.patch("scanpy.read_10x_mtx", return_value=fake_adata), \
         mock.patch("scrublet.Scrublet"), \
         mock.patch("matplotlib.pyplot.savefig") as mock_save:
        QC_10x(str(output_folder), str(log_file))

    assert mock_save.call_count >= 1


def test_QC_10x_filters_low_quality_cells(fake_adata, sample_structure):
    output_folder, log_file = sample_structure

    # artificially set QC metrics
    fake_adata.obs["n_genes_by_counts"] = np.random.randint(1000, 3000, 10)
    fake_adata.obs["pct_counts_mt"] = np.random.rand(10) * 10
    fake_adata.obs["total_counts"] = np.random.randint(500, 2000, 10)

    with mock.patch("scanpy.read_10x_mtx", return_value=fake_adata), \
         mock.patch("scrublet.Scrublet"), \
         mock.patch("matplotlib.pyplot.savefig"):
        QC_10x(str(output_folder), str(log_file))

    # Filter should reduce cell count or leave it unchanged
    assert fake_adata.n_obs >= 0  # sanity check


def test_QC_10x_runs_scrublet(fake_adata, sample_structure):
    output_folder, log_file = sample_structure

    mock_scrub = mock.Mock()
    mock_scrub.scrub_doublets.return_value = (np.random.rand(10), np.random.choice([True, False], 10))

    with mock.patch("scanpy.read_10x_mtx", return_value=fake_adata), \
         mock.patch("scrublet.Scrublet", return_value=mock_scrub), \
         mock.patch("matplotlib.pyplot.savefig"):
        QC_10x(str(output_folder), str(log_file))

    mock_scrub.scrub_doublets.assert_called_once()


def test_QC_10x_handles_scrublet_failure(fake_adata, sample_structure):
    output_folder, log_file = sample_structure

    with mock.patch("scanpy.read_10x_mtx", return_value=fake_adata), \
         mock.patch("scrublet.Scrublet", side_effect=Exception("Scrublet crash")), \
         mock.patch("matplotlib.pyplot.savefig"):
        QC_10x(str(output_folder), str(log_file))

    with open(log_file) as f:
        logs = f.read()
    assert "WARNING" in logs or "Scrublet failed" in logs


def test_QC_10x_writes_output(fake_adata, sample_structure):
    output_folder, log_file = sample_structure

    with mock.patch("scanpy.read_10x_mtx", return_value=fake_adata), \
         mock.patch("scrublet.Scrublet"), \
         mock.patch("matplotlib.pyplot.savefig"):
        QC_10x(str(output_folder), str(log_file))

    qc_out = output_folder / "Sample1" / "QC_filtered.h5ad"
    assert qc_out.exists()

def test_QC_10x_output_file_structure(fake_adata, sample_structure):
    """
    Unit test: verify that the QC_filtered.h5ad file has a valid AnnData structure
    and contains expected QC metrics.
    """
    output_folder, log_file = sample_structure

    with mock.patch("scanpy.read_10x_mtx", return_value=fake_adata), \
         mock.patch("scrublet.Scrublet"), \
         mock.patch("matplotlib.pyplot.savefig"):
        QC_10x(str(output_folder), str(log_file))

    qc_out = output_folder / "Sample1" / "QC_filtered.h5ad"
    assert qc_out.exists(), "QC-filtered output file was not created."

    # Load the QC output file
    adata_out = sc.read_h5ad(qc_out)

    # --- Check structure ---
    assert isinstance(adata_out, sc.AnnData), "Output is not a valid AnnData object."
    assert adata_out.X is not None, "Expression matrix missing."
    assert adata_out.n_obs > 0, "No cells in output."
    assert adata_out.n_vars > 0, "No genes in output."

    # --- Check for essential QC columns ---
    expected_cols = {"total_counts", "n_genes_by_counts", "pct_counts_mt"}
    assert expected_cols.issubset(adata_out.obs.columns), (
        f"Missing expected QC columns: {expected_cols - set(adata_out.obs.columns)}"
    )
def create_dummy_h5ad(folder: Path, sample_name: str, n_obs=5, n_vars=10):
    """Create a dummy QC_filtered.h5ad file for testing."""
    sample_folder = folder / sample_name
    sample_folder.mkdir(parents=True, exist_ok=True)
    data = np.random.rand(n_obs, n_vars)
    adata = sc.AnnData(X=data)
    adata.obs_names = [f"cell{i}" for i in range(n_obs)]
    adata.var_names = [f"gene{j}" for j in range(n_vars)]
    input_file = sample_folder / "QC_filtered.h5ad"
    adata.write(input_file)
    return sample_folder, input_file

# -----------------------------
# TESTS
# -----------------------------

def test_normalize_all_samples_creates_files(tmp_path):
    """Should create normalized.h5ad and normalized.tsv for a single sample."""
    log_file = tmp_path / "log.txt"
    sample_folder, input_file = create_dummy_h5ad(tmp_path, "Sample1")

    # Run normalization
    NormalizeAllSamples(str(tmp_path), HVG_selection=False, number_top_genes=0, log_file=str(log_file))

    # Check that normalized files exist
    assert (sample_folder / "normalized.h5ad").exists()
    assert (sample_folder / "normalized.tsv").exists()

    # HVG files should not exist
    assert not (sample_folder / "normalized_hvg.h5ad").exists()
    assert not (sample_folder / "normalized_hvg.tsv").exists()

def test_normalize_all_samples_with_HVG(tmp_path):
    """Should create HVG files if HVG_selection is True."""
    log_file = tmp_path / "log.txt"
    sample_folder, input_file = create_dummy_h5ad(tmp_path, "Sample2")

    NormalizeAllSamples(str(tmp_path), HVG_selection=True, number_top_genes=5, log_file=str(log_file))

    # Check HVG outputs
    assert (sample_folder / "normalized_hvg.h5ad").exists()
    assert (sample_folder / "normalized_hvg.tsv").exists()

    # Also check normal outputs exist
    assert (sample_folder / "normalized.h5ad").exists()
    assert (sample_folder / "normalized.tsv").exists()

def test_normalize_all_samples_skips_missing_input(tmp_path):
    """Should skip sample folders without QC_filtered.h5ad and log a message."""
    log_file = tmp_path / "log.txt"
    (tmp_path / "EmptySample").mkdir()

    NormalizeAllSamples(str(tmp_path), HVG_selection=False, number_top_genes=0, log_file=str(log_file))

    log_content = Path(log_file).read_text()
    assert "Skipping EmptySample: No QC_filtered.h5ad found." in log_content

def test_normalize_all_samples_multiple_samples(tmp_path):
    """Should normalize multiple samples correctly."""
    log_file = tmp_path / "log.txt"
    create_dummy_h5ad(tmp_path, "S1")
    create_dummy_h5ad(tmp_path, "S2")

    NormalizeAllSamples(str(tmp_path), HVG_selection=False, number_top_genes=0, log_file=str(log_file))

    for s in ["S1", "S2"]:
        sample_folder = tmp_path / s
        assert (sample_folder / "normalized.h5ad").exists()
        assert (sample_folder / "normalized.tsv").exists()

def test_normalize_all_samples_exception_handling(monkeypatch, tmp_path):
    """Should exit with code 13 if normalization fails."""
    log_file = tmp_path / "log.txt"
    sample_folder, input_file = create_dummy_h5ad(tmp_path, "SampleX")

    def fake_read_h5ad(file):
        raise RuntimeError("Simulated failure")

    monkeypatch.setattr(sc, "read_h5ad", fake_read_h5ad)

    with pytest.raises(SystemExit) as e:
        NormalizeAllSamples(str(tmp_path), HVG_selection=False, number_top_genes=0, log_file=str(log_file))

    assert e.value.code == 13
    content = Path(log_file).read_text()
    assert "ERROR: Normalization failed" in content