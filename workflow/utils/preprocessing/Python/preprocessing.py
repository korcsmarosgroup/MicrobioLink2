# This module/script will take either FASTQ files from a scRNA-seq analysis
# or it will take a raw count matrix from a scRNA-seq analysis


from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
from time import strftime
import scrublet as scr
import pandas as pd
import scanpy as sc
import numpy as np
import subprocess
import shutil
import gzip
import yaml
import glob
import sys
import os


def get_log_file(config_folder):
    """
        Create a timestamped log filename.
        Example:
            preprocessing_2025-11-18_14-26-30.log
    """
    timestamp = strftime("%Y-%m-%d_%H-%M-%S")
    return f"{config_folder}/preprocessing_{timestamp}.log"


def _log(message, log_file):
    """
        Log a message with a timestamp.
    """
    with open(log_file, "a") as log:
        log.write(f"[{strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")


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


def backup_config(config_folder, output_folder):
    """
        Copy YAML files (*.yaml, *.yml) from config_folder into output_folder.
        Adds a timestamp to each file name instead of creating a backup folder.
        Example: config.yaml → config_2025-02-10_15-32-18.yaml
    """
    timestamp = strftime("%Y-%m-%d_%H-%M-%S")

    os.makedirs(output_folder, exist_ok=True)

    for filename in os.listdir(config_folder):
        if filename.endswith(".yaml") or filename.endswith(".yml"):
            src = os.path.join(config_folder, filename)

            # Split filename into name + extension
            name, ext = os.path.splitext(filename)

            # New filename with timestamp
            new_filename = f"{name}_{timestamp}{ext}"
            dst = os.path.join(output_folder, new_filename)

            shutil.copy2(src, dst)

    return output_folder


def CheckFASTQFiles(input_folder, Fastq_file_format, log_file):
    """
        Validate FASTQ inputs and return only samples with complete paired-end data.

        This function scans an input directory that contains either:
        • A flat collection of FASTQ files ("merged" layout), or
        • Multiple per-sample subdirectories ("subdir" layout).

        The function:
            - Identifies FASTQ files for each sample (R1, R2, and optionally I1)
            - Logs warnings when samples are incomplete (missing R1 or R2)
            - Logs errors and exits when the folder structure is invalid
            - Produces a dictionary containing **only samples that have both R1 and R2 FASTQs**

        The returned dictionary is suitable for downstream tools that require
        complete paired FASTQ inputs (e.g., STAR, trimming tools, quantification).

        Args:
            input_folder:
                Path to the main input directory containing FASTQ files or sample folders.
            Fastq_file_format:
                Expected layout of the FASTQ input. Must be:
                    - "merged": all FASTQs in one folder
                    - "subdir": each sample in its own subdirectory
            log_file:
                Path to the log file where warnings and information messages are written.

        Error codes:
            ERROR CODE 5: Invalid Fastq_file_format argument.
            ERROR CODE 6: No FASTQ files found in merged layout.
            ERROR CODE 7: No sample subdirectories found in subdir layout.

        Returns:
            dict:
                A dictionary mapping sample prefixes to dictionaries of file paths, e.g.:

                {
                    "SampleA": {
                        "R1": "/path/to/SampleA_R1.fastq.gz",
                        "R2": "/path/to/SampleA_R2.fastq.gz",
                        "I1": "/path/to/SampleA_I1.fastq.gz"   # optional
                    },
                    "SampleB": {
                        "R1": "...",
                        "R2": "..."
                    }
                }

                Only samples that contain **both R1 and R2** files are included.
    """

    # Validate Fastq_file_format
    if Fastq_file_format not in ("merged", "subdir"):
        sys.stderr.write(
            f"ERROR: Invalid Fastq_file_format '{Fastq_file_format}'. "
            "Must be 'merged' or 'subdir'.\n"
        )
        sys.exit(5)

    # ------------------------------------------------------------------
    # Internal pairing logic (merged version of _pair_fastqs)
    # ------------------------------------------------------------------
    def _collect_fastqs(input_dir, layout):
        fastq_files = []

        if layout == "merged":
            fastq_files = (
                glob.glob(os.path.join(input_dir, "*.fastq")) +
                glob.glob(os.path.join(input_dir, "*.fastq.gz"))
            )

        elif layout == "subdir":
            subdirs = [
                d for d in glob.glob(os.path.join(input_dir, "*"))
                if os.path.isdir(d)
            ]
            for sd in subdirs:
                fastq_files.extend(
                    glob.glob(os.path.join(sd, "*.fastq")) +
                    glob.glob(os.path.join(sd, "*.fastq.gz"))
                )

        samples = {}
        for f in fastq_files:
            base = os.path.basename(f)

            if "_R1" in base:
                prefix = base.split("_R1")[0]
                samples.setdefault(prefix, {})["R1"] = f

            elif "_R2" in base:
                prefix = base.split("_R2")[0]
                samples.setdefault(prefix, {})["R2"] = f

            elif "_I1" in base:
                prefix = base.split("_I1")[0]
                samples.setdefault(prefix, {})["I1"] = f

        return samples

    # ------------------------------------------------------------------
    # CASE 1: MERGED layout
    # ------------------------------------------------------------------
    if Fastq_file_format == "merged":
        fastq_files = (
            glob.glob(os.path.join(input_folder, "*.fastq")) +
            glob.glob(os.path.join(input_folder, "*.fastq.gz"))
        )

        if len(fastq_files) == 0:
            sys.stderr.write(f"ERROR: No FASTQ files found in {input_folder}\n")
            sys.exit(6)

        _log("Using merged FASTQ layout", log_file)

        samples = _collect_fastqs(input_folder, "merged")

        # Validate and log
        for sample, files in samples.items():
            found = {"R1": "R1" in files, "R2": "R2" in files}
            missing = [k for k, v in found.items() if not v]

            if missing:
                _log(
                    f"WARNING: {sample} is missing {', '.join(missing)} file(s).",
                    log_file
                )
                continue

            if len(files) < 2:
                _log(
                    f"WARNING: {sample} should have at least R1 and R2, "
                    f"but only {len(files)} FASTQs found.",
                    log_file
                )
                continue

    # ------------------------------------------------------------------
    # CASE 2: SUBDIR layout
    # ------------------------------------------------------------------
    elif Fastq_file_format == "subdir":
        subdirs = [
            d for d in glob.glob(os.path.join(input_folder, "*"))
            if os.path.isdir(d)
        ]

        if len(subdirs) == 0:
            sys.stderr.write(
                f"ERROR: No sample subdirectories found in {input_folder}\n"
            )
            sys.exit(7)

        _log("Using per-sample subdirectory layout.\n", log_file)

        samples = _collect_fastqs(input_folder, "subdir")

        # Validate and log
        for folder in subdirs:
            sample_name = os.path.basename(folder)
            fastqs = (
                glob.glob(os.path.join(folder, "*.fastq")) +
                glob.glob(os.path.join(folder, "*.fastq.gz"))
            )

            if len(fastqs) == 0:
                _log(f"WARNING: {sample_name} has no FASTQ files", log_file)
                continue

            found = {"R1": False, "R2": False}
            for f in fastqs:
                fname = os.path.basename(f)
                for tag in found:
                    if tag in fname:
                        found[tag] = True

            missing = [k for k, v in found.items() if not v]
            if missing:
                _log(
                    f"WARNING: {sample_name} is missing {', '.join(missing)} file(s).",
                    log_file
                )
                continue

            if len(fastqs) < 2:
                _log(
                    f"WARNING: {sample_name} should have at least 2 FASTQ files, "
                    f"but {len(fastqs)} found.",
                    log_file
                )
                continue

            _log(f"{sample_name} has the necessary FASTQ files (R1, R2)", log_file)

    # ------------------------------------------------------------------
    # FINAL STEP — RETURN ONLY SAMPLES WITH BOTH R1 & R2
    # ------------------------------------------------------------------
    complete_samples = {
        sample: files
        for sample, files in samples.items()
        if "R1" in files and "R2" in files
    }

    return complete_samples
    '''
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
            ERROR CODE 11: STAR genome generation failed.

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
            required_files = {"Genome", "SA", "SAindex"} # TODO: Check what file extensions needed to be inside a STAR index folder
            if required_files.issubset(set(existing_files)):
                _log(f"STAR genome index already exists at {genome_dir}. Skipping generation.", log_file)
                return

        # ------------------------------------------------------------
        # Create genome directory if missing
        # ------------------------------------------------------------
        else:
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
            sys.exit(11)
    '''


def PrepareSTARGenome(genome_dir, fasta_file, genome_index_params, splice_junction_params, log_file):
    """
        Preparing the STAR genome index directory using full genome indexing and splice junction parameters.

        Checks whether a valid STAR genome index already exists in the specified
        genome directory. If all required index files are present, it skips index generation.
        Otherwise, it generates a new STAR genome index using the provided FASTA, genome indexing
        parameters, and splice junction parameters.

        Args:
            genome_dir (str): Directory where the STAR genome index will be stored.
            fasta_file (str): Path to the reference genome FASTA file.
            genome_index_params (dict): Genome indexing parameters from YAML (17.5).
            splice_junction_params (dict): Splice junction parameters from YAML (17.6).
            log_file (str): Path of the log file.

        Error codes:
            ERROR CODE 11: STAR genome generation failed.

        Returns:
            None
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
    else:
        _log(f"STAR genome index not found. Generating new index at {genome_dir}", log_file)
        os.makedirs(genome_dir, exist_ok=True)

    # ------------------------------------------------------------
    # Build STAR genomeGenerate command
    # ------------------------------------------------------------
    cmd = [
        "STAR",
        "--runThreadN", "8",
        "--runMode", "genomeGenerate",
        "--genomeDir", genome_dir,
        "--genomeFastaFiles", fasta_file,
        "--sjdbGTFfile", splice_junction_params["sjdbGTFfile"],
        "--sjdbOverhang", str(splice_junction_params["sjdbOverhang"]),
        "--genomeChrBinNbits", str(genome_index_params["genomeChrBinNbits"]),
        "--genomeSAindexNbases", str(genome_index_params["genomeSAindexNbases"]),
        "--genomeSAsparseD", str(genome_index_params["genomeSAsparseD"]),
        "--genomeSuffixLengthMax", str(genome_index_params["genomeSuffixLengthMax"])
    ]

    # Optional genome transform
    if genome_index_params["genomeTransformType"] is not None:
        cmd += ["--genomeTransformType", genome_index_params["genomeTransformType"]]
        if genome_index_params["genomeTransformVCF"]:
            cmd += ["--genomeTransformVCF", genome_index_params["genomeTransformVCF"]]

    # Optional splice junction file
    if splice_junction_params["sjdbFileChrStartEnd"]:
        cmd += ["--sjdbFileChrStartEnd", splice_junction_params["sjdbFileChrStartEnd"]]

    # Splice junction tags
    cmd += [
        "--sjdbGTFchrPrefix", splice_junction_params["sjdbGTFchrPrefix"],
        "--sjdbGTFfeatureExon", splice_junction_params["sjdbGTFfeatureExon"],
        "--sjdbGTFtagExonParentTranscript", splice_junction_params["sjdbGTFtagExonParentTranscript"],
        "--sjdbGTFtagExonParentGene", splice_junction_params["sjdbGTFtagExonParentGene"],
        "--sjdbGTFtagExonParentGeneName", splice_junction_params["sjdbGTFtagExonParentGeneName"],
        "--sjdbGTFtagExonParentGeneType", splice_junction_params["sjdbGTFtagExonParentGeneType"],
        "--sjdbScore", str(splice_junction_params["sjdbScore"]),
        "--sjdbInsertSave", splice_junction_params["sjdbInsertSave"]
    ]

    # ------------------------------------------------------------
    # Run STAR genome generation
    # ------------------------------------------------------------
    try:
        _log(f"Running STAR genomeGenerate command: {' '.join(cmd)}", log_file)
        subprocess.run(cmd, check=True)
        _log(f"STAR genome index successfully generated at {genome_dir}", log_file)

    except subprocess.CalledProcessError as e:
        error_msg = f"ERROR: STAR genome generation failed: {e}"
        _log(error_msg, log_file)
        sys.stderr.write(f"ERROR MESSAGE: STAR genome generation failed: {e}\n")
        sys.exit(11)
        

def RunSTARUnified(configuration, log_file):
    """
        Run STAR or STARsolo depending on the experimental platform and configuration settings.

    This function executes the STAR aligner for microwell-based experiments or STARsolo
    for droplet-based single-cell RNA-seq data. It automatically detects input FASTQ pairs,
    constructs appropriate STAR command-line arguments, and logs progress and errors.

    The function supports both gzipped and uncompressed FASTQ files, and dynamically
    includes custom STAR/STARsolo parameters defined in the configuration.

    Args:
        configuration (dict): Dictionary containing run parameters and paths. Expected keys:
            - "platform" (str): Sequencing platform ("10x", "SmartSeq" or "DropSeq").
            - "input_dir" (str): Directory containing FASTQ files.
            - "genome_dir" (str): Directory of the STAR genome index.
            - "threads" (int, optional): Number of CPU threads for STAR. Defaults to 8.
            - "Fastq_file_format" (str, optional): Layout of FASTQ files ("merged" or "subdir").
            - "STAR_outdir" (str, optional): Output directory for STAR alignments.
            - "STAR_params" (dict, optional): Additional STAR command-line parameters.
            - "STARsolo_outdir" (str, optional): Output directory for STARsolo.
            - "STARsolo_params" (dict, optional): Additional STARsolo command-line parameters.

        log_file (str): Path to a log file for recording messages and errors.
        
        solo (bool, optional): 
            - If True: Runs STARsolo (droplet-based single-cell data).
            - If False: Runs standard STAR alignment (microwell-based data).
            Defaults to False.
            
            Returns:
                

    Notes:
        - Requires STAR to be installed and accessible in the system PATH.
        - FASTQ files must contain '_R1' and '_R2' in their filenames.
        - Skips samples missing either R1 or R2.
        - Use   

    """

    platform = configuration.get("platform", "").strip()
    input_dir = configuration.get("input_dir")
    output_dir = configuration.get("output_dir")
    genome_dir = configuration.get("genome_dir")
    threads = configuration.get("threads", 4)
    layout = configuration.get("Fastq_file_format", "merged")

    # ============================================================
    # Determine STAR mode based on platform
    # ============================================================

    if platform == "SmartSeq":
        mode = "STAR"
        params = configuration.get("star_params", {})
        output_dir = os.path.join(output_dir, "STAR_out")

    elif platform == "10x":
        mode = "STARsolo"
        params = configuration.get("STARsolo_params_10x", {})
        output_dir = os.path.join(output_dir, "STARsolo_10x_out")

    elif platform == "DropSeq":
        mode = "STARsolo"
        params = configuration.get("STARsolo_params_DropSeq", {})
        output_dir = os.path.join(output_dir, "STARsolo_DropSeq_out")

    else:
        _log(f"ERROR: Platform '{platform}' not recognized. Use SmartSeq, 10x, or DropSeq.", log_file)
        return

    _log(f"Selected mode: {mode} for platform: {platform}", log_file)

    # ============================================================
    # Validate FASTQ files
    # ============================================================
    samples = CheckFASTQFiles(input_dir, layout, log_file)

    # ============================================================
    # Run STAR/STARsolo
    # ============================================================

    for sample, files in samples.items():

        sample_out = os.path.join(output_dir, sample)
        os.makedirs(sample_out, exist_ok=True)

        cmd = [
            "STAR",
            "--runThreadN", str(threads),
            "--genomeDir", genome_dir,
            "--readFilesIn", files["R1"], files["R2"],
            "--outFileNamePrefix", os.path.join(sample_out, "")
        ]

        # Gzipped reads
        if files["R1"].endswith(".gz") or files["R2"].endswith(".gz"):
            cmd += ["--readFilesCommand", "zcat"]

        # Append STAR or STARsolo parameters
        for param, value in params.items():
            cmd.append(f"--{param}")
            # Only append non-empty values
            if value not in (None, "", "None"):
                cmd.append(str(value))

        # Log and run
        try:
            _log(f"Running {mode} for {sample}: {' '.join(cmd)}", log_file)
            subprocess.run(cmd, check=True)
            _log(f"{mode} completed for {sample}. Output: {sample_out}", log_file)

        except subprocess.CalledProcessError as e:
            _log(f"ERROR: {mode} failed for {sample}: {e}", log_file)
            sys.stderr.write(f"ERROR MESSAGE: {mode} failed for {sample}: {e}\n")
            continue

def QC_10x(output_folder, n_genes_by_counts_thres, pct_counts_mt_thres, total_counts_thres, log_file):
    """
        Perform comprehensive QC on 10x Genomics data:
        - Load raw & filtered matrices
        - Compute QC metrics
        - Visualize QC (violin & bar plots)
        - Filter poor-quality cells
        - Detect potential doublets using Scrublet
    """

    files_to_gzip = [
        ("matrix.mtx", "matrix.mtx.gz"),
        ("barcodes.tsv", "barcodes.tsv.gz"),
        ("features.tsv", "features.tsv.gz")
    ]

    checking_directory = os.path.join(output_folder, "STARsolo_out")
    for sample in os.listdir(checking_directory):
        sample_path = os.path.join(checking_directory, sample)
        raw_dir = os.path.join(sample_path, "Solo.out", "Gene", "raw")
        filtered_dir = os.path.join(sample_path, "Solo.out", "Gene", "filtered")

        if not os.path.isdir(filtered_dir):
            _log(f"Skipping {sample}: no filtered directory found.", log_file)
            continue

        for in_name, out_name in files_to_gzip:
            in_path = os.path.join(filtered_dir, in_name)
            out_path = os.path.join(filtered_dir, out_name)

            # If gz file exists already, skip
            if os.path.isfile(out_path):
                continue

            # If ungzipped file exists, compress it
            if os.path.isfile(in_path):
                _log(f"Gzipping {in_name} for {sample}", log_file)
                try:
                    with open(in_path, "rb") as f_in, gzip.open(out_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                except Exception as e:
                    _log(f"ERROR gzipping {in_name} for {sample}: {e}", log_file)
                    continue
            else:
                _log(f"ERROR: {in_name} and {out_name} both missing for {sample}", log_file)
                continue

        # --- Load data ---
        try:
            _log(f"Loading filtered 10x data for {sample}", log_file)
            adata = sc.read_10x_mtx(filtered_dir, var_names="gene_symbols", cache=True, make_unique=True)
        except Exception as e:
            _log(f"ERROR loading data for {sample}: {e}", log_file)
            continue

        adata.layers["counts"] = adata.X.copy()

        # --- Basic QC metrics ---
        _log(f"Calculating QC metrics for {sample}", log_file)
        adata.var["mt"] = adata.var_names.str.startswith("MT-")  # Mitochondrial genes
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

        # --- Violin plots ---
        qc_plot_dir = os.path.join(sample_path, "QC_plots")
        os.makedirs(qc_plot_dir, exist_ok=True)
        _log(f"Generating QC violin plots for {sample}", log_file)

        sc.settings.figdir = qc_plot_dir
        sc.pl.violin(
            adata,
            ["total_counts", "n_genes_by_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            save=f"_{sample}_qc_violin.png"
        )

        # --- Barplot summary ---
        plt.figure(figsize=(6, 4))
        adata.obs["pct_counts_mt"].hist(bins=50)
        plt.title(f"{sample}: Mitochondrial % per Cell")
        plt.xlabel("% mitochondrial counts")
        plt.ylabel("Cell count")
        plt.tight_layout()
        barplot_path = os.path.join(qc_plot_dir, f"{sample}_mt_barplot.png")
        plt.savefig(barplot_path)
        plt.close()
        _log(f"Saved barplot: {barplot_path}", log_file)

        # --- Filtering poor-quality cells ---
        _log(f"Filtering low-quality cells for {sample}", log_file)
        before = adata.n_obs
        adata = adata[
            (adata.obs["n_genes_by_counts"] > n_genes_by_counts_thres)
            & (adata.obs["pct_counts_mt"] < pct_counts_mt_thres)
            & (adata.obs["total_counts"] > total_counts_thres)
        ].copy()
        after = adata.n_obs
        _log(f"{sample}: filtered {before - after} cells; {after} remain.", log_file)
        # handle empty Annadata after filtering
        if adata.n_obs == 0:
            _log(f"{sample}: all cells filtered out, skipping save.", log_file)
            continue
        # --- Doublet detection using Scrublet ---
        _log(f"Running Scrublet for doublet detection on {sample}", log_file)
        if adata.n_obs <= 50:
            _log(
                f"Skipping Scrublet for {sample}: too few cells after filtering (n={adata.n_obs})",
                log_file,
            )

        else:

            try:
                scrub = scr.Scrublet(adata.X)
                doublet_scores, predicted_doublets = scrub.scrub_doublets()

                adata.obs["doublet_score"] = doublet_scores
                adata.obs["predicted_doublet"] = predicted_doublets

                scrub.plot_histogram()
                scrub_hist = os.path.join(qc_plot_dir, f"{sample}_scrublet_hist.png")
                plt.savefig(scrub_hist)
                plt.close()
                _log(f"Scrublet histogram saved to {scrub_hist}", log_file)

                doublet_rate = (predicted_doublets.sum() / len(predicted_doublets)) * 100
                _log(f"{sample}: estimated doublet rate = {doublet_rate:.2f}%", log_file)
            except Exception as e:
                _log(f"WARNING: Scrublet failed for {sample}: {e}", log_file)

        # --- Save QC results ---
        qc_out = os.path.join(sample_path, "QC_filtered.h5ad")
        adata.write(qc_out)
        _log(f"QC complete for {sample}. Saved to {qc_out}", log_file)


def QC_smartseq2(output_folder, log_file, platform="smartseq2"):
    """
        Perform QC and visualization for Smart-seq2 or Smart-seq3 data.
    """

    _log(f"Running QC for {platform} samples...", log_file)

    for sample in os.listdir(output_folder):
        sample_path = os.path.join(output_folder, sample)
        count_file = os.path.join(sample_path, "counts.tsv")

        # --- Detect matrix type ---
        if os.path.exists(os.path.join(sample_path, "matrix.mtx")):
            adata = sc.read_10x_mtx(sample_path, var_names="gene_symbols", cache=True, make_unique=True)
        elif os.path.exists(count_file):
            adata = sc.read_text(count_file).T
        else:
            _log(f"No count matrix found for {sample}. Skipping.", log_file)
            continue

        adata.layers["counts"] = adata.X.copy()

        # --- Calculate QC metrics ---
        sc.pp.calculate_qc_metrics(adata, qc_vars=["MT"], percent_top=None, log1p=False, inplace=True)

        qc_dir = os.path.join(sample_path, "QC")
        os.makedirs(qc_dir, exist_ok=True)

        # --- Violin plots ---
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_MT"],
            jitter=0.4,
            multi_panel=True,
            save=f"_{sample}.png",
            show=False
        )

        # --- Bar plot: total counts ---
        plt.figure()
        plt.bar(range(len(adata.obs_names)), adata.obs["total_counts"])
        plt.title(f"{sample} - Total Counts per Cell")
        plt.xlabel("Cells")
        plt.ylabel("Total Counts")
        plt.tight_layout()
        plt.savefig(os.path.join(qc_dir, f"{sample}_barplot_total_counts.png"))
        plt.close()

        # --- Filtering ---
        adata = adata[adata.obs.n_genes_by_counts < 10000, :]
        adata = adata[adata.obs.pct_counts_MT < 15, :]

        _log(f"Filtered {sample}: retained {adata.n_obs} cells.", log_file)

        # Save results
        adata.write_h5ad(os.path.join(qc_dir, f"{sample}_filtered_qc.h5ad"))
        _log(f"QC results saved for {sample}", log_file)
        
        
def NormalizeAllSamples(output_dir, HVG_selection, number_top_genes, log_file):
    """
        Normalize all QC-filtered .h5ad files from each sample directory.

        For each sample in the output_dir:
            - Loads its 'QC_filtered.h5ad' file
            - Performs total-count normalization (counts per 10,000)
            - Applies log1p transformation
            - Saves the normalized file as 'normalized.h5ad' in the same folder

        Args:
            output_dir (str): Path to the main output directory containing per-sample subfolders.
            log_file (str): Path to the central log file.

        Error Codes:
            ERROR CODE 13: Normalization failed.
    """

    output_dir = os.path.join(output_dir, "STARsolo_out")
    _log(f"Starting normalization of all samples in {output_dir}", log_file)

    try:
        # Loop through all sample folders
        for sample in os.listdir(output_dir):
            sample_path = os.path.join(output_dir, sample)
            if not os.path.isdir(sample_path):
                continue

            input_h5ad = os.path.join(sample_path, "QC_filtered.h5ad")
            output_h5ad = os.path.join(sample_path, "normalized.h5ad")
            output_tsv = os.path.join(sample_path, "normalized.tsv")
            output_hvg_h5ad = os.path.join(sample_path, "normalized_hvg.h5ad")
            output_hvg_tsv = os.path.join(sample_path, "normalized_hvg.tsv")

            if not os.path.exists(input_h5ad):
                _log(f"Skipping {sample}: No QC_filtered.h5ad found.", log_file)
                continue

            _log(f"Normalizing sample: {sample}", log_file)
            adata = sc.read_h5ad(input_h5ad)

            if HVG_selection:
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=number_top_genes,
                    flavor='seurat_v3'
                )

                adata_hvg = adata[:, adata.var['highly_variable']].copy()
                adata = adata[:, adata.var['highly_variable']]
                adata_hvg.write(output_hvg_h5ad)
                _log(f"Saved normalized HVG data for {sample} to {output_hvg_h5ad}", log_file)
                
                norm_hvg_count_matix = pd.DataFrame(
                    adata_hvg.X.toarray() if hasattr(adata_hvg.X, "toarray") else adata_hvg.X,
                    index = adata_hvg.obs_names,
                    columns = adata_hvg.var_names
                )
                norm_hvg_count_matix.to_csv(output_hvg_tsv, sep = '\t')
                _log(f"Saved normalized data for {sample} to {output_hvg_tsv}", log_file)
                
            # Normalization steps
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            sc.pp.filter_genes(adata, min_variance=1e-8)
            sc.pp.scale(adata, max_value=10)

            # Save normalized file
            adata.write(output_h5ad)
            _log(f"Saved normalized data for {sample} to {output_h5ad}", log_file)

            norm_count_matix = pd.DataFrame(
                adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
                index = adata.obs_names,
                columns = adata.var_names
            )
            norm_count_matix.to_csv(output_tsv, sep = '\t')
            _log(f"Saved normalized data for {sample} to {output_tsv}", log_file)

        _log("All sample normalizations completed successfully.", log_file)

    except Exception as e:
        msg = f"ERROR: Normalization failed due to: {e}"
        _log(msg, log_file)
        sys.stderr.write(msg + "\n")
        sys.exit(13)
        
        
def main():
    """
        Main function for the script.
        
        ERROR CODE 4: Invalid or missing 'input_dir' in configuration file
        ERROR CODE 9:  Missing 'Fastq_file_format' in configuration file.
        ERROR CODE 10:  Missing required genome preparation parameters in configuration file
        ERROR CODE 12: Platform '{platform}' not recognized. Use: SmartSeq, 10x, or DropSeq. Exiting
    """

    # TODO: Implement chmod 777 for the input folder to avoid any kind of file permission errors

    config_folder = os.path.dirname(os.path.abspath(__file__))
    logfile = get_log_file(config_folder)
    _log("Log file created", logfile)
    
    configuration_file_path = CheckingConfiguration(config_folder)
    _log("Configuration file is good, ready to start", logfile)
    
    configuration = LoadingConfiguration(configuration_file_path, logfile)
    _log("Configuration file successfully loaded.", logfile)
    output_dir = configuration.get("output_dir")

    # Run backup
    result_path = backup_config(config_folder, output_dir)

    print(f"Backup complete! Files saved to: {result_path}")

    input_dir = configuration.get("input_dir")
    Fastq_file_format = configuration.get("Fastq_file_format")  # must be "merged" or "subdir"

    if not input_dir or not os.path.isdir(input_dir):
        sys.stderr.write("ERROR: Invalid or missing 'input_dir' in configuration file.\n")
        _log("ERROR: Invalid or missing 'input_dir' in configuration file.", logfile)
        sys.exit(4)

    if not Fastq_file_format:
        sys.stderr.write("ERROR: Missing 'Fastq_file_format' in configuration file.\n")
        _log("ERROR: Missing 'Fastq_file_format' in configuration file.", logfile)
        sys.exit(9)

    _log(f"Checking FASTQ files in: {input_dir} (mode: {Fastq_file_format})", logfile)
    CheckFASTQFiles(input_dir, Fastq_file_format, logfile)
    _log("FASTQ file structure successfully validated.", logfile)

    genome_dir = configuration.get("genome_dir")
    fasta_file = configuration.get("fasta_file")
    genome_index_params = configuration.get("genome_index_params", {})
    splice_junction_params = configuration.get("splice_junction_params", {})
    HVG_selection = configuration.get("HVG_selection")
    number_top_genes = configuration.get("number_top_genes")
    n_genes_by_counts_thres = configuration.get("n_genes_by_counts_thres")
    pct_counts_mt_thres = configuration.get("pct_counts_mt_thres")
    total_counts_thres = configuration.get("total_counts_thres")
    
    # Validate required genome parameters
    if not all([genome_dir, fasta_file, genome_index_params, splice_junction_params]):
        sys.stderr.write("ERROR: Missing required genome preparation parameters in configuration file.\n")
        _log("ERROR: Missing required genome preparation parameters in configuration file.", logfile)
        sys.exit(10)

    # Prepare STAR genome index using new function signature
    PrepareSTARGenome(
        genome_dir=genome_dir,
        fasta_file=fasta_file,
        genome_index_params=genome_index_params,
        splice_junction_params=splice_junction_params,
        log_file=logfile
    )
    platform = configuration.get("platform", "").strip()

    if platform in ("10x", "DropSeq"):
        # STARsolo for both
        RunSTARUnified(configuration, logfile)
        QC_10x(output_dir, n_genes_by_counts_thres, pct_counts_mt_thres, total_counts_thres, logfile)

    elif platform == "SmartSeq":
        # Standard STAR
        RunSTARUnified(configuration, logfile)
        QC_smartseq2(output_dir, logfile)

    else:
        msg = f"Platform '{platform}' not recognized. Use: SmartSeq, 10x, or DropSeq."
        _log(msg, logfile)
        sys.stderr.write(msg + "\n")
        sys.exit(12)
        
    # ✅ Normalization step runs for both platforms
    NormalizeAllSamples(output_dir, HVG_selection, number_top_genes, logfile)
    _log("Normalization completed for all samples.", logfile)
    

if __name__ == '__main__':
    main()
