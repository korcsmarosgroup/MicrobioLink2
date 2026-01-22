import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import yaml
import sys
import os
from time import strftime
from scipy import sparse

############################################################
# LOGGING
############################################################
    """
    Generate a timestamped log file path inside the provided configuration folder.

    Parameters
    ----------
    config_folder : str
        Path to the folder containing the analysis configuration.

    Returns
    -------
    str
        Full file path to a new log file, named as analysis_<timestamp>.log.
    """
def get_log_file(config_folder):
    ts = strftime("%Y-%m-%d_%H-%M-%S")
    return f"{config_folder}/analysis_{ts}.log"

def _log(message, log_file):
    with open(log_file, "a") as f:
        f.write(f"[{strftime('%Y-%m-%d %H:%M:%S')}] {message}\n")

############################################################
# CSV → AnnData Loader
############################################################
    """
    Load a tab-separated normalized expression matrix (CSV/TSV) into an AnnData object.

    Parameters
    ----------
    path : str
        Path to the TSV file containing gene expression data.
    sample : str
        Sample name assigned to all observations in the AnnData object.
    log_file : str
        Log file path for logging status messages.

    Returns
    -------
    AnnData
        AnnData object containing expression matrix, genes in `.var_names`,
        cells in `.obs_names`, and sample annotation in `adata.obs["sample"]`.

    Notes
    -----
    - Expects tab-separated file with genes as columns and cells as rows.
    - Logs success/failure into log_file.
    """
def load_csv_to_adata(path, sample, log_file):
    _log(f"Loading CSV file {path}", log_file)

    df = pd.read_csv(path, sep="\t", index_col=0)

    # Create AnnData
    adata = ad.AnnData(sparse.csr_matrix(df.values))
    adata.obs_names = df.index
    adata.var_names = df.columns

    # Add sample name
    adata.obs["sample"] = sample

    _log(f"Successfully loaded CSV for sample {sample}", log_file)
    return adata

############################################################
# Configuration loading
############################################################

def CheckingConfiguration(config_folder):
    config_path = f"{config_folder}/analysis_config.yaml"
    if not os.path.isfile(config_path):
        sys.stderr.write(f"ERROR MESSAGE: Missing analysis_config.yaml at {config_path}\n")
        sys.exit(1)
    return config_path

def LoadingConfiguration(config_path, log_file):
    try:
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
    except Exception as e:
        _log(f"ERROR reading configuration: {e}", log_file)
        sys.stderr.write(f"ERROR MESSAGE: Cannot read configuration: {e}\n")
        sys.exit(2)
    return config

############################################################
# LOAD NORMALIZED OUTPUTS
############################################################

def load_all_samples(input_dir, log_file):
    """
    Loads:
      - normalized.h5ad   (preferred)
      - normalized.tsv    (fallback)
      - normalized_hvg.*  (if PCA/Harmony set to use HVG)
    """
    if not os.path.isdir(input_dir):
        _log(f"ERROR: input_dir {input_dir} does not exist", log_file)
        sys.exit(4)

    samples = []

    for sample in os.listdir(input_dir):
        spath = os.path.join(input_dir, sample)
        if not os.path.isdir(spath):
            continue

        # Preferred input file
        h5 = os.path.join(spath, "normalized.h5ad")
        tsv = os.path.join(spath, "normalized.tsv")

        if os.path.isfile(h5):
            _log(f"Loading H5AD for sample {sample}", log_file)
            adata = sc.read_h5ad(h5)
            adata.obs["sample"] = sample
            samples.append(adata)
            continue

        if os.path.isfile(tsv):
            adata = load_csv_to_adata(tsv, sample, log_file)
            samples.append(adata)
            continue

        _log(f"WARNING: Sample {sample} contains no usable normalized output", log_file)

    if len(samples) == 0:
        sys.stderr.write("ERROR MESSAGE: No samples found for analysis.\n")
        sys.exit(5)

    return samples

############################################################
# PCA
############################################################

def run_pca(adata, n_comps, log_file):
    _log(f"Running PCA with {n_comps} components", log_file)
    sc.pp.scale(adata, max_value=10)          # Ensure scaling (important)
    sc.tl.pca(adata, n_comps=n_comps)
    return adata

############################################################
# Harmony
############################################################

def batch_harmony(adata, batch_key, log_file):
    try:
        import harmonypy as hm
    except ImportError:
        sys.stderr.write("ERROR MESSAGE: harmonypy is required for Harmony.\n")
        sys.exit(6)

    _log(f"Running Harmony batch correction using key '{batch_key}'", log_file)
    harmony_out = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_key)
    adata.obsm["X_harmony"] = harmony_out.Z_corr.T
    return adata

############################################################
#UMAP + tSNE
############################################################

def run_umap(adata, log_file):
    _log("Running UMAP", log_file)
    sc.tl.umap(adata)
    return adata

def run_tsne(adata, log_file):
    _log("Running t-SNE", log_file)
    sc.tl.tsne(adata)
    return adata

############################################################
# Clustering
############################################################

def run_clustering(adata, method, resolution, log_file):
    _log(f"Clustering using {method} at resolution {resolution}", log_file)
    if method == "leiden":
        sc.tl.leiden(adata, resolution=resolution)
    elif method == "louvain":
        sc.tl.louvain(adata, resolution=resolution)
    return adata

############################################################
# Marker genes
############################################################
    """
    Run rank_genes_groups to compute differential expression markers per cluster.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing clustering labels.
    groupby : str
        Column in `adata.obs` defining groups for comparison (e.g., "leiden").
    log_file : str
        Log file path.

    Returns
    -------
    AnnData
        AnnData with marker gene results stored in `adata.uns["rank_genes_groups"]`.

    Notes
    -----
    Uses Wilcoxon rank-sum test by default.
    """
def run_marker_genes(adata, groupby, log_file):
    _log(f"Running marker gene analysis by {groupby}", log_file)
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon")
    return adata

############################################################
# Main
############################################################
    """
    Run the complete single-cell analysis workflow:

    Steps performed
    ---------------
    1. Load configuration (YAML)
    2. Initialize logging
    3. Load normalized sample data (H5AD or TSV)
    4. Merge samples into a combined AnnData
    5. Run PCA
    6. Optionally run Harmony batch correction
    7. Compute neighbors
    8. Compute embeddings (UMAP, t-SNE)
    9. Perform clustering (Leiden or Louvain)
    10. Run marker gene analysis
    11. Save output AnnData to disk
    """

def main():
    config_folder = os.path.dirname(os.path.abspath(__file__))
    log_file = get_log_file(config_folder)
    _log("Log file created", log_file)

    cfg_path = CheckingConfiguration(config_folder)
    config = LoadingConfiguration(cfg_path, log_file)

    input_dir = config.get("input_dir")
    output_file = config.get("output_adata", "analysis_output.h5ad")

    # Load normalized data (CSV + H5AD supported)
    samples = load_all_samples(input_dir, log_file)

    # Merge into one AnnData
    _log("Merging samples", log_file)
    adata = ad.concat(samples, join="outer", label="sample")

    # PCA
    n_comps = config.get("pca_components", 50)
    adata = run_pca(adata, n_comps, log_file)

    # Batch correction
    if config.get("batch_correction", {}).get("use_harmony", False):
        batch_key = config["batch_correction"]["batch_key"]
        adata = batch_harmony(adata, batch_key, log_file)
        rep = "X_harmony"
    else:
        rep = "X_pca"

    # embeddings

    if config.get("embedding", {}).get("umap", True):
        run_umap(adata, log_file)
    if config.get("embedding", {}).get("tsne", False):
        run_tsne(adata, log_file)

    # Clustering
    cl = config.get("clustering", {})
    adata = run_clustering(
        adata,
        cl.get("method", "leiden"),
        cl.get("resolution", 1.0),
        log_file
    )

    # Marker genes
    if config.get("marker_genes", {}).get("run", True):
        groupby = config["marker_genes"].get("groupby", "leiden")
        run_marker_genes(adata, groupby, log_file)

    adata.write(output_file)
    _log(f"Analysis complete. Saved {output_file}", log_file)


if __name__ == "__main__":
    main()