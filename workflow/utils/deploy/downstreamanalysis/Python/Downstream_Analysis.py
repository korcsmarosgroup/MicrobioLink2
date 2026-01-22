import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import yaml
import sys
import os
from time import strftime
from scipy import sparse
import matplotlib.pyplot as plt
import shutil



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

############################################################
# LOAD NORMALIZED OUTPUTS
############################################################

def load_all_samples(input_dir, log_file):
    """
    Loads:
      - any .h5ad file (preferred)
      - any .tsv file (fallback)
    """
    if not os.path.isdir(input_dir):
        _log(f"ERROR: input_dir {input_dir} does not exist", log_file)
        sys.exit(4)

    samples = []

    for sample in os.listdir(input_dir):
        spath = os.path.join(input_dir, sample)
        if not os.path.isdir(spath):
            continue

        # Find files by extension (not name)
        h5ad_files = sorted(
            f for f in os.listdir(spath) if f.lower().endswith(".h5ad")
        )
        tsv_files = sorted(
            f for f in os.listdir(spath) if f.lower().endswith(".tsv")
        )

        if h5ad_files:
            h5_path = os.path.join(spath, h5ad_files[0])
            _log(f"Loading H5AD for {sample}: {h5ad_files[0]}", log_file)
            adata = sc.read_h5ad(h5_path)
            adata.obs["sample"] = sample
            samples.append(adata)
            continue

        if tsv_files:
            tsv_path = os.path.join(spath, tsv_files[0])
            _log(f"Loading TSV for sample {sample}: {tsv_files[0]}", log_file)
            adata = load_csv_to_adata(tsv_path, sample, log_file)
            samples.append(adata)
            continue

        _log(f"WARNING: Sample {sample} contains no usable normalized output", log_file)

    if len(samples) == 0:
        sys.stderr.write("ERROR MESSAGE: No samples found for analysis.\n")
        sys.exit(5)

    return samples
############################################################
# SCVI
############################################################

def batch_scvi(adata, batch_key, n_latent, max_epochs, log_file):
    try:
        import scvi
    except ImportError:
        sys.stderr.write("ERROR MESSAGE: scvi-tools is required for scVI.\n")
        sys.exit(6)

    _log(f"Running scVI with batch key '{batch_key}', n_latent={n_latent}, max_epochs={max_epochs}", log_file)

    # Setup AnnData for SCVI
    if "counts" in adata.layers:
        scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
    else:
        _log("WARNING: No 'counts' layer found; using adata.X as counts (not ideal).", log_file)
        scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)

    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(max_epochs=max_epochs)
    
    # Latent space for clustering / UMAP
    adata.obsm["X_scVI"] = model.get_latent_representation()

    # scVI-normalized expression
    adata.layers["scvi_norm"] = model.get_normalized_expression(
        library_size=1e4
        )
    
    adata.layers["scvi_log"] = np.log1p(adata.layers["scvi_norm"])
   
    return adata
############################################################
# PCA
############################################################

def run_pca(adata, n_comps, log_file):
    _log(f"Running PCA with {n_comps} components", log_file)
    sc.pp.scale(adata, max_value=10)          # Ensure scaling
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

    # ---- FIX: ensure contiguous array (no negative strides) ----
    X_pca = adata.obsm["X_pca"].copy()

    harmony_out = hm.run_harmony(X_pca, adata.obs, batch_key)
    adata.obsm["X_harmony"] = harmony_out.Z_corr

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
# UMAP plotting utilities
############################################################

def compute_neighbors_and_umap(adata, rep, log_file):
    _log(f"Computing neighbors using {rep}", log_file)
    sc.pp.neighbors(adata, use_rep=rep)
    sc.tl.umap(adata)
    return adata

def save_umap_plot(adata, color, output_dir, filename, log_file):
    os.makedirs(output_dir, exist_ok=True)

    _log(f"Saving UMAP plot: {filename}", log_file)

    sc.settings.figdir = output_dir

    sc.pl.umap(
        adata,
        color=color,
        show=False,
        save=f"_{filename}"
    )
############################################################
# Clustering
############################################################

def run_clustering(adata, method, resolution, rep, log_file):
    _log(f"Clustering using {method} at resolution {resolution}", log_file)

    # Build neighbors graph using the chosen latent space
    sc.pp.neighbors(
        adata,
        use_rep=rep,
        n_neighbors=15,
        n_pcs=None
    )

    if method == "leiden":
        sc.tl.leiden(
            adata,
            resolution=resolution,
            flavor="igraph",
            directed=False,
            n_iterations=2
        )
    elif method == "louvain":
        sc.tl.louvain(
            adata,
            resolution=resolution
        )

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
def run_marker_genes(adata, groupby, log_file, batch_method=None):
    _log(f"Running marker gene analysis by {groupby}", log_file)

    if batch_method == "scvi":
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method="wilcoxon",
            layer="scvi_log",
            use_raw=False
        )
    else:
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method="wilcoxon",
            use_raw=True  # default for Harmony/PCA
        )

    return adata

############################################################
# Cell annotation
############################################################

def annotate_cells(adata, config, log_file, batch_method="pca"):
    """
    Annotate cells depending on batch correction:
        - Harmony or PCA -> Scanpy ingest
        - scVI -> scANVI
    Fills adata.obs['celltype'] and optionally adata.obs['prediction_confidence'].

    Parameters
    ----------
    adata : AnnData
        Query dataset to annotate.
    config : dict
        Configuration dictionary containing:
            - cell_annotation_file : path to reference .h5ad or marker CSV
            - cell_annotation_column : name of column in reference with cell types
            - batch_correction : dict (contains batch_key if needed)
            - cell_annotation : dict (max_epochs if scANVI)
    log_file : str
        Path to log file.
    batch_method : str
        One of 'pca', 'harmony', 'scvi'.
    """

    annotation_file = config.get("cell_annotation_file", None)
    celltype_col = config.get("cell_annotation_column", "celltype")
    batch_col = config.get("batch_correction", {}).get("batch_key", None)
    max_epochs = config.get("cell_annotation", {}).get("max_epochs", 50)

    # ------------------ No annotation file ------------------
    if not annotation_file:
        _log("No annotation file provided. Using Leiden clusters as cell types.", log_file)
        adata.obs[celltype_col] = adata.obs.get(
            config.get("marker_genes", {}).get("groupby", "leiden")
        )
        return adata

    # ------------------ Reference h5ad ------------------
    if annotation_file.endswith(".h5ad") and os.path.isfile(annotation_file):
        ref = sc.read_h5ad(annotation_file)
        if celltype_col not in ref.obs.columns:
            _log(f"Reference file missing '{celltype_col}'. Using Leiden clusters as fallback.", log_file)
            adata.obs[celltype_col] = adata.obs.get(
                config.get("marker_genes", {}).get("groupby", "leiden")
            )
            return adata

        # Keep only shared genes
        shared_genes = ref.var_names.intersection(adata.var_names)
        ref = ref[:, shared_genes]
        adata = adata[:, shared_genes]

        if batch_method in ["pca", "harmony"]:
            # ------------------ Use Scanpy ingest ------------------
            _log(f"Using Scanpy ingest for label transfer ({batch_method} batch correction).", log_file)
            rep = "X_harmony" if "X_harmony" in adata.obsm else "X_pca"
            sc.tl.ingest(adata, ref, obs=celltype_col, embedding_method=rep)
            _log("Label transfer complete with ingest.", log_file)
            return adata

        elif batch_method == "scvi":
            # ------------------ Use scANVI ------------------
            _log("Using scANVI for label transfer (scVI batch correction).", log_file)
            import scvi

            # If reference has only one batch, batch_key can be None
            if batch_col not in ref.obs.columns:
                batch_col = None
                _log("Reference has no batch column; scANVI will run without batch info.", log_file)

            # Setup scANVI on reference
            scvi.model.SCANVI.setup_anndata(ref, labels_key=celltype_col, batch_key=batch_col)

            # Load query into scANVI
            scanvi_model = scvi.model.SCANVI.load_query_data(ref, adata, labels_key=celltype_col)

            # Train scANVI
            _log(f"Training scANVI for {max_epochs} epochs.", log_file)
            scanvi_model.train(max_epochs=max_epochs, plan_kwargs={"weight_decay": 0.0})

            # Predict labels
            pred_labels, pred_probs = scanvi_model.predict(adata, return_proba=True)
            adata.obs[celltype_col] = pred_labels
            adata.obs["prediction_confidence"] = pred_probs.max(axis=1)

            _log("Label transfer complete with scANVI.", log_file)
            return adata

        else:
            _log(f"Unknown batch_method '{batch_method}'. Using Leiden clusters as fallback.", log_file)
            adata.obs[celltype_col] = adata.obs.get(
                config.get("marker_genes", {}).get("groupby", "leiden")
            )
            return adata

    # ------------------ Marker gene list fallback ------------------
    if os.path.isfile(annotation_file):
        _log(f"Using marker gene list for annotation: {annotation_file}", log_file)
        with open(annotation_file) as f:
            markers = [line.strip() for line in f if line.strip()]

        celltypes = []
        for cell in adata.obs_names:
            expr = (
                adata[cell, markers].X.toarray().flatten()
                if sparse.issparse(adata.X)
                else adata[cell, markers].X.flatten()
            )
            if np.any(expr > 0):
                celltypes.append(markers[np.argmax(expr)])
            else:
                celltypes.append("Unknown")
        adata.obs[celltype_col] = celltypes
        _log("Cell annotation using marker list complete.", log_file)
        return adata

    # ------------------ Fallback to Leiden ------------------
    _log("Annotation file invalid. Using Leiden clusters as fallback.", log_file)
    adata.obs[celltype_col] = adata.obs.get(
        config.get("marker_genes", {}).get("groupby", "leiden")
    )
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
    11. cell annotations
    11. Save output AnnData to disk
    """

def main():

    config_folder = os.path.dirname(os.path.abspath(__file__))
    log_file = get_log_file(config_folder)
    _log("Log file created", log_file)

    cfg_path = CheckingConfiguration(config_folder)
    config = LoadingConfiguration(cfg_path, log_file)
    output_dir = config.get("output_dir")

    # Run backup
    result_path = backup_config(config_folder, output_dir)

    print(f"Backup complete! Files saved to: {result_path}")
    
    input_dir = config.get("input_dir")
    output_file = config.get("output_adata", "analysis_output.h5ad")

    plot_cfg = config.get("plots", {}).get("umap", {})
    plot_umap = plot_cfg.get("run", False)
    plot_dir = plot_cfg.get("output_dir", "umap_plots")

    # ------------------ Load samples ------------------
    samples = load_all_samples(input_dir, log_file)

    # Merge samples
    _log("Merging samples", log_file)
    adata = ad.concat(samples, join="outer", label="sample")

    # ------------------ Batch correction ------------------
    batch_cfg = config.get("batch_correction", {})
    use_scvi = batch_cfg.get("use_scvi", False)
    use_harmony = batch_cfg.get("use_harmony", False)

    if use_scvi:
        # scVI does NOT use PCA
        if "batch_key" not in batch_cfg:
            sys.stderr.write("ERROR: batch_key must be provided for scVI.\n")
            sys.exit(1)
        n_latent = batch_cfg.get("n_latent", 30)
        max_epochs = batch_cfg.get("max_epochs", 400)
        adata = batch_scvi(adata, batch_cfg["batch_key"], n_latent, max_epochs, log_file)
        rep_for_downstream = "X_scVI"

    elif use_harmony:
        # PCA is required for Harmony
        n_comps = config.get("pca_components", 50)
        adata = run_pca(adata, n_comps, log_file)

        if "batch_key" not in batch_cfg:
            sys.stderr.write("ERROR: batch_key must be provided for Harmony.\n")
            sys.exit(1)
        batch_key = batch_cfg["batch_key"]
        adata = batch_harmony(adata, batch_key, log_file)
        rep_for_downstream = "X_harmony"

    else:
        # No batch correction, optionally run PCA
        n_comps = config.get("pca_components", 50)
        adata = run_pca(adata, n_comps, log_file)
        rep_for_downstream = "X_pca"

    # ------------------ Neighbors + UMAP ------------------
    if plot_umap:
        adata = compute_neighbors_and_umap(adata, rep_for_downstream, log_file)
        save_umap_plot(
            adata,
            color="sample",
            output_dir=plot_dir,
            filename="umap_after_batch.png",
            log_file=log_file
        )

    # ------------------ Clustering ------------------
    cl = config.get("clustering", {})
    adata = run_clustering(
        adata,
        cl.get("method", "leiden"),
        cl.get("resolution", 1.0),
        rep_for_downstream,
        log_file
    )

    if plot_umap:
        save_umap_plot(
            adata,
            color=cl.get("method", "leiden"),
            output_dir=plot_dir,
            filename="umap_after_clustering.png",
            log_file=log_file
        )
        # Determine batch method from config
        if use_scvi:
            batch_method = "scvi"
        elif use_harmony:
            batch_method = "harmony"
        else:
            batch_method = "pca"
                    
    # ------------------ Marker genes ------------------
    if config.get("marker_genes", {}).get("run", True):
        groupby = config["marker_genes"].get("groupby", "leiden")
        run_marker_genes(adata, groupby, log_file, batch_method=batch_method)

        # ------------------ Cell annotation ------------------

            # Annotate cells
        adata = annotate_cells(adata, config, log_file, batch_method=batch_method)
        
        if plot_umap:
            save_umap_plot(
                adata,
                color="celltype",
                output_dir=plot_dir,
                filename="umap_after_cell_annotation.png",
                log_file=log_file
            )

    # ------------------ Save output ------------------
    adata.write(output_file)
    _log(f"Analysis complete. Saved {output_file}", log_file)


if __name__ == "__main__":
    main()
