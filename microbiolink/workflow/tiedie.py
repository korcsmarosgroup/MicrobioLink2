"""TieDie tied-diffusion network propagation: build inputs, run TieDie, format the output network."""

import contextlib
import io
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from ..utils import id_resolution

UPSTREAM_SIGN_DEFAULT = "-"  # DMI/DDI tables carry no sign column (Decisions 2, 5)
HMI_PAIR_COLUMNS = ["human_uniprot_id", "bacterial_uniprot_id"]
STIMULATION_INTERACTION = "stimulates>"
INHIBITION_INTERACTION = "inhibits>"
CAUSAL_SIF_COLUMNS = ["Source.node", "Relationship", "Target.node"]
NETWORK_COLUMNS = ["Target.node", "Source.node", "Relationship", "layer"]
NODE_TABLE_COLUMNS = [
    "node",
    "all_nodes",
    "bacteria_layer",
    "bindingprot_layer",
    "ppi_layer",
    "tf_layer",
    "deg_layer",
]

# Per-layer node membership: (edge layer, node column, annotation column, label).
_NODE_LAYER_SPECS = [
    ("bacteria-bindingprot", "Source.node", "bacteria_layer", "bacteria"),
    ("bacteria-bindingprot", "Target.node", "bindingprot_layer", "bindingprot"),
    ("bindingprot-tf", "Source.node", "ppi1_layer", "bindingprot and/or protein"),
    ("bindingprot-tf", "Target.node", "ppi2_layer", "protein and/or tf"),
    ("tf-deg", "Source.node", "tf_layer", "tf"),
    ("tf-deg", "Target.node", "deg_layer", "deg"),
]


# --- OmniPath network fetch (live; Decision 4) -----------------------------------------------------


def fetch_omnipath_networks() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fetch the OmniPath PPI and CollecTRI TF-target networks (lazy omnipath import).

    Returns:
        A (ppi, tf_target) pair of data frames: the OmniPath signalling PPI and the CollecTRI
        transcription-factor to target-gene network, both with gene-symbol columns.
    """

    import omnipath as op

    ppi = op.interactions.OmniPath.get(genesymbols=1)
    tf_target = op.interactions.Transcriptional.get(
        databases="CollecTRI", genesymbols=1
    )
    return ppi, tf_target


# --- Step 1: build TieDie inputs -------------------------------------------------------------------


def combine_host_microbe_interactions(
    dmi_table: pd.DataFrame,
    ddi_table: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Union the DMI and (optional) DDI tables into distinct host-microbe protein pairs (Decision 3).

    Projects each table to HMI_PAIR_COLUMNS (human_uniprot_id, bacterial_uniprot_id), concatenates,
    and drops duplicates, so a pair present in both tables is kept once (equal-weight union, not
    sum). With ddi_table=None this is just the DMI table's distinct pairs.

    Args:
        dmi_table: A Module 6/7/8 DMI table.
        ddi_table: An optional Module 5 DDI table.

    Returns:
        A data frame of distinct (human_uniprot_id, bacterial_uniprot_id) host-microbe pairs.
    """

    frames = [dmi_table[HMI_PAIR_COLUMNS]]
    if ddi_table is not None:
        frames.append(ddi_table[HMI_PAIR_COLUMNS])

    combined = pd.concat(frames, ignore_index=True)
    return combined.drop_duplicates(ignore_index=True)


def contextualise_networks(
    ppi: pd.DataFrame,
    tf_target: pd.DataFrame,
    expressed_genes: list[str],
    endpoint_genes: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Restrict the PPI to expressed<->expressed edges and the TF-target network to expressed TFs of DEGs.

    Args:
        ppi: The OmniPath PPI network.
        tf_target: The CollecTRI TF-target network.
        expressed_genes: Gene symbols with non-zero expression (the contextualisation universe).
        endpoint_genes: The DEG table (its first column holds the target gene symbols).

    Returns:
        A (contextualised_ppi, contextualised_tf_target) pair. The TF-target table is also the
        step-3 'contextualised regulator-target network'.
    """

    contextualised_ppi = ppi[
        ppi["source_genesymbol"].isin(expressed_genes)
        & ppi["target_genesymbol"].isin(expressed_genes)
    ]

    endpoint_symbols = endpoint_genes.iloc[:, 0]
    contextualised_tf_target = tf_target[
        tf_target["source_genesymbol"].isin(expressed_genes)
        & tf_target["target_genesymbol"].isin(endpoint_symbols)
    ].drop_duplicates()

    return contextualised_ppi, contextualised_tf_target


def build_pathway_sif(contextualised_ppi: pd.DataFrame) -> pd.DataFrame:
    """Build the headerless source/<direction>/target SIF (consensus_stimulation -> stimulates>/inhibits>).

    Args:
        contextualised_ppi: The expressed-gene-contextualised OmniPath PPI.

    Returns:
        A three-column data frame [source, direction, target] of distinct interactions.
    """

    sif = contextualised_ppi[["source", "consensus_stimulation", "target"]].copy()
    # np.where (not .loc assignment) replaces the whole column at once, so a 1/0 int column
    # upcasts cleanly to the string labels; True/False bool values resolve too (True == 1).
    sif["consensus_stimulation"] = np.where(
        sif["consensus_stimulation"] == 1,
        STIMULATION_INTERACTION,
        INHIBITION_INTERACTION,
    )
    return sif.drop_duplicates().rename(columns={"consensus_stimulation": "direction"})


def build_upstream_heats(hmi_pairs: pd.DataFrame) -> pd.DataFrame:
    """Build the upstream heat table from the combined HMI pairs.

    One row per human_uniprot_id, heat = the number of distinct bacterial partners targeting it
    (Decision 5, over the deduplicated DMI union DDI pair set), sign = '-' (Decision 5).

    Args:
        hmi_pairs: Distinct host-microbe pairs from combine_host_microbe_interactions.

    Returns:
        A headerless [protein, heat, sign] data frame.
    """

    partner_counts = hmi_pairs.groupby("human_uniprot_id").size()
    return pd.DataFrame(
        {
            "protein": partner_counts.index,
            "heat": partner_counts.to_numpy(),
            "sign": UPSTREAM_SIGN_DEFAULT,
        }
    )


def build_downstream_heats(
    contextualised_tf_target: pd.DataFrame,
    endpoint_genes: pd.DataFrame,
    value_column: int,
) -> pd.DataFrame:
    """Build the downstream heat table: per TF, the sign-corrected mean DEG log2FC and its sign.

    For each contextualised TF, each DEG target's log2FC is signed by the TF->target consensus
    (stimulation keeps the sign, inhibition flips it) and averaged across the TF's targets
    (Decision 6). Ported from the ground-truth transformation.

    Args:
        contextualised_tf_target: The contextualised TF-target network.
        endpoint_genes: The DEG table (first column = gene symbol).
        value_column: 1-based column of the log2FC value in endpoint_genes.

    Returns:
        A headerless [tf, heat, sign] data frame.
    """

    endpoint = endpoint_genes.rename(
        columns={endpoint_genes.columns[0]: "target_genesymbol"}
    )
    value_name = str(endpoint.columns[value_column - 1])

    merged = contextualised_tf_target.merge(endpoint, on="target_genesymbol")
    merged["exp_sign"] = np.where(
        merged["consensus_stimulation"] == True,  # noqa: E712 - matches ground-truth truthiness
        merged[value_name],
        -merged[value_name],
    )

    grouped = merged.groupby("source")["exp_sign"]
    heats = grouped.mean().reset_index(name="heat")
    heats["sign"] = np.where(heats["heat"] >= 0, "+", "-")
    return heats.rename(columns={"source": "tf"})[["tf", "heat", "sign"]]


# --- Step 2: run TieDie (file-in / folder-out; Decision 1 + 9) -------------------------------------


def _write_tiedie_inputs(
    upstream: pd.DataFrame,
    downstream: pd.DataFrame,
    pathway_sif: pd.DataFrame,
    work_dir: Path,
) -> tuple[Path, Path, Path]:
    """Write the three TieDie inputs (tab-separated, headerless) and return their paths."""

    paths = (
        work_dir / "upstream.input",
        work_dir / "downstream.input",
        work_dir / "pathway.sif",
    )
    for table, path in zip((upstream, downstream, pathway_sif), paths):
        table.to_csv(path, sep="\t", index=False, header=False)
    return paths


def _tiedie_arg_list(
    pathway_path: Path,
    upstream_path: Path,
    downstream_path: Path,
    output_folder: Path,
    size: float,
    alpha: float | None,
    depth: int,
    permute: int,
    use_pagerank: bool,
) -> list[str]:
    """Build the tiedie.cli.main argument list from the input paths and algorithm parameters."""

    args = [
        "-n", str(pathway_path),
        "-u", str(upstream_path),
        "-d", str(downstream_path),
        "--output_folder", str(output_folder),
        "-s", str(size),
        "-c", str(depth),
        "-p", str(permute),
    ]  # fmt: skip
    if alpha is not None:
        args += ["-a", str(alpha)]
    if use_pagerank:
        args.append("--pagerank")
    return args


def _invoke_tiedie(args: list[str]) -> None:
    """Run tiedie.cli.main in-process, re-raising its SystemExit as a RuntimeError with stderr."""

    import tiedie.cli

    captured_stderr = io.StringIO()
    try:
        with contextlib.redirect_stderr(captured_stderr):
            tiedie.cli.main(args)
    except SystemExit as exit_error:
        raise RuntimeError(
            f"TieDie failed (exit {exit_error.code}):\n{captured_stderr.getvalue()}"
        ) from exit_error


def run_tiedie(
    upstream: pd.DataFrame,
    downstream: pd.DataFrame,
    pathway_sif: pd.DataFrame,
    work_dir: Path,
    size: float = 1.0,
    alpha: float | None = None,
    depth: int = 3,
    permute: int = 1000,
    use_pagerank: bool = False,
) -> Path:
    """Write the three inputs into work_dir, run tiedie.cli.main in-process, return the output folder.

    Args:
        upstream: Upstream heat table [protein, heat, sign].
        downstream: Downstream heat table [tf, heat, sign].
        pathway_sif: Pathway SIF [source, direction, target].
        work_dir: Directory to write the three input files and TieDie's output folder into.
        size: TieDie network size-control factor.
        alpha: Optional linker-cutoff override (overrides size when set).
        depth: Causal-path search depth.
        permute: Number of permutations for the significance test.
        use_pagerank: Diffuse with Personalized PageRank instead of the heat kernel.

    Returns:
        The path to TieDie's output folder (work_dir / "TieDIE").

    Raises:
        RuntimeError: If TieDie fails (its internal sys.exit is caught with the captured stderr).
    """

    upstream_path, downstream_path, pathway_path = _write_tiedie_inputs(
        upstream, downstream, pathway_sif, work_dir
    )
    output_folder = work_dir / "TieDIE"
    args = _tiedie_arg_list(
        pathway_path,
        upstream_path,
        downstream_path,
        output_folder,
        size,
        alpha,
        depth,
        permute,
        use_pagerank,
    )
    _invoke_tiedie(args)
    return output_folder


# --- Step 3: format the output network -------------------------------------------------------------


def _hmi_edges(hmi_pairs: pd.DataFrame, causal_sources: pd.Series) -> pd.DataFrame:
    """Build the bacteria-bindingprot edge layer from the combined HMI pairs.

    Keeps only edges whose human target is a causal-network source and whose bacterial protein is
    not itself one (ground-truth process_edges filter). Every HMI edge is signed 'stimulates>'.
    """

    hbps = hmi_pairs.rename(
        columns={
            "bacterial_uniprot_id": "Source.node",
            "human_uniprot_id": "Target.node",
        }
    ).copy()
    hbps["Relationship"] = STIMULATION_INTERACTION
    hbps["layer"] = "bacteria-bindingprot"
    hbps = hbps[NETWORK_COLUMNS]
    return hbps[
        ~hbps["Source.node"].isin(causal_sources)
        & hbps["Target.node"].isin(causal_sources)
    ]


def _tf_deg_edges(
    contextualised_tf_target: pd.DataFrame, causal_targets: pd.Series
) -> pd.DataFrame:
    """Build the tf-deg edge layer from the contextualised TF-target network."""

    tf_deg = contextualised_tf_target.copy()
    tf_deg["consensus_stimulation"] = np.where(
        tf_deg["consensus_stimulation"] == 1,
        STIMULATION_INTERACTION,
        INHIBITION_INTERACTION,
    )
    tf_deg = tf_deg[tf_deg["source"].isin(causal_targets)]
    tf_deg["layer"] = "tf-deg"
    return tf_deg[["source", "target", "consensus_stimulation", "layer"]].rename(
        columns={
            "source": "Source.node",
            "target": "Target.node",
            "consensus_stimulation": "Relationship",
        }
    )


def assemble_network(
    causal_sif: pd.DataFrame,
    hmi_pairs: pd.DataFrame,
    contextualised_tf_target: pd.DataFrame,
) -> pd.DataFrame:
    """Combine the three layers (bacteria-bindingprot, bindingprot-tf, tf-deg) into one edge table.

    The bacteria-bindingprot layer is built from the combined HMI pairs (DMI union DDI), so
    DDI-derived host-microbe edges appear alongside DMI ones (Decision 3). Ported from the
    ground-truth process_edges.

    Args:
        causal_sif: TieDie's causal subnetwork (tiedie.cn.sif) as [Source.node, Relationship, Target.node].
        hmi_pairs: Distinct host-microbe pairs from combine_host_microbe_interactions.
        contextualised_tf_target: The contextualised TF-target network.

    Returns:
        The final network edge table with columns NETWORK_COLUMNS.
    """

    rec_tf = causal_sif.copy()
    rec_tf["layer"] = "bindingprot-tf"
    rec_tf = rec_tf[NETWORK_COLUMNS]

    hbps = _hmi_edges(hmi_pairs, rec_tf["Source.node"])
    tf_deg = _tf_deg_edges(contextualised_tf_target, rec_tf["Target.node"])

    return pd.concat([hbps, rec_tf, tf_deg])


def _node_layer_frames(whole_net: pd.DataFrame) -> list[tuple[str, pd.DataFrame]]:
    """Build the per-layer node-membership frames, each a distinct [node, <annotation>] table."""

    frames = []
    for layer, node_column, annotation_column, label in _NODE_LAYER_SPECS:
        members = whole_net.loc[whole_net["layer"] == layer, [node_column]].copy()
        members = members.rename(columns={node_column: "node"})
        members[annotation_column] = label
        frames.append((annotation_column, members.drop_duplicates()))
    return frames


def _annotate_nodes(
    nodes: pd.DataFrame, layer_frames: list[tuple[str, pd.DataFrame]]
) -> pd.DataFrame:
    """Left-join every layer-membership frame onto a node list and collapse the two PPI columns."""

    annotated = nodes.copy()
    for _, frame in layer_frames:
        annotated = annotated.merge(frame, how="left", on="node")

    # Collapse the two PPI-role columns into one, keeping only the roles a node actually holds.
    # A node in neither role gets "NA" (the same absent sentinel every other layer column uses),
    # so all_nodes drops it instead of carrying the ground truth's "nan,nan" noise.
    annotated["ppi_layer"] = annotated[["ppi1_layer", "ppi2_layer"]].apply(
        lambda roles: ",".join(role for role in roles if pd.notna(role)) or "NA", axis=1
    )
    annotated = annotated.drop(columns=["ppi1_layer", "ppi2_layer"])
    return annotated.fillna("NA")


def _node_membership_table(whole_net: pd.DataFrame) -> pd.DataFrame:
    """Build the per-node layer-membership table (node + all_nodes + one column per layer)."""

    layer_frames = _node_layer_frames(whole_net)

    source_nodes = whole_net[["Source.node"]].rename(columns={"Source.node": "node"})
    target_nodes = whole_net[["Target.node"]].rename(columns={"Target.node": "node"})
    merged = pd.concat(
        [
            _annotate_nodes(source_nodes, layer_frames),
            _annotate_nodes(target_nodes, layer_frames),
        ],
        ignore_index=True,
    )

    layer_columns = [
        "bacteria_layer",
        "bindingprot_layer",
        "ppi_layer",
        "tf_layer",
        "deg_layer",
    ]
    merged["all_nodes"] = merged[layer_columns].apply(
        lambda values: ",".join(value for value in values if value != "NA"), axis=1
    )
    return merged[NODE_TABLE_COLUMNS].drop_duplicates()


def assemble_node_table(
    whole_net: pd.DataFrame,
    heats: pd.DataFrame,
    endpoint_genes: pd.DataFrame,
    value_column: int,
) -> pd.DataFrame:
    """Build the per-node annotation table: layer membership, gene symbol, linker heat, and DEG log2FC.

    Ported from the ground-truth process_node_table / fetch_gene_symbols, re-pointed at
    id_resolution.translate_uniprot_to_gene_symbols.

    Args:
        whole_net: The final network edge table from assemble_network.
        heats: TieDie's linker heats as a [node, LinkerHeats] data frame.
        endpoint_genes: The DEG table (first column = gene symbol).
        value_column: 1-based column of the log2FC value in endpoint_genes.

    Returns:
        The node annotation table (one row per node).
    """

    node_table = _node_membership_table(whole_net)

    symbols = id_resolution.translate_uniprot_to_gene_symbols(
        node_table["node"].tolist()
    )
    node_table["gene_symbol"] = (
        node_table["node"].map(symbols).fillna(node_table["node"])
    )

    node_table = node_table.merge(heats, how="left", on="node")

    value_name = str(endpoint_genes.columns[value_column - 1])
    expression = endpoint_genes.iloc[:, [0, value_column - 1]].rename(
        columns={endpoint_genes.columns[0]: "gene_symbol"}
    )
    node_table = node_table.merge(
        expression[["gene_symbol", value_name]], how="left", on="gene_symbol"
    )

    return node_table


# --- Orchestrator (the public entry point) ---------------------------------------------------------


def _read_tiedie_outputs(output_folder: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read TieDie's causal subnetwork and linker heats back into data frames."""

    causal_sif = pd.read_csv(
        output_folder / "tiedie.cn.sif",
        sep="\t",
        header=None,
        names=CAUSAL_SIF_COLUMNS,
    )

    heats = pd.read_csv(output_folder / "heats.NA", sep="=")
    heats = heats.reset_index()
    heats["index"] = heats["index"].str.replace(" ", "", regex=False)
    heats = heats.rename(columns={"index": "node"})

    return causal_sif, heats


def _prepare_tiedie_inputs(
    dmi_table: pd.DataFrame,
    ddi_table: pd.DataFrame | None,
    expressed_genes: list[str],
    endpoint_genes: pd.DataFrame,
    endpoint_value_column: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build step-1 outputs: (hmi_pairs, contextualised_tf_target, pathway_sif, upstream, downstream)."""

    hmi_pairs = combine_host_microbe_interactions(dmi_table, ddi_table)
    ppi, tf_target = fetch_omnipath_networks()
    contextualised_ppi, contextualised_tf_target = contextualise_networks(
        ppi, tf_target, expressed_genes, endpoint_genes
    )

    pathway_sif = build_pathway_sif(contextualised_ppi)
    upstream = build_upstream_heats(hmi_pairs)
    downstream = build_downstream_heats(
        contextualised_tf_target, endpoint_genes, endpoint_value_column
    )
    return hmi_pairs, contextualised_tf_target, pathway_sif, upstream, downstream


def _run_tiedie_in_workdir(
    upstream: pd.DataFrame,
    downstream: pd.DataFrame,
    pathway_sif: pd.DataFrame,
    work_dir: str | None,
    size: float,
    alpha: float | None,
    depth: int,
    permute: int,
    use_pagerank: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run TieDie in a working directory (a TemporaryDirectory if work_dir is None) and read it back.

    Returns the (causal_sif, heats) pair, read while the directory is still alive.
    """

    if work_dir is None:
        work_context = tempfile.TemporaryDirectory()
    else:
        Path(work_dir).mkdir(parents=True, exist_ok=True)
        work_context = contextlib.nullcontext(work_dir)

    with work_context as active_dir:
        output_folder = run_tiedie(
            upstream,
            downstream,
            pathway_sif,
            Path(active_dir),
            size=size,
            alpha=alpha,
            depth=depth,
            permute=permute,
            use_pagerank=use_pagerank,
        )
        return _read_tiedie_outputs(output_folder)


def run_tiedie_pipeline(
    dmi_table: pd.DataFrame,
    endpoint_genes: pd.DataFrame,
    expressed_genes: list[str],
    endpoint_value_column: int,
    ddi_table: pd.DataFrame | None = None,
    work_dir: str | None = None,
    size: float = 1.0,
    alpha: float | None = None,
    depth: int = 3,
    permute: int = 1000,
    use_pagerank: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run all three TieDie steps and return (final_network, node_table).

    Combines the DMI table with the optional DDI table into the host-microbe pair set, fetches and
    contextualises the OmniPath networks, builds the three TieDie inputs, runs TieDie in a working
    directory, and assembles the final network and node tables. The DEG table is used as supplied
    (no p-value filtering; Decision 11).

    Args:
        dmi_table: A Module 6/7/8 DMI table (only human_uniprot_id/bacterial_uniprot_id are read).
        endpoint_genes: The DEG table, filtered to significant genes (first column = gene symbol).
        expressed_genes: Gene symbols with non-zero expression (contextualisation universe).
        endpoint_value_column: 1-based column of the log2FC value in endpoint_genes.
        ddi_table: An optional Module 5 DDI table, merged into the pair set (Decision 3).
        work_dir: Working directory for TieDie I/O; a TemporaryDirectory if None.
        size, alpha, depth, permute, use_pagerank: TieDie algorithm parameters (passed through).

    Returns:
        (final_network, node_table) as data frames. Raises RuntimeError if TieDie fails.
    """

    hmi_pairs, tf_net, sif, upstream, downstream = _prepare_tiedie_inputs(
        dmi_table, ddi_table, expressed_genes, endpoint_genes, endpoint_value_column
    )

    causal_sif, heats = _run_tiedie_in_workdir(
        upstream, downstream, sif, work_dir, size, alpha, depth, permute, use_pagerank
    )

    final_network = assemble_network(causal_sif, hmi_pairs, tf_net)
    node_table = assemble_node_table(
        final_network, heats, endpoint_genes, endpoint_value_column
    )
    return final_network, node_table
