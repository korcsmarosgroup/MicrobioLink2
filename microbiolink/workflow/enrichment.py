"""Functional enrichment analysis: reduce a network/HMI table to human targets, run Enrichr, plot."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ..utils import id_resolution

if TYPE_CHECKING:
    from matplotlib.figure import Figure

HMI_HUMAN_COLUMN = "human_uniprot_id"
NETWORK_SOURCE_COLUMN = "Source.node"
NETWORK_TARGET_COLUMN = "Target.node"
DEFAULT_DATABASE = "Reactome_2022"
SIGNIFICANCE_CUTOFF = 0.05
TOP_N = 20
HIGHLIGHT_COLOUR = "#FF8800"
RESULT_COLUMNS = [
    "rank", "path_name", "p_val", "z_score",
    "combined_score", "overlapping_genes", "adj_p_val", "database",
]  # fmt: skip


# --- Target extraction (one per analysis level) ----------------------------------------------------


def extract_hmi_targets(hmi_table: pd.DataFrame) -> list[str]:
    """Distinct human UniProt targets from a Module 5/6/7/8 table (the human_uniprot_id column).

    Read by column name, not tab position. The human_uniprot_id column holds the human target of
    the microbial protein in both DMI directions (Module 9 established this), so this equals the
    ground-truth target set (its bacterial column provably drops in the human-only translation).

    Args:
        hmi_table: A Module 5 DDI table or a Module 6/7/8 DMI table.

    Returns:
        The distinct human UniProt accessions.
    """

    return hmi_table[HMI_HUMAN_COLUMN].drop_duplicates().tolist()


def extract_tiedie_nodes(network: pd.DataFrame) -> list[str]:
    """Every node of a Module 9 final network: the union of Source.node and Target.node.

    Fixes the ground-truth Target-only bug (it read Target.node plus Relationship, silently
    dropping the Source-only nodes). Bacterial source nodes are included and drop out naturally in
    the human-only translation.

    Args:
        network: A Module 9 final network table (columns Target.node, Source.node, ...).

    Returns:
        The distinct UniProt accessions across both node columns.
    """

    nodes = pd.concat(
        [network[NETWORK_SOURCE_COLUMN], network[NETWORK_TARGET_COLUMN]],
        ignore_index=True,
    )
    return nodes.drop_duplicates().tolist()


# --- Enrichment ------------------------------------------------------------------------------------


def run_enrichment(
    target_symbols: list[str],
    background_symbols: list[str],
    database: str = DEFAULT_DATABASE,
) -> pd.DataFrame:
    """Run gget.enrichr against `database` with the background list; filter to adj_p_val < 0.05.

    Lazy `import gget`; called with plot=False, species="human". gget is a thin client over the
    live Enrichr API, so this makes a network request.

    Args:
        target_symbols: Human gene symbols to test for over-representation.
        background_symbols: The expressed-gene-symbol universe used as the enrichment background.
        database: Enrichr gene-set library (or a gget shortcut alias).

    Returns:
        The significant results table (RESULT_COLUMNS schema); may be empty if nothing passes the
        adj_p_val < 0.05 cutoff.
    """

    import gget

    results = gget.enrichr(
        target_symbols,
        database=database,
        background_list=background_symbols,
        species="human",
        plot=False,
    )
    # gget.enrichr is typed to also return None / a JSON list; plot=False always yields a frame,
    # but guard defensively so a None (no Enrichr hit) collapses to an empty, correctly-typed table.
    if not isinstance(results, pd.DataFrame):
        return pd.DataFrame(columns=RESULT_COLUMNS)
    return results[results["adj_p_val"] < SIGNIFICANCE_CUTOFF]


# --- Plot (object-oriented Figure API; no pyplot) --------------------------------------------------


def _plot_combined_score(results: pd.DataFrame, ax) -> None:
    """Draw the ground-truth twin-axis plot: bars by combined score, scatter of -log10(adj_p_val)."""

    # Live gget returns combined_score as an object column mixing floats and the string "inf", so
    # coerce before comparing/sorting (the ground truth's inf rows). np.isfinite drops inf and NaN.
    scores = pd.to_numeric(results["combined_score"], errors="coerce")
    finite = results[np.isfinite(scores)].copy()
    finite["combined_score"] = pd.to_numeric(finite["combined_score"])
    top = finite.nlargest(TOP_N, "combined_score").sort_values("combined_score")

    ax.barh(y=top["path_name"], width=top["combined_score"])
    ax.set_xlabel("Combined score")
    ax.set_title("Top 20 enriched pathways by combined score")

    significance_ax = ax.twiny()
    significance_ax.scatter(
        x=-np.log10(top["adj_p_val"]), y=top["path_name"], color=HIGHLIGHT_COLOUR
    )
    significance_ax.set_xlabel("-log 10(Adjusted p-value)", color=HIGHLIGHT_COLOUR)
    significance_ax.tick_params("x", colors=HIGHLIGHT_COLOUR)


def _plot_adj_p(results: pd.DataFrame, ax) -> None:
    """Draw our own horizontal bar of -log10(adj_p_val) for the top-20 terms by rank."""

    top = results.head(TOP_N)
    top = top.assign(neg_log_adj_p=-np.log10(top["adj_p_val"])).sort_values(
        "neg_log_adj_p"
    )

    ax.barh(y=top["path_name"], width=top["neg_log_adj_p"])
    ax.set_xlabel("-log 10(Adjusted p-value)")
    ax.set_title("Top 20 enriched pathways by adjusted p-value")


def plot_enrichment(
    results: pd.DataFrame, ranking: str = "combined_score"
) -> Figure | None:
    """Build the top-20 enrichment plot; return a matplotlib Figure, or None if `results` is empty.

    Lazy `from matplotlib.figure import Figure` (no pyplot, no backend selection). The results are
    sorted by the chosen ranking and the top 20 are shown.

    Args:
        results: The significant-enrichment table from run_enrichment.
        ranking: 'combined_score' (ground-truth twin-axis bars + -log10 adj_p scatter, inf rows
            excluded) or 'adj_p' (horizontal bars of -log10(adj_p_val)).

    Returns:
        A matplotlib Figure, or None when `results` is empty.

    Raises:
        ValueError: If ranking is not 'combined_score' or 'adj_p'.
    """

    if results.empty:
        return None

    from matplotlib.figure import Figure

    figure = Figure()
    ax = figure.subplots()

    if ranking == "combined_score":
        _plot_combined_score(results, ax)
    elif ranking == "adj_p":
        _plot_adj_p(results, ax)
    else:
        raise ValueError(
            f"ranking must be 'combined_score' or 'adj_p', got {ranking!r}"
        )

    return figure


# --- Orchestrator (the public entry point) ---------------------------------------------------------


def run_enrichment_analysis(
    target_table: pd.DataFrame,
    background_symbols: list[str],
    analysis_level: str,
    database: str = DEFAULT_DATABASE,
    ranking: str = "combined_score",
) -> tuple[pd.DataFrame, Figure | None]:
    """Extract targets for the given analysis_level, translate to symbols, enrich, and plot.

    Args:
        target_table: A Module 5/6/7/8 HMI table (analysis_level='HMI') or a Module 9 final
            network (analysis_level='TieDIE').
        background_symbols: Expressed-gene-symbol universe for the enrichment background.
        analysis_level: 'HMI' (distinct human_uniprot_id) or 'TieDIE' (all network nodes).
        database: Enrichr gene-set library (default Reactome_2022).
        ranking: 'combined_score' or 'adj_p' — how the returned plot is ranked.

    Returns:
        (results_table, figure) — the significant-enrichment table and its top-20 plot; figure is
        None when no term passes the significance cutoff.

    Raises:
        ValueError: If analysis_level is not 'HMI' or 'TieDIE'.
    """

    if analysis_level == "HMI":
        targets = extract_hmi_targets(target_table)
    elif analysis_level == "TieDIE":
        targets = extract_tiedie_nodes(target_table)
    else:
        raise ValueError(
            f"analysis_level must be 'HMI' or 'TieDIE', got {analysis_level!r}"
        )

    symbols = list(id_resolution.translate_uniprot_to_gene_symbols(targets).values())

    results = run_enrichment(symbols, background_symbols, database=database)
    figure = plot_enrichment(results, ranking=ranking)
    return results, figure
