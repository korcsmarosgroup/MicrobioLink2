"""Console-script entry point for the TieDie tied-diffusion pipeline. Argparse boilerplate only."""

import argparse

import pandas as pd


def _add_tiedie_source_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the DMI/DDI tables, DEG list, and transcriptomics source."""

    parser.add_argument(
        "--dmi_file",
        required=True,
        help="DMI predictions CSV (microbiolink-dmi/-idr-filter/-monte-carlo output).",
    )
    parser.add_argument(
        "--ddi_file",
        help="Optional DDI predictions CSV (microbiolink-ddi output), merged into the pair set.",
    )
    parser.add_argument(
        "--deg_file",
        required=True,
        help="Endpoint DEG table, pre-filtered to significant genes (first column = gene symbol).",
    )
    parser.add_argument(
        "--deg_sep",
        default=",",
        help="Field separator for --deg_file.",
    )
    parser.add_argument(
        "--deg_value_column",
        type=int,
        required=True,
        help="One-based log2FC/expression value column in --deg_file.",
    )
    parser.add_argument(
        "--transcriptomics_file",
        required=True,
        help="Expressed-gene source (count matrix; first column = gene symbol).",
    )
    parser.add_argument(
        "--transcriptomics_sep",
        default=",",
        help="Field separator for --transcriptomics_file.",
    )
    parser.add_argument(
        "--transcriptomics_value_column",
        type=int,
        required=True,
        help="One-based expression value column used for the non-zero expression filter.",
    )


def _add_tiedie_algorithm_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the TieDie algorithm parameters and working directory."""

    parser.add_argument(
        "--size",
        type=float,
        default=1.0,
        help="TieDie network size-control factor.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        help="Optional TieDie linker-cutoff override (overrides --size when set).",
    )
    parser.add_argument(
        "--depth",
        type=int,
        default=3,
        help="TieDie causal-path search depth.",
    )
    parser.add_argument(
        "--permute",
        type=int,
        default=1000,
        help="Number of permutations for the TieDie significance test.",
    )
    parser.add_argument(
        "--pagerank",
        action="store_true",
        help="Diffuse with Personalized PageRank instead of the heat kernel.",
    )
    parser.add_argument(
        "--work_dir",
        help="Keep TieDie's intermediate files in this directory instead of a temp dir.",
    )


def _build_tiedie_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the tiedie CLI command."""

    parser = argparse.ArgumentParser(
        description="Run the three-step TieDie tied-diffusion network propagation end to end.",
    )
    _add_tiedie_source_arguments(parser)
    _add_tiedie_algorithm_arguments(parser)
    parser.add_argument(
        "--network_output",
        required=True,
        help="Output TSV for the final host-microbe signalling network.",
    )
    parser.add_argument(
        "--node_output",
        required=True,
        help="Output TSV for the per-node annotation table.",
    )
    return parser


def _read_expressed_genes(file_path: str, sep: str, value_column: int) -> list[str]:
    """Read gene symbols with present, non-NaN, non-zero expression (the non-zero expression filter)."""

    matrix = pd.read_csv(file_path, sep=sep)
    values = pd.to_numeric(matrix.iloc[:, value_column - 1], errors="coerce")
    mask = values.notna() & (values != 0)
    return matrix.loc[mask, matrix.columns[0]].tolist()


def main() -> int:
    """Run the TieDie tied-diffusion pipeline and write the network and node tables."""

    from ..workflow import tiedie as tiedie_workflow

    args = _build_tiedie_parser().parse_args()
    dmi_table = pd.read_csv(args.dmi_file)
    ddi_table = pd.read_csv(args.ddi_file) if args.ddi_file else None
    endpoint_genes = pd.read_csv(args.deg_file, sep=args.deg_sep)
    expressed_genes = _read_expressed_genes(
        args.transcriptomics_file,
        args.transcriptomics_sep,
        args.transcriptomics_value_column,
    )

    network, node_table = tiedie_workflow.run_tiedie_pipeline(
        dmi_table,
        endpoint_genes,
        expressed_genes,
        args.deg_value_column,
        ddi_table=ddi_table,
        work_dir=args.work_dir,
        size=args.size,
        alpha=args.alpha,
        depth=args.depth,
        permute=args.permute,
        use_pagerank=args.pagerank,
    )
    network.to_csv(args.network_output, sep="\t", index=False)
    node_table.to_csv(args.node_output, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
