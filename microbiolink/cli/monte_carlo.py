"""Console-script entry point for the Monte Carlo over-representation filter. Argparse boilerplate only."""

import argparse

import pandas as pd


def _add_monte_carlo_source_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the interaction table and its source FASTA files."""

    parser.add_argument(
        "--interaction_file",
        required=True,
        help="IDR-filtered interaction CSV file (microbiolink-idr-filter output).",
    )
    parser.add_argument(
        "--human_fasta_file",
        help="Human protein FASTA file (microbiolink-download-fasta output). Required if the table has 'forward' rows.",
    )
    parser.add_argument(
        "--bacterial_fasta_file",
        help="Bacterial protein FASTA file (microbiolink-download-fasta output). Required if the table has 'reverse' rows.",
    )


def _add_monte_carlo_disorder_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the disorder method, cutoff, and device."""

    parser.add_argument(
        "--method",
        required=True,
        choices=["iupred", "aiupred"],
        help="Disorder prediction method (locates each motif's disordered region and the pool).",
    )
    parser.add_argument(
        "--disorder_cutoff",
        required=True,
        type=float,
        help="Per-residue disorder score at/above which a residue is disordered.",
    )
    parser.add_argument(
        "--force-cpu",
        dest="force_cpu",
        action="store_true",
        help="Force CPU inference (aiupred only; ignored for iupred).",
    )
    parser.add_argument(
        "--gpu-num",
        dest="gpu_num",
        type=int,
        default=0,
        help="GPU index to use (aiupred only; ignored for iupred).",
    )


def _add_monte_carlo_test_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the Monte Carlo test itself."""

    parser.add_argument(
        "--iterations",
        type=int,
        default=1000,
        help="Samples per motif instance. The MC p-value floor is 1/(iterations + 1); for m instances at FDR alpha, resolving the BH tail needs iterations >~ m/alpha, so raise this for large runs.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Target false discovery rate; a motif passes if its BH-corrected q-value <= alpha.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Random seed for reproducible sampling.",
    )


def _build_monte_carlo_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the monte-carlo CLI command."""

    parser = argparse.ArgumentParser(
        description="Filter domain-motif interactions by a Monte Carlo over-representation test against pooled disorder.",
    )
    _add_monte_carlo_source_arguments(parser)
    _add_monte_carlo_disorder_arguments(parser)
    _add_monte_carlo_test_arguments(parser)
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for the filtered domain-motif interactions.",
    )
    return parser


def main() -> int:
    """Filter domain-motif interactions by a Monte Carlo over-representation test."""

    from ..utils import fasta
    from ..workflow import monte_carlo

    args = _build_monte_carlo_parser().parse_args()
    interaction_table = pd.read_csv(args.interaction_file)
    human_sequences = (
        fasta.read_fasta_sequences(args.human_fasta_file)
        if args.human_fasta_file
        else None
    )
    bacterial_sequences = (
        fasta.read_fasta_sequences(args.bacterial_fasta_file)
        if args.bacterial_fasta_file
        else None
    )

    result = monte_carlo.filter_by_monte_carlo(
        interaction_table,
        human_sequences=human_sequences,
        bacterial_sequences=bacterial_sequences,
        method=args.method,
        disorder_cutoff=args.disorder_cutoff,
        iterations=args.iterations,
        alpha=args.alpha,
        seed=args.seed,
        force_cpu=args.force_cpu,
        gpu_num=args.gpu_num,
    )
    result.to_csv(args.output_file, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
