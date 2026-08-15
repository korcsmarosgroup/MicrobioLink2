"""Console-script entry point for the IDR (disorder/binding) filter. Argparse boilerplate only."""

import argparse

import pandas as pd


def _add_idr_filter_source_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the DMI table and its source FASTA files."""

    parser.add_argument(
        "--dmi_file",
        required=True,
        help="DMI predictions CSV file (microbiolink-dmi output).",
    )
    parser.add_argument(
        "--human_fasta_file",
        help="Human protein FASTA file (microbiolink-download-fasta output). Required if dmi_file has 'forward' rows.",
    )
    parser.add_argument(
        "--bacterial_fasta_file",
        help="Bacterial protein FASTA file (microbiolink-download-fasta output). Required if dmi_file has 'reverse' rows.",
    )


def _add_idr_filter_scoring_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the disorder/binding method, cutoffs, and device."""

    parser.add_argument(
        "--method",
        required=True,
        choices=["iupred", "aiupred"],
        help="Disorder/binding prediction method.",
    )
    parser.add_argument(
        "--disorder_cutoff",
        required=True,
        type=float,
        help="Minimum per-residue disorder score required across the whole motif window.",
    )
    parser.add_argument(
        "--binding_cutoff",
        required=True,
        type=float,
        help="Minimum per-residue binding score required across the whole motif window.",
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


def _build_idr_filter_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the idr-filter CLI command."""

    parser = argparse.ArgumentParser(
        description="Filter domain-motif interactions to those in a disordered, binding-prone region.",
    )
    _add_idr_filter_source_arguments(parser)
    _add_idr_filter_scoring_arguments(parser)
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for the filtered domain-motif interactions.",
    )
    return parser


def main() -> int:
    """Filter domain-motif interactions to those in a disordered, binding-prone region."""

    from ..utils import fasta
    from ..workflow import idr_filter as idr_filter_module

    args = _build_idr_filter_parser().parse_args()
    dmi_table = pd.read_csv(args.dmi_file)
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

    result = idr_filter_module.filter_by_disorder(
        dmi_table,
        human_sequences=human_sequences,
        bacterial_sequences=bacterial_sequences,
        method=args.method,
        disorder_cutoff=args.disorder_cutoff,
        binding_cutoff=args.binding_cutoff,
        force_cpu=args.force_cpu,
        gpu_num=args.gpu_num,
    )
    result.to_csv(args.output_file, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
