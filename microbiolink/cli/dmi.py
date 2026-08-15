"""Console-script entry point for domain-motif interaction prediction. Argparse boilerplate only."""

import argparse

from ._common import _read_domain_mapping


def _build_dmi_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the dmi CLI command."""

    parser = argparse.ArgumentParser(
        description="Predict domain-motif interactions between bacterial and human proteins.",
    )
    parser.add_argument(
        "-m",
        "--mode",
        required=True,
        choices=["forward", "reverse", "both"],
        help="forward: human motif -> bacterial domain. reverse: bacterial motif -> human domain.",
    )
    parser.add_argument(
        "-hf",
        "--human_fasta_file",
        help="Human protein FASTA file (microbiolink-download-fasta output). Required for forward/both.",
    )
    parser.add_argument(
        "-b",
        "--bacterial_domain_file",
        help="Bacterial Pfam/Entries domain TSV file (microbiolink-download-domains output). Required for forward/both.",
    )
    parser.add_argument(
        "-bf",
        "--bacterial_fasta_file",
        help="Bacterial protein FASTA file (microbiolink-download-fasta output). Required for reverse/both.",
    )
    parser.add_argument(
        "-hu",
        "--human_domain_file",
        help="Human Pfam/Entries domain TSV file (microbiolink-download-domains output). Required for reverse/both.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for predicted domain-motif interactions.",
    )
    return parser


def main() -> int:
    """Predict domain-motif interactions between bacterial and human proteins."""

    from ..utils import fasta
    from ..workflow import dmi as dmi_workflow

    args = _build_dmi_parser().parse_args()
    human_sequences = (
        fasta.read_fasta_sequences(args.human_fasta_file)
        if args.human_fasta_file
        else None
    )
    bacterial_domains = (
        _read_domain_mapping(args.bacterial_domain_file)
        if args.bacterial_domain_file
        else None
    )
    bacterial_sequences = (
        fasta.read_fasta_sequences(args.bacterial_fasta_file)
        if args.bacterial_fasta_file
        else None
    )
    human_domains = (
        _read_domain_mapping(args.human_domain_file) if args.human_domain_file else None
    )

    result = dmi_workflow.predict_domain_motif_interactions(
        args.mode,
        human_sequences=human_sequences,
        bacterial_domains=bacterial_domains,
        bacterial_sequences=bacterial_sequences,
        human_domains=human_domains,
    )
    result.to_csv(args.output_file, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
