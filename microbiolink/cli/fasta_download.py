"""Console-script entry point for downloading protein FASTA sequences. Argparse boilerplate only."""

import argparse

from ._common import (
    _add_human_identifier_arguments,
    _add_microbial_identifier_arguments,
    _resolve_human_identifiers,
    _resolve_microbial_identifiers,
)


def _build_fasta_download_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the fasta-download CLI command."""

    parser = argparse.ArgumentParser(
        description="Download human and/or microbial protein sequences as FASTA files.",
    )
    _add_human_identifier_arguments(parser)
    _add_microbial_identifier_arguments(parser)
    parser.add_argument(
        "-o",
        "--output_folder",
        required=True,
        help="Folder to write human_proteins.fasta / microbial_proteins.fasta into.",
    )
    return parser


def main() -> int:
    """Download human and/or microbial protein sequences as FASTA files."""

    from ..workflow import fasta_download

    args = _build_fasta_download_parser().parse_args()
    human_identifiers = _resolve_human_identifiers(args)
    microbial_identifiers = _resolve_microbial_identifiers(args)

    fasta_download.download_fasta(
        args.output_folder,
        human_identifiers=human_identifiers,
        human_id_type=args.human_id_type,
        microbial_identifiers=microbial_identifiers,
        microbial_id_type=args.microbial_id_type,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
