"""Console-script entry point for the membrane/secreted protein filter. Argparse boilerplate only."""

import argparse

import pandas as pd


def _add_identifier_source_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing where identifiers come from."""

    parser.add_argument(
        "-i",
        "--input_file",
        required=True,
        help="Identifier list file, or a count matrix if --from-count-matrix is set.",
    )
    parser.add_argument(
        "--from-count-matrix",
        action="store_true",
        help="Treat --input_file as a gene count matrix and extract gene symbols from it.",
    )
    parser.add_argument(
        "-sep",
        "--sep",
        default=",",
        help="Field separator for a plain identifier list (ignored with --from-count-matrix).",
    )
    parser.add_argument(
        "-col",
        "--id_column",
        type=int,
        default=1,
        help="One-based identifier column for a plain list (ignored with --from-count-matrix).",
    )
    parser.add_argument(
        "--no-header",
        action="store_true",
        help="Treat the identifier list as having no header row (ignored with --from-count-matrix).",
    )


def _add_membrane_filter_target_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the filter and output."""

    parser.add_argument(
        "-id",
        "--id_type",
        required=True,
        choices=["uniprot", "genesymbol", "proteome"],
        help="Identifier type in the input file.",
    )
    parser.add_argument(
        "-sp",
        "--species",
        required=True,
        choices=["human", "microbial"],
        help="Species the identifiers belong to.",
    )
    parser.add_argument(
        "-lfl",
        "--location_filters",
        required=True,
        nargs="+",
        help="Location categories (human) or location substrings (microbial) to keep.",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        required=True,
        help="Output CSV file for the filtered results.",
    )


def _build_membrane_filter_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the membrane-filter CLI command."""

    parser = argparse.ArgumentParser(
        description="Filter human or microbial proteins down to membrane/secreted proteins.",
    )
    _add_identifier_source_arguments(parser)
    _add_membrane_filter_target_arguments(parser)
    return parser


def _resolve_membrane_filter_identifiers(args: argparse.Namespace) -> list[str]:
    """Resolve the identifier list from CLI args, per --from-count-matrix."""

    from ..utils import uniprot_client

    if args.from_count_matrix:
        count_matrix = pd.read_csv(args.input_file, index_col=0)
        return count_matrix.index.tolist()

    return uniprot_client.read_ids(
        args.input_file, args.sep, args.id_column, has_header=not args.no_header
    )


def main() -> int:
    """Filter human or microbial proteins to membrane/secreted proteins."""

    from ..workflow import membrane_filter as membrane_filter_module

    args = _build_membrane_filter_parser().parse_args()
    identifiers = _resolve_membrane_filter_identifiers(args)

    result = membrane_filter_module.filter_membrane_proteins(
        identifiers,
        id_type=args.id_type,
        species=args.species,
        location_filters=args.location_filters,
    )
    result.to_csv(args.output_file, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
