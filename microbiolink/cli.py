"""Console-script entry points for MicrobioLink. Argparse boilerplate only."""

import argparse


def zscore_filter() -> int:
    """Filter a gene/protein count matrix by a z-score cutoff."""

    from . import zscore_filter as zscore_filter_module

    parser = argparse.ArgumentParser(
        description='Gene expression filtration based on individual cell count and z-score.',
    )
    parser.add_argument(
        '-i',
        '--input_file',
        required=True,
        help='Input CSV file with gene expression data.',
    )
    parser.add_argument(
        '-zscore',
        '--zscore',
        required=True,
        type=float,
        help='Z-score cut-off to filter lowly expressed genes.',
    )
    parser.add_argument(
        '-o',
        '--output_file',
        required=True,
        help='Output CSV file for filtered results.',
    )
    args = parser.parse_args()

    zscore_filter_module.filter_count_matrix_file(
        args.input_file,
        zscore_threshold=args.zscore,
        output_file=args.output_file,
    )
    return 0


def _add_membrane_filter_source_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing where identifiers come from."""

    parser.add_argument(
        '-i',
        '--input_file',
        required=True,
        help='Identifier list file, or a count matrix if --from-count-matrix is set.',
    )
    parser.add_argument(
        '--from-count-matrix',
        action='store_true',
        help='Treat --input_file as a gene count matrix and extract gene symbols from it.',
    )
    parser.add_argument(
        '-sep',
        '--sep',
        default=',',
        help='Field separator for a plain identifier list (ignored with --from-count-matrix).',
    )
    parser.add_argument(
        '-col',
        '--id_column',
        type=int,
        default=1,
        help='One-based identifier column for a plain list (ignored with --from-count-matrix).',
    )


def _add_membrane_filter_target_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments describing the filter and output."""

    parser.add_argument(
        '-id',
        '--id_type',
        required=True,
        choices=['uniprot', 'genesymbol', 'proteome'],
        help='Identifier type in the input file.',
    )
    parser.add_argument(
        '-sp',
        '--species',
        required=True,
        choices=['human', 'microbial'],
        help='Species the identifiers belong to.',
    )
    parser.add_argument(
        '-lfl',
        '--location_filters',
        required=True,
        nargs='+',
        help='Location categories (human) or location substrings (microbial) to keep.',
    )
    parser.add_argument(
        '-o',
        '--output_file',
        required=True,
        help='Output CSV file for the filtered results.',
    )


def _build_membrane_filter_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the membrane-filter CLI command."""

    parser = argparse.ArgumentParser(
        description='Filter human or microbial proteins down to membrane/secreted proteins.',
    )
    _add_membrane_filter_source_arguments(parser)
    _add_membrane_filter_target_arguments(parser)
    return parser


def _resolve_membrane_filter_identifiers(args: argparse.Namespace) -> list[str]:
    """Resolve the identifier list from CLI args, per --from-count-matrix."""

    from . import gene_matrix
    from . import uniprot_client
    from . import zscore_filter

    if args.from_count_matrix:
        count_matrix = zscore_filter.read_count_matrix(args.input_file)
        return gene_matrix.extract_gene_symbols_from_count_matrix(count_matrix)

    return uniprot_client.read_ids(args.input_file, args.sep, args.id_column)


def membrane_filter() -> int:
    """Filter human or microbial proteins to membrane/secreted proteins."""

    from . import membrane_filter as membrane_filter_module

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
