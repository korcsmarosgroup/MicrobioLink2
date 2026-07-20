"""Console-script entry points for MicrobioLink. Argparse boilerplate only."""

import argparse

import pandas as pd


def zscore_filter() -> int:
    """Filter a gene/protein count matrix by a z-score cutoff."""

    from .workflow import zscore_filter as zscore_filter_module

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


def _add_identifier_source_arguments(parser: argparse.ArgumentParser) -> None:
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
    _add_identifier_source_arguments(parser)
    _add_membrane_filter_target_arguments(parser)
    return parser


def _resolve_membrane_filter_identifiers(args: argparse.Namespace) -> list[str]:
    """Resolve the identifier list from CLI args, per --from-count-matrix."""

    from .utils import uniprot_client

    if args.from_count_matrix:
        count_matrix = pd.read_csv(args.input_file, index_col=0)
        return count_matrix.index.tolist()

    return uniprot_client.read_ids(args.input_file, args.sep, args.id_column)


def membrane_filter() -> int:
    """Filter human or microbial proteins to membrane/secreted proteins."""

    from .workflow import membrane_filter as membrane_filter_module

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


def _add_human_fasta_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments for the human fasta identifiers (all optional)."""

    parser.add_argument(
        '-hi',
        '--human_input_file',
        help='Human identifier list file, or a count matrix if --human-from-count-matrix is set.',
    )
    parser.add_argument(
        '-hid',
        '--human_id_type',
        choices=['uniprot', 'genesymbol'],
        help='Identifier type in --human_input_file.',
    )
    parser.add_argument(
        '--human-from-count-matrix',
        action='store_true',
        help='Treat --human_input_file as a gene count matrix and extract expressed gene symbols from it.',
    )
    parser.add_argument(
        '-hsep',
        '--human_sep',
        default=',',
        help='Field separator for a plain human identifier list.',
    )
    parser.add_argument(
        '-hcol',
        '--human_id_column',
        type=int,
        default=1,
        help='One-based identifier column for a plain human identifier list.',
    )


def _add_microbial_fasta_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments for the microbial fasta identifiers (all optional)."""

    parser.add_argument(
        '-mi',
        '--microbial_input_file',
        help='Microbial identifier list file.',
    )
    parser.add_argument(
        '-mid',
        '--microbial_id_type',
        choices=['uniprot', 'proteome'],
        help='Identifier type in --microbial_input_file.',
    )
    parser.add_argument(
        '-msep',
        '--microbial_sep',
        default=',',
        help='Field separator for the microbial identifier list.',
    )
    parser.add_argument(
        '-mcol',
        '--microbial_id_column',
        type=int,
        default=1,
        help='One-based identifier column for the microbial identifier list.',
    )


def _build_fasta_download_parser() -> argparse.ArgumentParser:
    """Build the argument parser for the fasta-download CLI command."""

    parser = argparse.ArgumentParser(
        description='Download human and/or microbial protein sequences as FASTA files.',
    )
    _add_human_fasta_arguments(parser)
    _add_microbial_fasta_arguments(parser)
    parser.add_argument(
        '-o',
        '--output_folder',
        required=True,
        help='Folder to write human_proteins.fasta / microbial_proteins.fasta into.',
    )
    return parser


def _resolve_human_fasta_identifiers(args: argparse.Namespace) -> list[str] | None:
    """Resolve the human identifier list from CLI args, or None if not supplied."""

    if args.human_input_file is None:
        return None

    from .utils import uniprot_client

    if args.human_from_count_matrix:
        count_matrix = pd.read_csv(args.human_input_file, index_col=0)
        expressed_mask = count_matrix.notna().any(axis=1) & (count_matrix != 0).any(axis=1)
        return count_matrix.index[expressed_mask].tolist()

    return uniprot_client.read_ids(args.human_input_file, args.human_sep, args.human_id_column)


def _resolve_microbial_fasta_identifiers(args: argparse.Namespace) -> list[str] | None:
    """Resolve the microbial identifier list from CLI args, or None if not supplied."""

    if args.microbial_input_file is None:
        return None

    from .utils import uniprot_client

    return uniprot_client.read_ids(args.microbial_input_file, args.microbial_sep, args.microbial_id_column)


def download_fasta() -> int:
    """Download human and/or microbial protein sequences as FASTA files."""

    from .workflow import fasta_download

    args = _build_fasta_download_parser().parse_args()
    human_identifiers = _resolve_human_fasta_identifiers(args)
    microbial_identifiers = _resolve_microbial_fasta_identifiers(args)

    fasta_download.download_fasta(
        args.output_folder,
        human_identifiers=human_identifiers,
        human_id_type=args.human_id_type,
        microbial_identifiers=microbial_identifiers,
        microbial_id_type=args.microbial_id_type,
    )
    return 0
