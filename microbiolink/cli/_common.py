"""CLI helpers shared across more than one command module."""

import argparse

import pandas as pd


def _add_human_identifier_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments for the human fasta identifiers (all optional)."""

    parser.add_argument(
        "-hi",
        "--human_input_file",
        help="Human identifier list file, or a count matrix if --human-from-count-matrix is set.",
    )
    parser.add_argument(
        "-hid",
        "--human_id_type",
        choices=["uniprot", "genesymbol"],
        help="Identifier type in --human_input_file.",
    )
    parser.add_argument(
        "--human-from-count-matrix",
        action="store_true",
        help="Treat --human_input_file as a gene count matrix and extract expressed gene symbols from it.",
    )
    parser.add_argument(
        "-hsep",
        "--human_sep",
        default=",",
        help="Field separator for a plain human identifier list.",
    )
    parser.add_argument(
        "-hcol",
        "--human_id_column",
        type=int,
        default=1,
        help="One-based identifier column for a plain human identifier list.",
    )
    parser.add_argument(
        "--human-no-header",
        action="store_true",
        help="Treat the human identifier list as having no header row.",
    )


def _add_microbial_identifier_arguments(parser: argparse.ArgumentParser) -> None:
    """Add CLI arguments for the microbial fasta identifiers (all optional)."""

    parser.add_argument(
        "-mi",
        "--microbial_input_file",
        help="Microbial identifier list file.",
    )
    parser.add_argument(
        "-mid",
        "--microbial_id_type",
        choices=["uniprot", "proteome"],
        help="Identifier type in --microbial_input_file.",
    )
    parser.add_argument(
        "-msep",
        "--microbial_sep",
        default=",",
        help="Field separator for the microbial identifier list.",
    )
    parser.add_argument(
        "-mcol",
        "--microbial_id_column",
        type=int,
        default=1,
        help="One-based identifier column for the microbial identifier list.",
    )
    parser.add_argument(
        "--microbial-no-header",
        action="store_true",
        help="Treat the microbial identifier list as having no header row.",
    )


def _resolve_human_identifiers(args: argparse.Namespace) -> list[str] | None:
    """Resolve the human identifier list from CLI args, or None if not supplied."""

    if args.human_input_file is None:
        return None

    from ..utils import uniprot_client

    if args.human_from_count_matrix:
        count_matrix = pd.read_csv(args.human_input_file, index_col=0)
        expressed_mask = count_matrix.notna().any(axis=1) & (count_matrix != 0).any(
            axis=1
        )
        return count_matrix.index[expressed_mask].tolist()

    return uniprot_client.read_ids(
        args.human_input_file,
        args.human_sep,
        args.human_id_column,
        has_header=not args.human_no_header,
    )


def _resolve_microbial_identifiers(args: argparse.Namespace) -> list[str] | None:
    """Resolve the microbial identifier list from CLI args, or None if not supplied."""

    if args.microbial_input_file is None:
        return None

    from ..utils import uniprot_client

    return uniprot_client.read_ids(
        args.microbial_input_file,
        args.microbial_sep,
        args.microbial_id_column,
        has_header=not args.microbial_no_header,
    )


def _read_domain_mapping(filename: str) -> dict[str, list[str]]:
    """Read a Pfam/Entries domain TSV into a pfam_id -> uniprot_ids mapping."""

    domain_table = pd.read_csv(filename, sep="\t")
    return {
        row["Pfam"]: [entry for entry in row["Entries"].split(";") if entry]
        for _, row in domain_table.iterrows()
    }
