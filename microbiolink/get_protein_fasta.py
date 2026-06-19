#!/usr/bin/env python

"""Download FASTA sequences from UniProt for any list of protein accessions."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests


UNIPROT_FASTA_URL = (
    'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=accession:('
)
DEFAULT_BATCH_SIZE = 100


def fetch_fasta_sequences(uniprot_ids: list[str]) -> str:
    """Fetch FASTA sequences for a list of UniProt accessions.

    Args:
        uniprot_ids: List of UniProt accession strings.

    Returns:
        FASTA-formatted sequence text.

    Raises:
        requests.HTTPError: If the UniProt request fails.
    """
    url = UNIPROT_FASTA_URL + '+OR+'.join(uniprot_ids) + ')'
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.text


def read_ids(
    filename: str | Path,
    separator: str,
    id_column: int,
) -> list[str]:
    """Read UniProt identifiers from a delimited file.

    Args:
        filename: Path to the identifier file.
        separator: Field separator used in the file.
        id_column: One-based column index containing the identifiers.

    Returns:
        List of identifier strings, one per data row.
    """
    ids: list[str] = []
    column_index = id_column - 1

    with open(Path(filename), encoding='utf-8-sig') as id_list:
        next(id_list, None)

        for line_number, line in enumerate(id_list, start=2):
            fields = line.strip().split(separator)

            if column_index >= len(fields):
                raise ValueError(
                    f'Column {id_column} is out of range on line {line_number}.',
                )

            ids.append(fields[column_index])

    return ids


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description='Download FASTA sequences for a list of UniProt accessions.',
    )
    parser.add_argument(
        '--id_list',
        required=True,
        help='Path to a file containing UniProt accessions.',
    )
    parser.add_argument(
        '--sep',
        required=True,
        help='Field separator in the identifier file.',
    )
    parser.add_argument(
        '--id_column',
        type=int,
        required=True,
        help='One-based column number containing the identifiers.',
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Path to the output FASTA file.',
    )
    parser.add_argument(
        '--batch_size',
        type=int,
        default=DEFAULT_BATCH_SIZE,
        help=f'Number of accessions per request (default: {DEFAULT_BATCH_SIZE}).',
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    """Run the FASTA download workflow.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Exit code.
    """
    args = parse_args(argv)
    ids = read_ids(args.id_list, args.sep, args.id_column)

    with open(Path(args.output), 'w', encoding='utf-8') as output_file:
        for start in range(0, len(ids), args.batch_size):
            batch = ids[start : start + args.batch_size]
            output_file.write(fetch_fasta_sequences(batch))

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
