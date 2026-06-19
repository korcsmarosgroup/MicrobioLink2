#!/usr/bin/env python

"""Download Pfam domain annotations from UniProt for any list of proteins."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests


UNIPROT_STREAM_BASE_URL = 'https://rest.uniprot.org/uniprotkb/stream?'
DEFAULT_UNIPROT_FIELDS = [
    'accession',
    'xref_pfam',
    'gene_names',
]
UNIPROT_BATCH_SIZE = 1000


def read_ids(
    filename: str | Path,
    separator: str,
    id_column: int,
) -> list[str]:
    """Read UniProt or proteome identifiers from a delimited file.

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


def build_uniprot_accession_query(uniprot_ids: list[str]) -> str:
    """Build a percent-encoded UniProt query for a list of accessions.

    Args:
        uniprot_ids: List of UniProt accession strings.

    Returns:
        A percent-encoded query string joining accessions with OR.
    """
    clauses = [f'%28accession%3A{uniprot_id}%29' for uniprot_id in uniprot_ids]
    return '%28' + '+OR+'.join(clauses) + '%29'


def build_uniprot_stream_url(
    query: str,
    fields: list[str] | None = None,
) -> str:
    """Build a UniProt stream endpoint URL.

    Args:
        query: Percent-encoded UniProt query string.
        fields: Return fields. Defaults to DEFAULT_UNIPROT_FIELDS.

    Returns:
        A complete UniProt stream URL.
    """
    selected_fields = fields or DEFAULT_UNIPROT_FIELDS
    encoded_fields = '%2C'.join(selected_fields)
    return (
        f'{UNIPROT_STREAM_BASE_URL}'
        f'fields={encoded_fields}&format=tsv&query={query}'
    )


def download_proteome_with_fields(
    proteome_id: str,
    fields: list[str] | None = None,
) -> str:
    """Download a full proteome domain table from UniProt with selected fields.

    Args:
        proteome_id: UniProt proteome identifier (e.g. UP000000625).
        fields: Return fields. Defaults to DEFAULT_UNIPROT_FIELDS.

    Returns:
        TSV response text from UniProt.
    """
    url = build_uniprot_stream_url(
        f'%28%28proteome%3A{proteome_id}%29%29',
        fields=fields,
    )
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.text


def download_protein_list_with_fields(
    uniprot_ids: list[str],
    fields: list[str] | None = None,
) -> str:
    """Download a domain table for a list of UniProt accessions with selected fields.

    Args:
        uniprot_ids: List of UniProt accession strings.
        fields: Return fields. Defaults to DEFAULT_UNIPROT_FIELDS.

    Returns:
        TSV response text from UniProt.
    """
    url = build_uniprot_stream_url(
        build_uniprot_accession_query(uniprot_ids),
        fields=fields,
    )
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.text


def download_proteome(proteome_id: str) -> str:
    """Download a full proteome domain table from UniProt.

    Args:
        proteome_id: UniProt proteome identifier.

    Returns:
        TSV response text from UniProt.
    """
    return download_proteome_with_fields(proteome_id)


def download_protein_list(uniprot_ids: list[str]) -> str:
    """Download a domain table for a list of UniProt accessions.

    Args:
        uniprot_ids: List of UniProt accession strings.

    Returns:
        TSV response text from UniProt.
    """
    return download_protein_list_with_fields(uniprot_ids)


def parse_args(argv: list[str]) -> argparse.Namespace:
    """Parse command-line arguments.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Parsed argument namespace.
    """
    parser = argparse.ArgumentParser(
        description='Download Pfam domain annotations for a list of UniProt accessions.',
    )
    parser.add_argument(
        '--id_list',
        required=True,
        help='Path to a file containing UniProt accessions or proteome IDs.',
    )
    parser.add_argument(
        '--sep',
        required=True,
        help='Field separator in the identifier file.',
    )
    parser.add_argument(
        '--id_type',
        choices=['uniprot', 'UP'],
        required=True,
        help='Identifier type: UniProt accessions or UniProt proteomes.',
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
        help='Path to the output TSV file.',
    )
    return parser.parse_args(argv)


def main(argv: list[str]) -> int:
    """Run the UniProt domain download workflow.

    Args:
        argv: Argument list (typically sys.argv[1:]).

    Returns:
        Exit code.
    """
    args = parse_args(argv)
    ids = read_ids(args.id_list, args.sep, args.id_column)
    header_written = False

    with open(Path(args.output), 'w', encoding='utf-8') as output_file:
        for start in range(0, len(ids), UNIPROT_BATCH_SIZE):
            batch = ids[start : start + UNIPROT_BATCH_SIZE]

            if args.id_type == 'UP':
                for proteome_id in batch:
                    lines = download_proteome(proteome_id).rstrip().split('\n')

                    if not lines or not lines[0]:
                        continue

                    if not header_written:
                        output_file.write(lines[0] + '\tProteome_ID\n')
                        header_written = True

                    for row in lines[1:]:
                        if row:
                            output_file.write(f'{row}\t{proteome_id}\n')
                continue

            lines = download_protein_list(batch).rstrip().split('\n')

            if not lines or not lines[0]:
                continue

            if not header_written:
                output_file.write(lines[0] + '\n')
                header_written = True

            for row in lines[1:]:
                if row:
                    output_file.write(row + '\n')

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
