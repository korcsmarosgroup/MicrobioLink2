#!/usr/bin/env python

"""Organism-agnostic UniProt REST client primitives."""

from __future__ import annotations

from pathlib import Path

import requests


PathLike = str | Path

UNIPROT_STREAM_BASE_URL = 'https://rest.uniprot.org/uniprotkb/stream?'
DEFAULT_UNIPROT_FIELDS = [
    'accession',
    'xref_pfam',
    'gene_names',
]
UNIPROT_BATCH_SIZE = 1000

UNIPROT_FASTA_URL = (
    'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=accession:('
)
DEFAULT_BATCH_SIZE = 100

UNIPROT_FASTA_PROTEOME_URL = (
    'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query='
    '%28proteome%3A{proteome_id}%29'
)


def read_ids(
    filename: PathLike,
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

    with open(Path(filename), encoding = 'utf-8-sig') as id_list:
        next(id_list, None)

        for line_number, line in enumerate(id_list, start = 2):
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
        fields = fields,
    )
    response = requests.get(url, timeout = 60)
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
        fields = fields,
    )
    response = requests.get(url, timeout = 60)
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
    response = requests.get(url, timeout = 60)
    response.raise_for_status()
    return response.text


def fetch_proteome_fasta(proteome_id: str) -> str:
    """Download FASTA sequences for a complete UniProt proteome.

    Args:
        proteome_id: UniProt proteome identifier (e.g. UP000000625).

    Returns:
        FASTA-formatted sequence text for the entire proteome.

    Raises:
        requests.HTTPError: If the UniProt request fails.
    """

    url = UNIPROT_FASTA_PROTEOME_URL.format(proteome_id = proteome_id)
    response = requests.get(url, timeout = 120)
    response.raise_for_status()
    return response.text
