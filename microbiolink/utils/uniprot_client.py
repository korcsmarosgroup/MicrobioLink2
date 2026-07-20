"""Shared low-level UniProt REST client for identifier and proteome lookups."""

from io import StringIO
from pathlib import Path
from typing import Union

import pandas as pd
import requests

PathLike = Union[str, Path]

UNIPROT_STREAM_BASE_URL = 'https://rest.uniprot.org/uniprotkb/stream?'
DEFAULT_UNIPROT_FIELDS = ['accession', 'xref_pfam', 'gene_names']
UNIPROT_BATCH_SIZE = 1000


def read_ids(
    filename: PathLike,
    separator: str,
    id_column: int,
) -> list[str]:
    """Read UniProt or proteome identifiers from a delimited file.

    Args:
        filename: Path to the identifier file. The first line is treated as
            a header and skipped.
        separator: Field separator used in the file.
        id_column: One-based column number containing the identifiers.

    Returns:
        A list of identifiers in file order.
    """

    column_index = id_column - 1
    ids = []

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

    encoded_fields = '%2C'.join(fields or DEFAULT_UNIPROT_FIELDS)
    return f'{UNIPROT_STREAM_BASE_URL}fields={encoded_fields}&format=tsv&query={query}'


def _parse_uniprot_response(response_text: str) -> pd.DataFrame | None:
    """Parse a UniProt TSV response into a data frame.

    Args:
        response_text: Raw TSV text returned by the UniProt stream endpoint.

    Returns:
        A parsed data frame, or None if the response had no data rows.
    """

    if not response_text.strip():
        return None

    return pd.read_csv(StringIO(response_text), sep='\t')


def _fetch_protein_batch(
    identifiers: list[str],
    fields: list[str] | None,
) -> list[pd.DataFrame]:
    """Fetch one UniProt protein batch, retrying by splitting on HTTP 400.

    The UniProt stream endpoint rejects very long accession queries with
    HTTP 400. A failed batch is recursively split into smaller requests
    rather than failing the whole fetch.

    Args:
        identifiers: UniProt accessions for this batch.
        fields: Return fields. Defaults to DEFAULT_UNIPROT_FIELDS.

    Returns:
        A list of parsed data frames (zero or one, unless the batch was
        split by a retry).
    """

    url = build_uniprot_stream_url(
        build_uniprot_accession_query(identifiers),
        fields=fields,
    )
    response = requests.get(url, timeout=60)

    try:
        response.raise_for_status()
    except requests.HTTPError:
        if response.status_code == 400 and len(identifiers) > 1:
            midpoint = len(identifiers) // 2
            left = _fetch_protein_batch(identifiers[:midpoint], fields)
            right = _fetch_protein_batch(identifiers[midpoint:], fields)
            return left + right
        raise

    frame = _parse_uniprot_response(response.text)
    return [frame] if frame is not None else []


def fetch_protein_table(
    identifiers: list[str],
    fields: list[str] | None = None,
) -> pd.DataFrame:
    """Fetch a UniProt annotation table for a batch of protein accessions.

    Args:
        identifiers: UniProt accession strings.
        fields: Return fields. Defaults to DEFAULT_UNIPROT_FIELDS.

    Returns:
        A data frame with one row per accession, columns per fields.
    """

    frames = []

    for start in range(0, len(identifiers), UNIPROT_BATCH_SIZE):
        batch = identifiers[start : start + UNIPROT_BATCH_SIZE]
        frames.extend(_fetch_protein_batch(batch, fields))

    if not frames:
        return pd.DataFrame(columns=fields or DEFAULT_UNIPROT_FIELDS)

    return pd.concat(frames, ignore_index=True)


def fetch_proteome_table(
    proteome_id: str,
    fields: list[str] | None = None,
) -> pd.DataFrame:
    """Fetch a full proteome annotation table from UniProt.

    Args:
        proteome_id: UniProt proteome identifier (e.g. UP000000625).
        fields: Return fields. Defaults to DEFAULT_UNIPROT_FIELDS.

    Returns:
        A data frame with one row per protein in the proteome, plus a
        Proteome_ID column.
    """

    url = build_uniprot_stream_url(
        f'%28%28proteome%3A{proteome_id}%29%29',
        fields=fields,
    )
    response = requests.get(url, timeout=60)
    response.raise_for_status()

    frame = _parse_uniprot_response(response.text)
    if frame is None:
        frame = pd.DataFrame(columns=fields or DEFAULT_UNIPROT_FIELDS)

    frame['Proteome_ID'] = proteome_id
    return frame
