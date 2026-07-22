"""Shared low-level UniProt REST client for identifier and proteome lookups."""

from io import StringIO
from pathlib import Path
from typing import Union

import pandas as pd
import requests

PathLike = Union[str, Path]

UNIPROT_STREAM_BASE_URL = 'https://rest.uniprot.org/uniprotkb/stream?'
DEFAULT_UNIPROT_FIELDS = ['accession', 'xref_pfam', 'gene_names']
UNIPROT_BATCH_SIZE = 100
FASTA_BATCH_SIZE = 100


def read_ids(
    filename: PathLike,
    separator: str,
    id_column: int,
    has_header: bool = True,
) -> list[str]:
    """Read UniProt or proteome identifiers from a delimited file.

    Args:
        filename: Path to the identifier file.
        separator: Field separator used in the file.
        id_column: One-based column number containing the identifiers.
        has_header: Whether the first line is a header to skip. Defaults to
            True; set False for a file with no header row.

    Returns:
        A list of identifiers in file order.
    """

    column_index = id_column - 1
    ids = []

    with open(Path(filename), encoding='utf-8-sig') as id_list:
        if has_header:
            next(id_list, None)

        for line_number, line in enumerate(id_list, start=2 if has_header else 1):
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
    format: str = 'tsv',
    fields: list[str] | None = None,
) -> str:
    """Build a UniProt stream endpoint URL.

    Args:
        query: Percent-encoded UniProt query string.
        format: Response format, 'tsv' or 'fasta'. FASTA doesn't support
            field selection, so `fields` is ignored when format='fasta'.
        fields: Return fields for 'tsv' format. Defaults to
            DEFAULT_UNIPROT_FIELDS.

    Returns:
        A complete UniProt stream URL.
    """

    if format == 'fasta':
        return f'{UNIPROT_STREAM_BASE_URL}format=fasta&query={query}'

    encoded_fields = '%2C'.join(fields or DEFAULT_UNIPROT_FIELDS)
    return f'{UNIPROT_STREAM_BASE_URL}fields={encoded_fields}&format=tsv&query={query}'


def _get(url: str) -> requests.Response:
    """GET a UniProt REST URL, raising on any HTTP error."""

    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response


def _chunked(items: list[str], size: int):
    """Yield successive size-length chunks of items."""

    for start in range(0, len(items), size):
        yield items[start : start + size]


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

    for batch in _chunked(identifiers, UNIPROT_BATCH_SIZE):
        url = build_uniprot_stream_url(build_uniprot_accession_query(batch), fields=fields)
        frame = _parse_uniprot_response(_get(url).text)
        if frame is not None:
            frames.append(frame)

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
    response = _get(url)

    frame = _parse_uniprot_response(response.text)
    if frame is None:
        frame = pd.DataFrame(columns=fields or DEFAULT_UNIPROT_FIELDS)

    frame['Proteome_ID'] = proteome_id
    return frame


def fetch_fasta_sequences(identifiers: list[str]) -> str:
    """Fetch FASTA sequence text for a batch of UniProt accessions.

    Args:
        identifiers: UniProt accession strings.

    Returns:
        Concatenated raw FASTA text for all matched accessions.
    """

    texts = []
    for batch in _chunked(identifiers, FASTA_BATCH_SIZE):
        url = build_uniprot_stream_url(build_uniprot_accession_query(batch), format='fasta')
        texts.append(_get(url).text)

    return ''.join(texts)
