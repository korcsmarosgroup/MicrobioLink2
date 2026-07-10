#!/usr/bin/env python

"""Microbiome-side identifier and bacterial domain utilities."""

from __future__ import annotations

from io import StringIO
from pathlib import Path

import pandas as pd
import requests

from microbiolink.download_bacterial_proteins import download_protein_list_with_fields
from microbiolink.download_bacterial_proteins import download_proteome_with_fields
from microbiolink.download_bacterial_proteins import read_ids


PathLike = str | Path
LOCATION_FIELD = 'cc_subcellular_location'
LOCATION_COLUMN = 'Subcellular location [CC]'
UNIPROT_PROTEIN_BATCH_SIZE = 100


def read_microbiome_identifiers(
    filename: PathLike,
    separator: str,
    id_column: int,
) -> list[str]:
    """Read microbiome-side identifiers from a file.

    Args:
        filename: Path to a delimited file containing UniProt protein IDs or
            UniProt proteome IDs.
        separator: Field separator used in the input file.
        id_column: One-based column number containing the identifiers.

    Returns:
        A list of identifiers in file order.
    """

    return read_ids(filename, separator, id_column)


def fetch_bacterial_domain_table_from_ids(
    identifiers: list[str],
    id_type: str,
    include_location: bool = False,
) -> pd.DataFrame:
    """Fetch bacterial protein domain annotations from UniProt.

    Args:
        identifiers: UniProt protein accessions or UniProt proteome IDs.
        id_type: Identifier type. Use `'Uniprot'` for protein accessions or
            `'UP'` for proteome IDs.
        include_location: Whether to request the UniProt subcellular location
            annotation column.

    Returns:
        A data frame containing the downloaded UniProt table.
    """

    frames: list[pd.DataFrame] = []
    fields = None

    if include_location:
        fields = [
            'accession',
            'xref_pfam',
            'gene_names',
            LOCATION_FIELD,
        ]

    if id_type == 'UP':
        for proteome_id in identifiers:
            response_text = download_proteome_with_fields(
                proteome_id,
                fields = fields,
            )
            frame = pd.read_csv(
                StringIO(response_text),
                sep = '\t',
            )
            frame['Proteome_ID'] = proteome_id
            frames.append(frame)
    else:
        frames.extend(
            _download_protein_identifier_frames(
                identifiers,
                fields = fields,
            ),
        )

    if not frames:
        return pd.DataFrame(columns = ['Entry', 'Pfam', 'Gene Names'])

    return pd.concat(frames, ignore_index = True)


def _download_protein_identifier_frames(
    identifiers: list[str],
    fields: list[str] | None,
) -> list[pd.DataFrame]:
    """Download UniProt protein annotation tables in request-safe batches.

    The UniProt stream endpoint rejects very long accession queries with HTTP
    400. This helper batches accessions proactively and, if needed, retries a
    failed batch by recursively splitting it into smaller requests.
    """

    frames: list[pd.DataFrame] = []

    for start in range(0, len(identifiers), UNIPROT_PROTEIN_BATCH_SIZE):
        batch = identifiers[start : start + UNIPROT_PROTEIN_BATCH_SIZE]
        frames.extend(
            _download_protein_identifier_batch(
                batch,
                fields = fields,
            ),
        )

    return frames


def _download_protein_identifier_batch(
    identifiers: list[str],
    fields: list[str] | None,
) -> list[pd.DataFrame]:
    """Download one UniProt protein batch with retry-by-splitting fallback."""

    try:
        response_text = download_protein_list_with_fields(
            identifiers,
            fields = fields,
        )
    except requests.HTTPError as error:
        response = getattr(error, 'response', None)
        if (
            response is not None
            and response.status_code == 400
            and len(identifiers) > 1
        ):
            midpoint = len(identifiers) // 2
            left_batch = _download_protein_identifier_batch(
                identifiers[:midpoint],
                fields = fields,
            )
            right_batch = _download_protein_identifier_batch(
                identifiers[midpoint:],
                fields = fields,
            )
            return left_batch + right_batch
        raise

    if not response_text.strip():
        return []

    frame = pd.read_csv(
        StringIO(response_text),
        sep = '\t',
    )
    return [frame]


def fetch_bacterial_domain_table_from_file(
    id_file: PathLike,
    id_type: str,
    separator: str,
    id_column: int,
    include_location: bool = False,
    output_file: PathLike | None = None,
) -> pd.DataFrame:
    """Read microbiome identifiers from a file and fetch domain annotations.

    Args:
        id_file: File containing UniProt protein IDs or proteome IDs.
        id_type: Identifier type. Use `'Uniprot'` for proteins or `'UP'` for
            proteomes.
        separator: Field separator used in the input file.
        id_column: One-based column containing the identifiers.
        include_location: Whether to request the UniProt subcellular location
            annotation column.
        output_file: Optional output path for the downloaded domain table.

    Returns:
        A data frame containing the UniProt domain table.
    """

    identifiers = read_microbiome_identifiers(
        id_file,
        separator = separator,
        id_column = id_column,
    )
    domain_table = fetch_bacterial_domain_table_from_ids(
        identifiers,
        id_type = id_type,
        include_location = include_location,
    )

    if output_file is not None:
        domain_table.to_csv(output_file, sep = '\t', index = False)

    return domain_table


def filter_bacterial_domain_table_by_location(
    domain_table: pd.DataFrame,
    location_filters: list[str],
) -> pd.DataFrame:
    """Filter a bacterial domain table by subcellular-location annotations.

    Args:
        domain_table: UniProt download result.
        location_filters: Case-insensitive substrings to match against the
            `Subcellular location [CC]` column.

    Returns:
        The filtered domain table.
    """

    if not location_filters:
        return domain_table

    if LOCATION_COLUMN not in domain_table.columns:
        raise ValueError(
            'Subcellular location annotations are not present in the domain '
            'table. Fetch the table with `include_location = True` first.',
        )

    normalized_filters = [location.lower() for location in location_filters]
    location_series = domain_table[LOCATION_COLUMN].fillna('').astype(str).str.lower()

    mask = location_series.apply(
        lambda value: any(location in value for location in normalized_filters),
    )
    return domain_table.loc[mask].reset_index(drop = True)


def bacterial_domain_dataframe_to_mapping(
    domain_table: pd.DataFrame,
) -> dict[str, list[str]]:
    """Convert a UniProt domain table to the DMI domain mapping format.

    Args:
        domain_table: UniProt download result with at least `Entry` and `Pfam`
            columns.

    Returns:
        A mapping from Pfam domain to bacterial proteins containing that domain.
    """

    if domain_table.empty:
        return {}

    mapping: dict[str, list[str]] = {}

    for _, row in domain_table.iterrows():
        protein = str(row['Entry'])
        pfam_field = row.get('Pfam')

        if pd.isna(pfam_field):
            continue

        for pfam_domain in str(pfam_field).split(';'):
            pfam_domain = pfam_domain.strip()
            if pfam_domain:
                mapping.setdefault(pfam_domain, []).append(protein)

    return mapping
