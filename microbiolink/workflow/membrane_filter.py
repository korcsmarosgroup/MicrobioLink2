"""Filtering of human and bacterial proteins to membrane/secreted proteins."""

import pandas as pd

from ..utils import uniprot_client

LOCATION_FIELD = 'cc_subcellular_location'
LOCATION_COLUMN = 'Subcellular location [CC]'

OUTPUT_COLUMNS = ['uniprot_id', 'location_annotation']


def filter_human_membrane_proteins(
    identifiers: list[str],
    id_type: str,
    location_filters: list[str],
) -> pd.DataFrame:
    """Filter human proteins to membrane/secreted proteins via OmniPath Intercell.

    Args:
        identifiers: UniProt accessions or human gene symbols.
        id_type: Identifier type of `identifiers`: 'uniprot' or 'genesymbol'.
        location_filters: OmniPath Intercell `parent` categories to keep
            (e.g. 'plasma_membrane_transmembrane', 'secreted').

    Returns:
        A data frame with columns uniprot_id and location_annotation, one
        row per input protein present in the Intercell table under one of
        the requested categories.
    """

    if id_type not in ('uniprot', 'genesymbol'):
        raise ValueError(f"id_type must be 'uniprot' or 'genesymbol' for human species, got {id_type!r}")

    import omnipath as op

    intercell_table = op.requests.Intercell.get(
        parent=location_filters,
        scope=['generic', 'specific'],
        source=['resource_specific', 'composite'],
        entity_type='protein',
    )

    match_column = 'uniprot' if id_type == 'uniprot' else 'genesymbol'
    matches = intercell_table[intercell_table[match_column].isin(identifiers)]

    grouped = matches.groupby('uniprot')['parent'].agg(lambda values: ';'.join(sorted(set(values))))

    return grouped.reset_index().rename(
        columns={'uniprot': 'uniprot_id', 'parent': 'location_annotation'},
    )


def filter_bacterial_membrane_proteins(
    identifiers: list[str],
    id_type: str,
    location_filters: list[str],
) -> pd.DataFrame:
    """Filter bacterial proteins to membrane/secreted proteins via UniProt.

    Args:
        identifiers: UniProt accessions, or UniProt proteome identifiers
            when id_type is 'proteome'.
        id_type: Identifier type of `identifiers`: 'uniprot' or 'proteome'.
        location_filters: Case-insensitive substrings to match against
            UniProt's Subcellular location [CC] annotation (e.g.
            'outer membrane', 'plasma membrane').

    Returns:
        A data frame with columns uniprot_id and location_annotation, one
        row per protein whose location annotation matched at least one of
        the requested substrings.
    """

    domain_table = _fetch_bacterial_location_table(identifiers, id_type)
    domain_table = domain_table.rename(
        columns={'Entry': 'uniprot_id', LOCATION_COLUMN: 'location_annotation'},
    )

    if domain_table.empty or 'location_annotation' not in domain_table.columns:
        return pd.DataFrame(columns=OUTPUT_COLUMNS)

    normalized_filters = [location.lower() for location in location_filters]
    location_series = domain_table['location_annotation'].fillna('').astype(str).str.lower()
    mask = location_series.apply(lambda value: any(term in value for term in normalized_filters))

    return domain_table.loc[mask, OUTPUT_COLUMNS].reset_index(drop=True)


def _fetch_bacterial_location_table(
    identifiers: list[str],
    id_type: str,
) -> pd.DataFrame:
    """Fetch a UniProt table with accession and location annotation columns."""

    fields = ['accession', LOCATION_FIELD]

    if id_type == 'uniprot':
        return uniprot_client.fetch_protein_table(identifiers, fields=fields)

    if id_type == 'proteome':
        frames = [
            uniprot_client.fetch_proteome_table(proteome_id, fields=fields)
            for proteome_id in identifiers
        ]
        if not frames:
            return pd.DataFrame(columns=fields)
        return pd.concat(frames, ignore_index=True)

    raise ValueError(f"id_type must be 'uniprot' or 'proteome' for microbial species, got {id_type!r}")


def filter_membrane_proteins(
    identifiers: list[str],
    id_type: str,
    species: str,
    location_filters: list[str],
) -> pd.DataFrame:
    """Filter protein identifiers to membrane/secreted proteins.

    Public dispatcher for Module 2: routes to the OmniPath Intercell-based
    filter for human proteins or the UniProt location-annotation filter for
    bacterial proteins.

    Args:
        identifiers: UniProt accessions, human gene symbols, or UniProt
            proteome identifiers.
        id_type: Identifier type of `identifiers`. 'uniprot' or 'genesymbol'
            for species='human'; 'uniprot' or 'proteome' for
            species='microbial'.
        species: 'human' or 'microbial'.
        location_filters: Location categories/substrings to keep, meaning
            depends on species (see filter_human_membrane_proteins and
            filter_bacterial_membrane_proteins).

    Returns:
        A data frame with columns uniprot_id and location_annotation.
    """

    if species == 'human':
        return filter_human_membrane_proteins(identifiers, id_type, location_filters)

    if species == 'microbial':
        return filter_bacterial_membrane_proteins(identifiers, id_type, location_filters)

    raise ValueError(f"species must be 'human' or 'microbial', got {species!r}")
