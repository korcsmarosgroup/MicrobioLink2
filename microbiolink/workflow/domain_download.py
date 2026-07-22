"""Downloading Pfam domain annotations for human and/or microbial proteins."""

import pandas as pd

from ..utils import id_resolution
from ..utils import uniprot_client


def domain_table_to_mapping(domain_table: pd.DataFrame) -> dict[str, list[str]]:
    """Convert a UniProt Pfam table into a pfam_id -> uniprot_ids mapping."""

    mapping: dict[str, list[str]] = {}
    for _, row in domain_table.iterrows():
        pfam = row['Pfam']
        if pd.isna(pfam):
            continue
        for pfam_id in (entry for entry in pfam.split(';') if entry):
            mapping.setdefault(pfam_id, []).append(row['Entry'])
    return mapping


def _download_species_domains(identifiers: list[str], id_type: str) -> dict[str, list[str]]:
    """Resolve identifiers to UniProt accessions and fetch their Pfam domains."""

    uniprot_ids = id_resolution.resolve_uniprot_ids(identifiers, id_type)
    domain_table = uniprot_client.fetch_protein_table(uniprot_ids)
    return domain_table_to_mapping(domain_table)


def download_domains(
    human_identifiers: list[str] | None = None,
    human_id_type: str | None = None,
    microbial_identifiers: list[str] | None = None,
    microbial_id_type: str | None = None,
) -> dict[str, dict[str, list[str]]]:
    """Download Pfam domains for human and/or microbial proteins.

    Args:
        human_identifiers: UniProt accessions or gene symbols, or None to
            skip the human side. human_id_type is required if given.
        human_id_type: 'uniprot' or 'genesymbol'.
        microbial_identifiers: UniProt accessions or proteome identifiers,
            or None to skip the microbial side. microbial_id_type is
            required if given.
        microbial_id_type: 'uniprot' or 'proteome'.

    Returns:
        Mapping of 'human'/'microbial' to a pfam_id -> uniprot_ids mapping.
    """

    if human_identifiers is None and microbial_identifiers is None:
        raise ValueError('At least one of human_identifiers or microbial_identifiers is required.')

    species_requests = [
        ('human', human_identifiers, human_id_type),
        ('microbial', microbial_identifiers, microbial_id_type),
    ]

    return {
        species: _download_species_domains(identifiers, id_type)
        for species, identifiers, id_type in species_requests
        if identifiers is not None
    }
