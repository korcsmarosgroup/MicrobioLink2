"""Downloading human and/or microbial protein sequences as FASTA files."""

from pathlib import Path
from typing import Union

from ..utils import id_resolution
from ..utils import uniprot_client

PathLike = Union[str, Path]

HUMAN_FASTA_FILENAME = 'human_proteins.fasta'
MICROBIAL_FASTA_FILENAME = 'microbial_proteins.fasta'


def _download_species_fasta(
    identifiers: list[str],
    id_type: str,
    output_path: Path,
) -> Path:
    """Resolve identifiers to UniProt accessions and write their FASTA to disk."""

    uniprot_ids = id_resolution.resolve_uniprot_ids(identifiers, id_type)
    fasta_text = uniprot_client.fetch_fasta_sequences(uniprot_ids)
    output_path.write_text(fasta_text)
    return output_path


def download_fasta(
    output_folder: PathLike,
    human_identifiers: list[str] | None = None,
    human_id_type: str | None = None,
    microbial_identifiers: list[str] | None = None,
    microbial_id_type: str | None = None,
) -> dict[str, Path]:
    """Download human and/or microbial protein sequences as FASTA files.

    Args:
        output_folder: Directory to write the FASTA file(s) into.
        human_identifiers: UniProt accessions or gene symbols, or None to
            skip the human side. human_id_type is required if given.
        human_id_type: 'uniprot' or 'genesymbol'.
        microbial_identifiers: UniProt accessions or proteome identifiers,
            or None to skip the microbial side. microbial_id_type is
            required if given.
        microbial_id_type: 'uniprot' or 'proteome'.

    Returns:
        Mapping of 'human'/'microbial' to the FASTA file path written.
    """

    if human_identifiers is None and microbial_identifiers is None:
        raise ValueError('At least one of human_identifiers or microbial_identifiers is required.')

    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    species_requests = [
        ('human', human_identifiers, human_id_type, HUMAN_FASTA_FILENAME),
        ('microbial', microbial_identifiers, microbial_id_type, MICROBIAL_FASTA_FILENAME),
    ]

    return {
        species: _download_species_fasta(identifiers, id_type, output_folder / filename)
        for species, identifiers, id_type, filename in species_requests
        if identifiers is not None
    }
