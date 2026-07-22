"""Reading and parsing FASTA protein sequence files."""

from pathlib import Path
from typing import Union

PathLike = Union[str, Path]


def read_fasta_sequences(filename: PathLike) -> dict[str, str]:
    """Read a FASTA file into a header-to-sequence mapping.

    Args:
        filename: Path to the FASTA file.

    Returns:
        A mapping from FASTA header (without the leading '>') to sequence.
    """

    sequences: dict[str, str] = {}
    current_header: str | None = None
    current_fragments: list[str] = []

    with open(Path(filename), encoding='utf-8') as fasta_file:
        for raw_line in fasta_file:
            line = raw_line.strip()

            if not line:
                continue

            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_fragments)

                current_header = line[1:]
                current_fragments = []
                continue

            current_fragments.append(line)

    if current_header is not None:
        sequences[current_header] = ''.join(current_fragments)

    return sequences


def extract_uniprot_id(fasta_header: str) -> str:
    """Extract the UniProt accession from a FASTA header.

    Args:
        fasta_header: FASTA description line without the leading '>'.

    Returns:
        The UniProt accession parsed from the header.

    Raises:
        ValueError: If the header does not follow UniProt's
            '>db|accession|entry_name ...' structure.
    """

    fields = fasta_header.split('|')
    if len(fields) < 2:
        raise ValueError(f'FASTA header does not contain a UniProt accession: {fasta_header}')
    return fields[1]
