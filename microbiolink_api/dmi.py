#!/usr/bin/env python

"""Domain-motif interaction utilities for the library-style MicrobioLink API."""

from __future__ import annotations

from dataclasses import asdict
from dataclasses import dataclass
import importlib.resources
from pathlib import Path
from typing import Iterable

import pandas as pd
import re

from microbiolink_api.exceptions import InputFormatError
import microbiolink_api.resources


PathLike = str | Path


@dataclass(frozen = True)
class DomainMotifInteraction:
    """Represent one predicted domain-motif interaction."""

    human_protein: str
    motif: str
    start: int
    end: int
    bacterial_domain: str
    bacterial_protein: str


@dataclass(frozen = True)
class DMIResourceBundle:
    """Collect the built-in or user-provided DMI resource tables."""

    elm_regex: dict[str, str]
    motif_domains: dict[str, list[str]]


def extract_uniprot_id(fasta_header: str) -> str:
    """Extract the UniProt accession from a FASTA header.

    Args:
        fasta_header: FASTA description line without the leading `>`.

    Returns:
        The UniProt accession parsed from the header.

    Raises:
        InputFormatError: If the header does not follow the expected UniProt
            structure.
    """

    fields = fasta_header.split('|')
    if len(fields) < 2:
        raise InputFormatError(
            f'FASTA header does not contain a UniProt accession: {fasta_header}',
        )
    return fields[1]


def read_fasta_sequences(filename: PathLike) -> dict[str, str]:
    """Read a FASTA file into a header-to-sequence mapping.

    Args:
        filename: Path to the FASTA file.

    Returns:
        A mapping from FASTA header to sequence string.
    """

    sequences: dict[str, str] = {}
    current_header: str | None = None
    current_fragments: list[str] = []

    with open(filename, encoding = 'utf-8') as fasta_file:
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


def write_fasta_sequences(
    sequences: dict[str, str],
    output_file: PathLike,
) -> None:
    """Write a FASTA mapping to disk.

    Args:
        sequences: Mapping from FASTA header to sequence.
        output_file: Path to the output FASTA file.
    """

    with open(output_file, 'w', encoding = 'utf-8') as fasta_file:
        for header, sequence in sequences.items():
            fasta_file.write(f'>{header}\n')
            fasta_file.write(f'{sequence}\n')


def _parse_elm_regex_lines(lines: Iterable[str]) -> dict[str, str]:
    """Parse ELM motif regex data from an iterable of lines."""

    elm_regex: dict[str, str] = {}

    iterator = iter(lines)
    next(iterator, None)

    for line in iterator:
        if not line or line.startswith('#'):
            continue

        fields = line.replace('"', '').strip().split('\t')
        if len(fields) > 4:
            elm_regex[fields[1]] = fields[4]

    return elm_regex


def _parse_motif_domain_lines(lines: Iterable[str]) -> dict[str, list[str]]:
    """Parse motif-domain data from an iterable of lines."""

    motif_domains: dict[str, list[str]] = {}

    iterator = iter(lines)
    next(iterator, None)

    for line in iterator:
        fields = line.replace('"', '').strip().split('\t')
        if len(fields) <= 1:
            continue

        motif_domains.setdefault(fields[0], []).append(fields[1])

    return motif_domains


def read_elm_regex_table(filename: PathLike) -> dict[str, str]:
    """Read an ELM motif regex table from disk.

    Args:
        filename: Path to the ELM table.

    Returns:
        A mapping from motif identifier to regular expression.
    """

    with open(filename, encoding = 'utf-8') as motif_table:
        return _parse_elm_regex_lines(motif_table)


def read_motif_domain_table(filename: PathLike) -> dict[str, list[str]]:
    """Read motif-domain relationships from a table on disk.

    Args:
        filename: Path to the motif-domain table.

    Returns:
        A mapping from motif identifier to a list of Pfam domains.
    """

    with open(filename, encoding = 'utf-8') as motif_domain_table:
        return _parse_motif_domain_lines(motif_domain_table)


def load_default_dmi_resource_bundle() -> DMIResourceBundle:
    """Load the packaged default DMI resource tables.

    Returns:
        The packaged ELM regex and motif-domain tables.
    """

    elm_regex_path = importlib.resources.files(
        microbiolink_api.resources,
    ).joinpath('elm_classes.tsv')
    motif_domain_path = importlib.resources.files(
        microbiolink_api.resources,
    ).joinpath('elm_interaction_domains.tsv')

    return DMIResourceBundle(
        elm_regex = read_elm_regex_table(str(elm_regex_path)),
        motif_domains = read_motif_domain_table(str(motif_domain_path)),
    )


def _resolve_dmi_resource_bundle(
    elm_regex_file: PathLike | None,
    motif_domain_file: PathLike | None,
    resource_bundle: DMIResourceBundle | None,
) -> DMIResourceBundle:
    """Resolve built-in and user-provided DMI resource tables."""

    resolved_bundle = resource_bundle or load_default_dmi_resource_bundle()

    return DMIResourceBundle(
        elm_regex = (
            read_elm_regex_table(elm_regex_file)
            if elm_regex_file is not None
            else resolved_bundle.elm_regex
        ),
        motif_domains = (
            read_motif_domain_table(motif_domain_file)
            if motif_domain_file is not None
            else resolved_bundle.motif_domains
        ),
    )


def read_bacterial_domain_table(filename: PathLike) -> dict[str, list[str]]:
    """Read bacterial proteins grouped by Pfam domain.

    Args:
        filename: Path to the bacterial domain table.

    Returns:
        A mapping from Pfam domain to bacterial proteins containing that domain.
    """

    bacterial_domains: dict[str, list[str]] = {}

    with open(filename, encoding = 'utf-8') as protein_domain_file:
        next(protein_domain_file, None)

        for line in protein_domain_file:
            fields = line.strip().split('\t')
            if len(fields) <= 1:
                continue

            for pfam_domain in fields[1].split(';'):
                bacterial_domains.setdefault(pfam_domain, []).append(fields[0])

    return bacterial_domains


def select_sequences_by_uniprot_ids(
    sequences: dict[str, str],
    uniprot_ids: list[str],
) -> dict[str, str]:
    """Select FASTA sequences whose UniProt accessions are in a target list.

    Args:
        sequences: Mapping from FASTA header to sequence.
        uniprot_ids: UniProt accessions to retain.

    Returns:
        A filtered mapping containing only the requested accessions.
    """

    wanted_ids = set(uniprot_ids)
    return {
        header: sequence
        for header, sequence in sequences.items()
        if extract_uniprot_id(header) in wanted_ids
    }


def _find_motif_matches(
    sequences: dict[str, str],
    elm_regex: dict[str, str],
) -> dict[str, list[tuple[str, int, int]]]:
    """Find motif matches for each human protein."""

    motif_matches: dict[str, list[tuple[str, int, int]]] = {}

    for header, sequence in sequences.items():
        uniprot_id = extract_uniprot_id(header)
        matches: list[tuple[str, int, int]] = []

        for motif_name, motif_pattern in elm_regex.items():
            for match in re.finditer(motif_pattern, sequence):
                matches.append((motif_name, match.start(), match.end()))

        if matches:
            motif_matches[uniprot_id] = matches

    return motif_matches


def predict_domain_motif_interactions_from_data(
    human_sequences: dict[str, str],
    elm_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    bacterial_domains: dict[str, list[str]],
) -> list[DomainMotifInteraction]:
    """Predict domain-motif interactions from in-memory inputs.

    Args:
        human_sequences: Mapping from FASTA header to human protein sequence.
        elm_regex: Mapping from motif identifier to regular expression.
        motif_domains: Mapping from motif identifier to compatible Pfam domains.
        bacterial_domains: Mapping from Pfam domain to bacterial proteins.

    Returns:
        A list of predicted interactions.
    """

    motif_matches = _find_motif_matches(human_sequences, elm_regex)
    interactions: list[DomainMotifInteraction] = []

    for motif_name, compatible_domains in motif_domains.items():
        motif_hits = [
            (human_protein, start, end)
            for human_protein, matches in motif_matches.items()
            for match_name, start, end in matches
            if match_name == motif_name
        ]

        for domain_name in compatible_domains:
            for bacterial_protein in bacterial_domains.get(domain_name, []):
                for human_protein, start, end in motif_hits:
                    interactions.append(
                        DomainMotifInteraction(
                            human_protein = human_protein,
                            motif = motif_name,
                            start = start,
                            end = end,
                            bacterial_domain = domain_name,
                            bacterial_protein = bacterial_protein,
                        ),
                    )

    return interactions


def predict_domain_motif_interactions(
    fasta_file: PathLike,
    bacterial_domain_file: PathLike,
    elm_regex_file: PathLike | None = None,
    motif_domain_file: PathLike | None = None,
    resource_bundle: DMIResourceBundle | None = None,
) -> list[DomainMotifInteraction]:
    """Predict domain-motif interactions from input files.

    Args:
        fasta_file: Path to the human FASTA file.
        bacterial_domain_file: Path to the bacterial domain table.
        elm_regex_file: Optional override for the ELM regex table. When omitted,
            the built-in default table is used.
        motif_domain_file: Optional override for the motif-domain table. When
            omitted, the built-in default table is used.
        resource_bundle: Optional in-memory resource bundle. This is useful when
            you want to override the built-in defaults without reading from
            extra files.

    Returns:
        A list of predicted interactions.
    """

    human_sequences = read_fasta_sequences(fasta_file)
    resolved_resources = _resolve_dmi_resource_bundle(
        elm_regex_file = elm_regex_file,
        motif_domain_file = motif_domain_file,
        resource_bundle = resource_bundle,
    )
    bacterial_domains = read_bacterial_domain_table(bacterial_domain_file)

    return predict_domain_motif_interactions_from_data(
        human_sequences = human_sequences,
        elm_regex = resolved_resources.elm_regex,
        motif_domains = resolved_resources.motif_domains,
        bacterial_domains = bacterial_domains,
    )


def interactions_to_dataframe(
    interactions: list[DomainMotifInteraction],
) -> pd.DataFrame:
    """Convert interaction records to a data frame.

    Args:
        interactions: Predicted interaction records.

    Returns:
        A tabular representation of the interactions.
    """

    columns = [
        'human_protein',
        'motif',
        'start',
        'end',
        'bacterial_domain',
        'bacterial_protein',
    ]

    return pd.DataFrame(
        [asdict(interaction) for interaction in interactions],
        columns = columns,
    )


def write_domain_motif_interactions(
    interactions: list[DomainMotifInteraction],
    output_file: PathLike,
    separator: str = ';',
) -> None:
    """Write predicted interactions to disk.

    Args:
        interactions: Predicted interaction records.
        output_file: Path to the output file.
        separator: Field separator for the output table.
    """

    interaction_frame = interactions_to_dataframe(interactions)
    interaction_frame.to_csv(
        output_file,
        sep = separator,
        index = False,
    )
