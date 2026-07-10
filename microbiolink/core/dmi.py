#!/usr/bin/env python

"""Domain-motif interaction utilities for the library-style MicrobioLink API."""

from __future__ import annotations

from dataclasses import asdict
from dataclasses import dataclass
from dataclasses import field
import importlib.resources
from pathlib import Path
import re
from typing import Iterable

import pandas as pd

from microbiolink.core.exceptions import InputFormatError
import microbiolink.core.resources


PathLike = str | Path
REVERSE_EXCLUDED_MOTIF_PREFIXES = ('CLV_',)


@dataclass(frozen = True)
class DomainMotifInteraction:
    """Represent one forward host-microbe domain-motif interaction."""

    human_protein: str
    motif: str
    start: int
    end: int
    bacterial_domain: str
    bacterial_protein: str
    resource: str = 'ELM'


@dataclass(frozen = True)
class BidirectionalDomainMotifInteraction:
    """Represent one directional host-microbe domain-motif interaction."""

    host_protein: str
    microbial_protein: str
    motif: str
    start: int
    end: int
    domain: str
    motif_protein_side: str
    domain_protein_side: str
    resource: str = 'ELM'


@dataclass(frozen = True)
class DMIResourceBundle:
    """Collect the built-in or user-provided DMI resource tables."""

    elm_regex: dict[str, str]
    motif_domains: dict[str, list[str]]
    motif_sources: dict[str, str] = field(default_factory = dict)


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
    """Parse motif regex data from an iterable of lines."""

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


def _unique_preserve_order(values: Iterable[str]) -> list[str]:
    """Return unique values while preserving input order."""

    return list(dict.fromkeys(values))


def _is_reverse_compatible_motif(motif_name: str) -> bool:
    """Return whether a motif should be retained in reverse-mode DMI."""

    return not motif_name.startswith(REVERSE_EXCLUDED_MOTIF_PREFIXES)


def _filter_reverse_motif_resources(
    motif_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    motif_sources: dict[str, str] | None,
) -> tuple[dict[str, str], dict[str, list[str]], dict[str, str] | None]:
    """Filter reverse-mode motif resources to exclude cleavage motifs."""

    filtered_regex = {
        motif_name: motif_pattern
        for motif_name, motif_pattern in motif_regex.items()
        if _is_reverse_compatible_motif(motif_name)
    }
    filtered_domains = {
        motif_name: compatible_domains
        for motif_name, compatible_domains in motif_domains.items()
        if _is_reverse_compatible_motif(motif_name)
    }

    if motif_sources is None:
        return filtered_regex, filtered_domains, None

    filtered_sources = {
        motif_name: source_name
        for motif_name, source_name in motif_sources.items()
        if _is_reverse_compatible_motif(motif_name)
    }

    return filtered_regex, filtered_domains, filtered_sources


def read_elm_regex_table(filename: PathLike) -> dict[str, str]:
    """Read a motif regex table from disk.

    Args:
        filename: Path to the motif table.

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


def merge_dmi_resource_bundles(
    *resource_bundles: DMIResourceBundle,
) -> DMIResourceBundle:
    """Merge multiple DMI resource bundles.

    Later bundles may add motif-domain relationships to existing motifs, but
    cannot redefine a motif regular expression.
    """

    merged_regex: dict[str, str] = {}
    merged_domains: dict[str, list[str]] = {}
    merged_sources: dict[str, str] = {}

    for bundle in resource_bundles:
        for motif_name, motif_pattern in bundle.elm_regex.items():
            existing_pattern = merged_regex.get(motif_name)

            if existing_pattern is not None and existing_pattern != motif_pattern:
                raise ValueError(
                    f'Motif {motif_name} has conflicting patterns: '
                    f'{existing_pattern} vs {motif_pattern}',
                )

            merged_regex[motif_name] = motif_pattern

        for motif_name, domains in bundle.motif_domains.items():
            merged_domains[motif_name] = _unique_preserve_order(
                [*merged_domains.get(motif_name, []), *domains],
            )

        for motif_name, source_name in bundle.motif_sources.items():
            existing_source = merged_sources.get(motif_name)

            if existing_source is not None and existing_source != source_name:
                raise ValueError(
                    f'Motif {motif_name} has conflicting sources: '
                    f'{existing_source} vs {source_name}',
                )

            merged_sources[motif_name] = source_name

    for motif_name in set(merged_regex) | set(merged_domains):
        merged_sources.setdefault(motif_name, 'custom')

    return DMIResourceBundle(
        elm_regex = merged_regex,
        motif_domains = merged_domains,
        motif_sources = merged_sources,
    )


def load_default_elm_dmi_resource_bundle() -> DMIResourceBundle:
    """Load the packaged ELM DMI resource tables."""

    elm_regex_path = importlib.resources.files(
        microbiolink.core.resources,
    ).joinpath('elm_classes.tsv')
    motif_domain_path = importlib.resources.files(
        microbiolink.core.resources,
    ).joinpath('elm_interaction_domains.tsv')
    elm_regex = read_elm_regex_table(str(elm_regex_path))
    motif_domains = read_motif_domain_table(str(motif_domain_path))

    return DMIResourceBundle(
        elm_regex = elm_regex,
        motif_domains = motif_domains,
        motif_sources = {
            motif_name: 'ELM'
            for motif_name in set(elm_regex) | set(motif_domains)
        },
    )


def load_default_3did_dmi_resource_bundle() -> DMIResourceBundle:
    """Load the packaged 3did structural DMI resource tables."""

    motif_regex_path = importlib.resources.files(
        microbiolink.core.resources,
    ).joinpath('3did_dmi_classes.tsv')
    motif_domain_path = importlib.resources.files(
        microbiolink.core.resources,
    ).joinpath('3did_dmi_interaction_domains.tsv')
    motif_regex = read_elm_regex_table(str(motif_regex_path))
    motif_domains = read_motif_domain_table(str(motif_domain_path))

    return DMIResourceBundle(
        elm_regex = motif_regex,
        motif_domains = motif_domains,
        motif_sources = {
            motif_name: '3did'
            for motif_name in set(motif_regex) | set(motif_domains)
        },
    )


def load_default_dmi_resource_bundle() -> DMIResourceBundle:
    """Load the packaged default DMI resource tables.

    Returns:
        The packaged ELM and 3did structural motif-domain tables merged into
        one bundle.
    """

    return merge_dmi_resource_bundles(
        load_default_elm_dmi_resource_bundle(),
        load_default_3did_dmi_resource_bundle(),
    )


def resolve_dmi_resource_bundle_by_name(resource_set: str) -> DMIResourceBundle:
    """Resolve a packaged DMI resource bundle by name.

    Args:
        resource_set: One of `'default'` (ELM and 3did merged), `'elm'`
            (ELM only), or `'3did'` (3did structural motifs only).

    Returns:
        The packaged resource bundle selected by `resource_set`.

    Raises:
        InputFormatError: If `resource_set` is not one of the supported
            names.
    """

    if resource_set == 'elm':
        return load_default_elm_dmi_resource_bundle()

    if resource_set == '3did':
        return load_default_3did_dmi_resource_bundle()

    if resource_set == 'default':
        return load_default_dmi_resource_bundle()

    raise InputFormatError(
        f"`resource_set` must be one of: default, elm, 3did. Got: {resource_set!r}",
    )


def _resolve_dmi_resource_bundle(
    elm_regex_file: PathLike | None,
    motif_domain_file: PathLike | None,
    resource_bundle: DMIResourceBundle | None,
) -> DMIResourceBundle:
    """Resolve built-in and user-provided DMI resource tables."""

    resolved_bundle = resource_bundle or load_default_dmi_resource_bundle()
    elm_regex = (
        read_elm_regex_table(elm_regex_file)
        if elm_regex_file is not None
        else resolved_bundle.elm_regex
    )
    motif_domains = (
        read_motif_domain_table(motif_domain_file)
        if motif_domain_file is not None
        else resolved_bundle.motif_domains
    )

    if elm_regex_file is not None or motif_domain_file is not None:
        motif_sources = {
            motif_name: 'custom'
            for motif_name in set(elm_regex) | set(motif_domains)
        }
    else:
        motif_sources = {
            motif_name: resolved_bundle.motif_sources.get(motif_name, 'custom')
            for motif_name in set(elm_regex) | set(motif_domains)
        }

    return DMIResourceBundle(
        elm_regex = elm_regex,
        motif_domains = motif_domains,
        motif_sources = motif_sources,
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


# Identical logic to `read_bacterial_domain_table`; DDI uses this name since
# it reads domain tables for both bacterial and human proteins.
read_protein_domain_table = read_bacterial_domain_table


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
    motif_regex: dict[str, str],
) -> dict[str, list[tuple[str, int, int]]]:
    """Find motif matches for each protein sequence."""

    motif_matches: dict[str, list[tuple[str, int, int]]] = {}

    for header, sequence in sequences.items():
        uniprot_id = extract_uniprot_id(header)
        matches: list[tuple[str, int, int]] = []

        for motif_name, motif_pattern in motif_regex.items():
            for match in re.finditer(motif_pattern, sequence):
                matches.append((motif_name, match.start(), match.end()))

        if matches:
            motif_matches[uniprot_id] = matches

    return motif_matches


def _predict_directional_domain_motif_interactions(
    motif_sequences: dict[str, str],
    motif_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    partner_domains: dict[str, list[str]],
    motif_sources: dict[str, str] | None,
    motif_protein_side: str,
) -> list[BidirectionalDomainMotifInteraction]:
    """Predict one directional set of domain-motif interactions."""

    motif_matches = _find_motif_matches(motif_sequences, motif_regex)
    interactions: list[BidirectionalDomainMotifInteraction] = []
    resolved_sources = motif_sources or {}

    for motif_name, compatible_domains in motif_domains.items():
        motif_hits = [
            (protein_id, start, end)
            for protein_id, matches in motif_matches.items()
            for match_name, start, end in matches
            if match_name == motif_name
        ]

        if not motif_hits:
            continue

        for domain_name in compatible_domains:
            for partner_protein in partner_domains.get(domain_name, []):
                for motif_protein, start, end in motif_hits:
                    if motif_protein_side == 'host':
                        host_protein = motif_protein
                        microbial_protein = partner_protein
                        domain_protein_side = 'microbe'
                    else:
                        host_protein = partner_protein
                        microbial_protein = motif_protein
                        domain_protein_side = 'host'

                    interactions.append(
                        BidirectionalDomainMotifInteraction(
                            host_protein = host_protein,
                            microbial_protein = microbial_protein,
                            motif = motif_name,
                            start = start,
                            end = end,
                            domain = domain_name,
                            motif_protein_side = motif_protein_side,
                            domain_protein_side = domain_protein_side,
                            resource = resolved_sources.get(motif_name, 'custom'),
                        ),
                    )

    return interactions


def predict_domain_motif_interactions_from_data(
    human_sequences: dict[str, str],
    elm_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    bacterial_domains: dict[str, list[str]],
    motif_sources: dict[str, str] | None = None,
) -> list[DomainMotifInteraction]:
    """Predict forward domain-motif interactions from in-memory inputs.

    Args:
        human_sequences: Mapping from human FASTA header to sequence.
        elm_regex: Mapping from motif identifier to regular expression.
        motif_domains: Mapping from motif identifier to compatible Pfam domains.
        bacterial_domains: Mapping from Pfam domain to bacterial proteins.
        motif_sources: Optional mapping from motif identifier to resource name.

    Returns:
        A list of forward host-microbe interactions.
    """

    bidirectional_interactions = _predict_directional_domain_motif_interactions(
        motif_sequences = human_sequences,
        motif_regex = elm_regex,
        motif_domains = motif_domains,
        partner_domains = bacterial_domains,
        motif_sources = motif_sources,
        motif_protein_side = 'host',
    )

    return [
        DomainMotifInteraction(
            human_protein = interaction.host_protein,
            motif = interaction.motif,
            start = interaction.start,
            end = interaction.end,
            bacterial_domain = interaction.domain,
            bacterial_protein = interaction.microbial_protein,
            resource = interaction.resource,
        )
        for interaction in bidirectional_interactions
    ]


def predict_reverse_domain_motif_interactions_from_data(
    bacterial_sequences: dict[str, str],
    elm_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    human_domains: dict[str, list[str]],
    motif_sources: dict[str, str] | None = None,
) -> list[BidirectionalDomainMotifInteraction]:
    """Predict reverse domain-motif interactions from in-memory inputs."""

    filtered_regex, filtered_domains, filtered_sources = _filter_reverse_motif_resources(
        motif_regex = elm_regex,
        motif_domains = motif_domains,
        motif_sources = motif_sources,
    )

    return _predict_directional_domain_motif_interactions(
        motif_sequences = bacterial_sequences,
        motif_regex = filtered_regex,
        motif_domains = filtered_domains,
        partner_domains = human_domains,
        motif_sources = filtered_sources,
        motif_protein_side = 'microbe',
    )


def predict_bidirectional_domain_motif_interactions_from_data(
    elm_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    motif_sources: dict[str, str] | None = None,
    human_sequences: dict[str, str] | None = None,
    bacterial_domains: dict[str, list[str]] | None = None,
    bacterial_sequences: dict[str, str] | None = None,
    human_domains: dict[str, list[str]] | None = None,
    mode: str = 'both',
) -> list[BidirectionalDomainMotifInteraction]:
    """Predict forward, reverse, or bidirectional domain-motif interactions."""

    if mode not in {'forward', 'reverse', 'both'}:
        raise InputFormatError(
            '`mode` must be one of: forward, reverse, both.',
        )

    interactions: list[BidirectionalDomainMotifInteraction] = []

    if mode in {'forward', 'both'}:
        if human_sequences is None or bacterial_domains is None:
            raise InputFormatError(
                'Forward DMI mode requires `human_sequences` and '
                '`bacterial_domains`.',
            )

        interactions.extend(
            _predict_directional_domain_motif_interactions(
                motif_sequences = human_sequences,
                motif_regex = elm_regex,
                motif_domains = motif_domains,
                partner_domains = bacterial_domains,
                motif_sources = motif_sources,
                motif_protein_side = 'host',
            ),
        )

    if mode in {'reverse', 'both'}:
        if bacterial_sequences is None or human_domains is None:
            raise InputFormatError(
                'Reverse DMI mode requires `bacterial_sequences` and '
                '`human_domains`.',
            )

        filtered_regex, filtered_domains, filtered_sources = _filter_reverse_motif_resources(
            motif_regex = elm_regex,
            motif_domains = motif_domains,
            motif_sources = motif_sources,
        )

        interactions.extend(
            _predict_directional_domain_motif_interactions(
                motif_sequences = bacterial_sequences,
                motif_regex = filtered_regex,
                motif_domains = filtered_domains,
                partner_domains = human_domains,
                motif_sources = filtered_sources,
                motif_protein_side = 'microbe',
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
    """Predict forward domain-motif interactions from input files."""

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
        motif_sources = resolved_resources.motif_sources,
    )


def predict_reverse_domain_motif_interactions(
    bacterial_fasta_file: PathLike,
    human_domain_file: PathLike,
    elm_regex_file: PathLike | None = None,
    motif_domain_file: PathLike | None = None,
    resource_bundle: DMIResourceBundle | None = None,
) -> list[BidirectionalDomainMotifInteraction]:
    """Predict reverse domain-motif interactions from input files."""

    bacterial_sequences = read_fasta_sequences(bacterial_fasta_file)
    human_domains = read_protein_domain_table(human_domain_file)
    resolved_resources = _resolve_dmi_resource_bundle(
        elm_regex_file = elm_regex_file,
        motif_domain_file = motif_domain_file,
        resource_bundle = resource_bundle,
    )

    return predict_reverse_domain_motif_interactions_from_data(
        bacterial_sequences = bacterial_sequences,
        elm_regex = resolved_resources.elm_regex,
        motif_domains = resolved_resources.motif_domains,
        human_domains = human_domains,
        motif_sources = resolved_resources.motif_sources,
    )


def predict_bidirectional_domain_motif_interactions(
    human_fasta_file: PathLike | None = None,
    bacterial_domain_file: PathLike | None = None,
    bacterial_fasta_file: PathLike | None = None,
    human_domain_file: PathLike | None = None,
    elm_regex_file: PathLike | None = None,
    motif_domain_file: PathLike | None = None,
    resource_bundle: DMIResourceBundle | None = None,
    mode: str = 'both',
) -> list[BidirectionalDomainMotifInteraction]:
    """Predict forward, reverse, or bidirectional domain-motif interactions."""

    resolved_resources = _resolve_dmi_resource_bundle(
        elm_regex_file = elm_regex_file,
        motif_domain_file = motif_domain_file,
        resource_bundle = resource_bundle,
    )

    human_sequences = (
        read_fasta_sequences(human_fasta_file)
        if human_fasta_file is not None
        else None
    )
    bacterial_domains = (
        read_protein_domain_table(bacterial_domain_file)
        if bacterial_domain_file is not None
        else None
    )
    bacterial_sequences = (
        read_fasta_sequences(bacterial_fasta_file)
        if bacterial_fasta_file is not None
        else None
    )
    human_domains = (
        read_protein_domain_table(human_domain_file)
        if human_domain_file is not None
        else None
    )

    return predict_bidirectional_domain_motif_interactions_from_data(
        elm_regex = resolved_resources.elm_regex,
        motif_domains = resolved_resources.motif_domains,
        motif_sources = resolved_resources.motif_sources,
        human_sequences = human_sequences,
        bacterial_domains = bacterial_domains,
        bacterial_sequences = bacterial_sequences,
        human_domains = human_domains,
        mode = mode,
    )


def interactions_to_dataframe(
    interactions: list[DomainMotifInteraction],
) -> pd.DataFrame:
    """Convert forward interaction records to a data frame."""

    columns = [
        'human_protein',
        'motif',
        'start',
        'end',
        'bacterial_domain',
        'bacterial_protein',
        'resource',
    ]

    return pd.DataFrame(
        [asdict(interaction) for interaction in interactions],
        columns = columns,
    )


def bidirectional_interactions_to_dataframe(
    interactions: list[BidirectionalDomainMotifInteraction],
) -> pd.DataFrame:
    """Convert bidirectional interaction records to a data frame."""

    columns = [
        'host_protein',
        'microbial_protein',
        'motif',
        'start',
        'end',
        'domain',
        'motif_protein_side',
        'domain_protein_side',
        'resource',
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
    """Write forward predicted interactions to disk."""

    interaction_frame = interactions_to_dataframe(interactions)
    interaction_frame.to_csv(
        output_file,
        sep = separator,
        index = False,
    )


def write_bidirectional_domain_motif_interactions(
    interactions: list[BidirectionalDomainMotifInteraction],
    output_file: PathLike,
    separator: str = ';',
) -> None:
    """Write bidirectional predicted interactions to disk."""

    interaction_frame = bidirectional_interactions_to_dataframe(interactions)
    interaction_frame.to_csv(
        output_file,
        sep = separator,
        index = False,
    )
