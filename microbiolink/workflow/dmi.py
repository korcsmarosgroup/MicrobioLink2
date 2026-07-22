"""Predicting domain-motif interactions between bacterial and human proteins."""

import functools
import importlib.resources
import re
from typing import NamedTuple

import pandas as pd

from .. import data
from ..utils import fasta

OUTPUT_COLUMNS = [
    'dmi_type',
    'bacterial_uniprot_id',
    'bacterial_annotation',
    'human_uniprot_id',
    'human_annotation',
    'start',
    'end',
    'resource',
]

REGEX_RESOURCE_FILES = {
    'ELM': 'elm_classes.tsv',
    '3did': '3did_dmi_classes.tsv',
}
DOMAIN_RESOURCE_FILES = {
    'ELM': 'elm_interaction_domains.tsv',
    '3did': '3did_dmi_interaction_domains.tsv',
}
REVERSE_EXCLUDED_MOTIF_PREFIXES = ('CLV_',)


def _skip_comments_and_header(lines):
    """Advance a line iterator past leading '#' comment lines and one header row."""

    for line in lines:
        if line.startswith('#'):
            continue
        break


def _read_motif_regex_table(filename: str) -> dict[str, str]:
    """Read a motif-class TSV (ELM or 3did) into a motif_id -> regex mapping."""

    resource_path = importlib.resources.files(data).joinpath(filename)
    motif_regex: dict[str, str] = {}

    with resource_path.open(encoding='utf-8') as regex_table:
        lines = iter(regex_table)
        _skip_comments_and_header(lines)

        for line in lines:
            fields = line.replace('"', '').strip().split('\t')
            if len(fields) > 4:
                motif_regex[fields[1]] = fields[4]

    return motif_regex


def _read_motif_domain_table(filename: str) -> dict[str, list[str]]:
    """Read a motif-domain TSV (ELM or 3did) into a motif_id -> Pfam ID list mapping."""

    resource_path = importlib.resources.files(data).joinpath(filename)
    motif_domains: dict[str, list[str]] = {}

    with resource_path.open(encoding='utf-8') as domain_table:
        lines = iter(domain_table)
        _skip_comments_and_header(lines)

        for line in lines:
            fields = line.replace('"', '').strip().split('\t')
            if len(fields) > 1:
                motif_domains.setdefault(fields[0], []).append(fields[1])

    return motif_domains


@functools.lru_cache(maxsize=None)
def _load_dmi_resources() -> tuple[dict[str, str], dict[str, list[str]], dict[str, str]]:
    """Load and merge the packaged ELM and 3did motif-regex/motif-domain resources.

    Returns:
        (motif_regex, motif_domains, motif_sources) — motif_sources maps
        motif_id to 'ELM' or '3did'.
    """

    motif_regex: dict[str, str] = {}
    motif_domains: dict[str, list[str]] = {}
    motif_sources: dict[str, str] = {}

    for source_name, regex_filename in REGEX_RESOURCE_FILES.items():
        source_regex = _read_motif_regex_table(regex_filename)
        source_domains = _read_motif_domain_table(DOMAIN_RESOURCE_FILES[source_name])

        motif_regex.update(source_regex)
        for motif_id, domains in source_domains.items():
            motif_domains.setdefault(motif_id, []).extend(domains)
        for motif_id in set(source_regex) | set(source_domains):
            motif_sources[motif_id] = source_name

    return motif_regex, motif_domains, motif_sources


def _filter_reverse_motifs(
    motif_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    motif_sources: dict[str, str],
) -> tuple[dict[str, str], dict[str, list[str]], dict[str, str]]:
    """Drop CLV_-prefixed (cleavage-site) motif classes for reverse-mode DMI."""

    def is_reverse_compatible(motif_id: str) -> bool:
        return not motif_id.startswith(REVERSE_EXCLUDED_MOTIF_PREFIXES)

    filtered_regex = {
        motif_id: pattern for motif_id, pattern in motif_regex.items() if is_reverse_compatible(motif_id)
    }
    filtered_domains = {
        motif_id: domains for motif_id, domains in motif_domains.items() if is_reverse_compatible(motif_id)
    }
    filtered_sources = {
        motif_id: source for motif_id, source in motif_sources.items() if is_reverse_compatible(motif_id)
    }

    return filtered_regex, filtered_domains, filtered_sources


class MotifResource(NamedTuple):
    """One motif class's regex, domain-relevance-filtered compatible domains, and source."""

    regex: str
    domains: list[str]
    source: str


def _build_motif_index(
    motif_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    motif_sources: dict[str, str],
    partner_domains: dict[str, list[str]],
) -> dict[str, MotifResource]:
    """Index motif classes with at least one compatible domain present in partner_domains."""

    motif_index: dict[str, MotifResource] = {}

    for motif_id, domains in motif_domains.items():
        pattern = motif_regex.get(motif_id)
        if pattern is None:
            continue

        relevant_domains = [domain for domain in domains if domain in partner_domains]
        if not relevant_domains:
            continue

        motif_index[motif_id] = MotifResource(
            regex=pattern,
            domains=relevant_domains,
            source=motif_sources.get(motif_id, 'ELM'),
        )

    return motif_index


def _find_motif_matches(
    sequences: dict[str, str],
    motif_index: dict[str, MotifResource],
) -> dict[str, list[tuple[str, int, int]]]:
    """Find every (motif_id, start, end) regex match per UniProt accession."""

    motif_matches: dict[str, list[tuple[str, int, int]]] = {}

    for header, sequence in sequences.items():
        uniprot_id = fasta.extract_uniprot_id(header)
        hits = [
            (motif_id, match.start(), match.end())
            for motif_id, resource in motif_index.items()
            for match in re.finditer(resource.regex, sequence)
        ]

        if hits:
            motif_matches[uniprot_id] = hits

    return motif_matches


def _predict_directional_interactions(
    motif_sequences: dict[str, str],
    motif_index: dict[str, MotifResource],
    partner_domains: dict[str, list[str]],
    motif_is_bacterial: bool,
) -> list[tuple]:
    """Predict one directional set of domain-motif interaction rows."""

    motif_matches = _find_motif_matches(motif_sequences, motif_index)
    dmi_type = 'reverse' if motif_is_bacterial else 'forward'
    rows: list[tuple] = []

    for motif_protein, hits in motif_matches.items():
        for motif_id, start, end in hits:
            resource = motif_index[motif_id]

            for domain_name in resource.domains:
                for partner_protein in partner_domains[domain_name]:
                    if motif_is_bacterial:
                        row = (
                            dmi_type,
                            motif_protein,
                            motif_id,
                            partner_protein,
                            domain_name,
                            start,
                            end,
                            resource.source,
                        )
                    else:
                        row = (
                            dmi_type,
                            partner_protein,
                            domain_name,
                            motif_protein,
                            motif_id,
                            start,
                            end,
                            resource.source,
                        )
                    rows.append(row)

    return rows


def _merge_duplicate_matches(rows: list[tuple]) -> list[tuple]:
    """Merge rows that share protein pair, domain, and motif position across resources."""

    groups: dict[tuple, dict[str, list[str]]] = {}
    order: list[tuple] = []

    for dmi_type, bacterial_id, bacterial_annotation, human_id, human_annotation, start, end, resource in rows:
        if dmi_type == 'forward':
            domain_annotation, motif_annotation = bacterial_annotation, human_annotation
        else:
            domain_annotation, motif_annotation = human_annotation, bacterial_annotation

        key = (dmi_type, bacterial_id, human_id, start, end, domain_annotation)
        if key not in groups:
            groups[key] = {'motif_annotations': [], 'resources': []}
            order.append(key)

        group = groups[key]
        if motif_annotation not in group['motif_annotations']:
            group['motif_annotations'].append(motif_annotation)
        if resource not in group['resources']:
            group['resources'].append(resource)

    merged_rows: list[tuple] = []
    for key in order:
        dmi_type, bacterial_id, human_id, start, end, domain_annotation = key
        group = groups[key]
        motif_annotation = '|'.join(group['motif_annotations'])
        resource = '|'.join(group['resources'])

        if dmi_type == 'forward':
            bacterial_annotation, human_annotation = domain_annotation, motif_annotation
        else:
            bacterial_annotation, human_annotation = motif_annotation, domain_annotation

        merged_rows.append(
            (dmi_type, bacterial_id, bacterial_annotation, human_id, human_annotation, start, end, resource),
        )

    return merged_rows


def predict_domain_motif_interactions(
    mode: str,
    human_sequences: dict[str, str] | None = None,
    bacterial_domains: dict[str, list[str]] | None = None,
    bacterial_sequences: dict[str, str] | None = None,
    human_domains: dict[str, list[str]] | None = None,
) -> pd.DataFrame:
    """Predict forward and/or reverse domain-motif interactions.

    Args:
        mode: 'forward', 'reverse', or 'both'.
        human_sequences: Human FASTA header -> sequence mapping (Module 3's
            output shape, read via utils.fasta.read_fasta_sequences).
            Required for forward/both.
        bacterial_domains: Mapping of Pfam ID to bacterial UniProt accessions
            carrying that domain (Module 4's output shape). Required for
            forward/both.
        bacterial_sequences: Bacterial FASTA header -> sequence mapping.
            Required for reverse/both.
        human_domains: Mapping of Pfam ID to human UniProt accessions
            carrying that domain. Required for reverse/both.

    Returns:
        A data frame with columns dmi_type, bacterial_uniprot_id,
        bacterial_annotation, human_uniprot_id, human_annotation, start, end,
        and resource, one row per predicted domain-motif interaction.
    """

    if mode not in {'forward', 'reverse', 'both'}:
        raise ValueError(f"mode must be 'forward', 'reverse', or 'both', got {mode!r}")

    motif_regex, motif_domains, motif_sources = _load_dmi_resources()
    rows: list[tuple] = []

    if mode in {'forward', 'both'}:
        if human_sequences is None or bacterial_domains is None:
            raise ValueError('forward mode requires human_sequences and bacterial_domains.')

        motif_index = _build_motif_index(motif_regex, motif_domains, motif_sources, bacterial_domains)
        rows.extend(
            _predict_directional_interactions(
                human_sequences,
                motif_index,
                bacterial_domains,
                motif_is_bacterial=False,
            ),
        )

    if mode in {'reverse', 'both'}:
        if bacterial_sequences is None or human_domains is None:
            raise ValueError('reverse mode requires bacterial_sequences and human_domains.')

        filtered_regex, filtered_domains, filtered_sources = _filter_reverse_motifs(
            motif_regex,
            motif_domains,
            motif_sources,
        )
        motif_index = _build_motif_index(filtered_regex, filtered_domains, filtered_sources, human_domains)
        rows.extend(
            _predict_directional_interactions(
                bacterial_sequences,
                motif_index,
                human_domains,
                motif_is_bacterial=True,
            ),
        )

    return pd.DataFrame(_merge_duplicate_matches(rows), columns=OUTPUT_COLUMNS)
