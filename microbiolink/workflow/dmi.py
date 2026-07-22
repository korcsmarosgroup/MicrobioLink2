"""Predicting domain-motif interactions between bacterial and human proteins."""

import functools
import importlib.resources
import io
import re
from typing import NamedTuple

import pandas as pd

from .. import data
from ..utils import fasta

OUTPUT_COLUMNS = [
    "dmi_type",
    "bacterial_uniprot_id",
    "bacterial_annotation",
    "human_uniprot_id",
    "human_annotation",
    "start",
    "end",
    "resource",
]

REGEX_RESOURCE_FILES = {
    "ELM": "elm_classes.tsv",
    "3did": "3did_dmi_classes.tsv",
}
DOMAIN_RESOURCE_FILES = {
    "ELM": "elm_interaction_domains.tsv",
    "3did": "3did_dmi_interaction_domains.tsv",
}
REVERSE_EXCLUDED_MOTIF_PREFIX = "CLV_"


def _read_tsv_table(filename: str) -> pd.DataFrame:
    """Read a packaged TSV resource, skipping any leading '#'-prefixed comment lines."""

    resource_path = importlib.resources.files(data).joinpath(filename)
    lines = resource_path.read_text(encoding="utf-8").splitlines()
    while lines and lines[0].startswith("#"):
        lines.pop(0)

    return pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", quotechar='"')


def _read_motif_regex_table(filename: str) -> dict[str, str]:
    """Read a motif-class TSV (ELM or 3did) into a motif_id -> regex mapping."""

    table = _read_tsv_table(filename)
    return dict(zip(table.iloc[:, 1], table.iloc[:, 4]))


def _read_motif_domain_table(filename: str) -> dict[str, list[str]]:
    """Read a motif-domain TSV (ELM or 3did) into a motif_id -> Pfam ID list mapping."""

    table = _read_tsv_table(filename)
    motif_domains: dict[str, list[str]] = {}

    for motif_id, domain in zip(table.iloc[:, 0], table.iloc[:, 1]):
        motif_domains.setdefault(motif_id, []).append(domain)

    return motif_domains


class MotifResource(NamedTuple):
    """One motif class's regex, compatible Pfam domains, and source."""

    regex: str
    domains: list[str]
    source: str


@functools.lru_cache(maxsize=None)
def _load_dmi_resources() -> dict[str, MotifResource]:
    """Load and merge the packaged ELM and 3did motif resources into one motif_id -> MotifResource index."""

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

    return {
        motif_id: MotifResource(
            regex=motif_regex[motif_id],
            domains=domains,
            source=motif_sources[motif_id],
        )
        for motif_id, domains in motif_domains.items()
        if motif_id in motif_regex
    }


def _build_motif_index(
    motif_resources: dict[str, MotifResource],
    partner_domains: dict[str, list[str]],
) -> dict[str, MotifResource]:
    """Restrict a motif index to classes with at least one domain present in partner_domains."""

    motif_index: dict[str, MotifResource] = {}

    for motif_id, resource in motif_resources.items():
        relevant_domains = [
            domain for domain in resource.domains if domain in partner_domains
        ]
        if relevant_domains:
            motif_index[motif_id] = resource._replace(domains=relevant_domains)

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
    """Predict one direction's role-based rows: (dmi_type, motif_protein, motif_id, partner_protein, domain_name, start, end, resource)."""

    motif_matches = _find_motif_matches(motif_sequences, motif_index)
    dmi_type = "reverse" if motif_is_bacterial else "forward"
    rows: list[tuple] = []

    for motif_protein, hits in motif_matches.items():
        for motif_id, start, end in hits:
            resource = motif_index[motif_id]

            for domain_name in resource.domains:
                for partner_protein in partner_domains[domain_name]:
                    rows.append(
                        (
                            dmi_type,
                            motif_protein,
                            motif_id,
                            partner_protein,
                            domain_name,
                            start,
                            end,
                            resource.source,
                        ),
                    )

    return rows


def _group_matches_by_key(rows: list[tuple]) -> dict[tuple, dict[str, list[str]]]:
    """Group role-based rows sharing protein pair, domain, and motif position, collecting motif_ids/resources per group."""

    groups: dict[tuple, dict[str, list[str]]] = {}

    for (
        dmi_type,
        motif_protein,
        motif_id,
        partner_protein,
        domain_name,
        start,
        end,
        resource,
    ) in rows:
        key = (dmi_type, motif_protein, partner_protein, domain_name, start, end)
        if key not in groups:
            groups[key] = {"motif_ids": [], "resources": []}

        group = groups[key]
        if motif_id not in group["motif_ids"]:
            group["motif_ids"].append(motif_id)
        if resource not in group["resources"]:
            group["resources"].append(resource)

    return groups


def _merge_duplicate_matches(rows: list[tuple]) -> list[tuple]:
    """Merge role-based rows sharing protein pair, domain, and motif position across resources."""

    groups = _group_matches_by_key(rows)

    merged_rows: list[tuple] = []
    for key, group in groups.items():
        dmi_type, motif_protein, partner_protein, domain_name, start, end = key
        motif_id = "|".join(group["motif_ids"])
        resource = "|".join(group["resources"])
        merged_rows.append(
            (
                dmi_type,
                motif_protein,
                motif_id,
                partner_protein,
                domain_name,
                start,
                end,
                resource,
            )
        )

    return merged_rows


def _to_output_row(row: tuple) -> tuple:
    """Translate one role-based row into the bacterial/human output schema (forward: motif=human; reverse: motif=bacterial)."""

    (
        dmi_type,
        motif_protein,
        motif_id,
        partner_protein,
        domain_name,
        start,
        end,
        resource,
    ) = row

    if dmi_type == "forward":
        bacterial_id, bacterial_annotation = partner_protein, domain_name
        human_id, human_annotation = motif_protein, motif_id
    else:
        bacterial_id, bacterial_annotation = motif_protein, motif_id
        human_id, human_annotation = partner_protein, domain_name

    return (
        dmi_type,
        bacterial_id,
        bacterial_annotation,
        human_id,
        human_annotation,
        start,
        end,
        resource,
    )


def _reverse_compatible_motif_resources(
    motif_resources: dict[str, MotifResource],
) -> dict[str, MotifResource]:
    """Drop CLV_-prefixed (cleavage-site) motif classes, incompatible with reverse-mode DMI."""

    return {
        motif_id: resource
        for motif_id, resource in motif_resources.items()
        if not motif_id.startswith(REVERSE_EXCLUDED_MOTIF_PREFIX)
    }


def _predict_direction(
    sequences: dict[str, str] | None,
    partner_domains: dict[str, list[str]] | None,
    motif_resources: dict[str, MotifResource],
    motif_is_bacterial: bool,
) -> list[tuple]:
    """Validate inputs, then predict one direction's role-based interaction rows."""

    if sequences is None or partner_domains is None:
        mode_name = "reverse" if motif_is_bacterial else "forward"
        raise ValueError(f"{mode_name} mode requires sequences and partner domains.")

    motif_index = _build_motif_index(motif_resources, partner_domains)
    return _predict_directional_interactions(
        sequences,
        motif_index,
        partner_domains,
        motif_is_bacterial=motif_is_bacterial,
    )


def _predict_rows(
    mode: str,
    human_sequences: dict[str, str] | None,
    bacterial_domains: dict[str, list[str]] | None,
    bacterial_sequences: dict[str, str] | None,
    human_domains: dict[str, list[str]] | None,
    motif_resources: dict[str, MotifResource],
) -> list[tuple]:
    """Predict combined role-based rows for whichever direction(s) mode selects."""

    rows: list[tuple] = []

    if mode in {"forward", "both"}:
        rows.extend(
            _predict_direction(
                human_sequences,
                bacterial_domains,
                motif_resources,
                motif_is_bacterial=False,
            )
        )

    if mode in {"reverse", "both"}:
        reverse_resources = _reverse_compatible_motif_resources(motif_resources)
        rows.extend(
            _predict_direction(
                bacterial_sequences,
                human_domains,
                reverse_resources,
                motif_is_bacterial=True,
            )
        )

    return rows


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
        human_sequences: Human FASTA header -> sequence mapping (Module 3's output). Required for forward/both.
        bacterial_domains: Pfam ID -> bacterial UniProt accessions (Module 4's output). Required for forward/both.
        bacterial_sequences: Bacterial FASTA header -> sequence mapping. Required for reverse/both.
        human_domains: Pfam ID -> human UniProt accessions. Required for reverse/both.

    Returns:
        A data frame with columns dmi_type, bacterial_uniprot_id, bacterial_annotation, human_uniprot_id, human_annotation, start, end, and resource, one row per predicted domain-motif interaction.

    Raises:
        ValueError: If mode is not 'forward', 'reverse', or 'both', or if a required input for the selected mode is missing.
    """

    if mode not in {"forward", "reverse", "both"}:
        raise ValueError(f"mode must be 'forward', 'reverse', or 'both', got {mode!r}")

    motif_resources = _load_dmi_resources()
    rows = _predict_rows(
        mode,
        human_sequences,
        bacterial_domains,
        bacterial_sequences,
        human_domains,
        motif_resources,
    )
    output_rows = [_to_output_row(row) for row in _merge_duplicate_matches(rows)]

    return pd.DataFrame(output_rows, columns=OUTPUT_COLUMNS)
