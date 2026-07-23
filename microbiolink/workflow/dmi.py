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


class MotifResource(NamedTuple):
    """One motif class's regex, compatible Pfam domains, and source."""

    regex: str
    domains: list[str]
    source: str


class RoleBasedMatch(NamedTuple):
    """One predicted interaction before translation into the bacterial/human output schema."""

    dmi_type: str
    motif_protein: str
    motif_id: str
    partner_protein: str
    domain_name: str
    start: int
    end: int
    resource: str


# SECTION 1: Resource loading and motif index building


def _read_tsv_table(filename: str) -> pd.DataFrame:
    """Read a packaged TSV resource, skipping any leading '#'-prefixed comment lines.

    Args:
        filename: Name of the packaged TSV resource file, relative to the data package.

    Returns:
        The parsed table, with any leading comment lines removed before parsing.
    """

    resource_path = importlib.resources.files(data).joinpath(filename)
    lines = resource_path.read_text(encoding="utf-8").splitlines()
    while lines and lines[0].startswith("#"):
        lines.pop(0)

    return pd.read_csv(io.StringIO("\n".join(lines)), sep="\t", quotechar='"')


def _read_motif_regex_table(filename: str) -> dict[str, str]:
    """Read a motif-class TSV (ELM or 3did) into a motif_id -> regex mapping.

    Args:
        filename: Name of the packaged motif-class TSV resource file.

    Returns:
        A dict mapping each motif_id to its regex pattern.
    """

    table = _read_tsv_table(filename)
    return dict(zip(table.iloc[:, 1], table.iloc[:, 4]))


def _read_motif_domain_table(filename: str) -> dict[str, list[str]]:
    """Read a motif-domain TSV (ELM or 3did) into a motif_id -> Pfam ID list mapping.

    Args:
        filename: Name of the packaged motif-domain TSV resource file.

    Returns:
        A dict mapping each motif_id to the list of Pfam domain IDs it is compatible with.
    """

    table = _read_tsv_table(filename)
    motif_domains: dict[str, list[str]] = {}

    for motif_id, domain in zip(table.iloc[:, 0], table.iloc[:, 1]):
        motif_domains.setdefault(motif_id, []).append(domain)

    return motif_domains


@functools.lru_cache(maxsize=None)
def _load_dmi_resources() -> dict[str, MotifResource]:
    """Load and merge the packaged ELM and 3did motif resources into one motif_id -> MotifResource index.

    Returns:
        A dict mapping each motif_id to its MotifResource (regex, compatible Pfam domains, and source).
    """

    motif_resources: dict[str, MotifResource] = {}

    for source_name, regex_filename in REGEX_RESOURCE_FILES.items():
        source_regex = _read_motif_regex_table(regex_filename)
        source_domains = _read_motif_domain_table(DOMAIN_RESOURCE_FILES[source_name])

        motif_resources.update(
            {
                motif_id: MotifResource(
                    regex=source_regex[motif_id], domains=domains, source=source_name
                )
                for motif_id, domains in source_domains.items()
                if motif_id in source_regex
            },
        )

    return motif_resources


# SECTION 2: Motif matching and role-based interaction prediction


def _build_motif_index(
    motif_resources: dict[str, MotifResource],
    partner_domains: dict[str, list[str]],
) -> dict[str, MotifResource]:
    """Restrict a motif index to classes with at least one domain present in partner_domains.

    Args:
        motif_resources: The full motif_id -> MotifResource index from _load_dmi_resources.
        partner_domains: Pfam ID -> partner protein accessions for this direction's partner species.

    Returns:
        A motif_id -> MotifResource index restricted to motifs with at least one domain present
        in partner_domains, with each MotifResource's domains list narrowed to just those domains.
    """

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
    """Find every (motif_id, start, end) regex match per UniProt accession.

    Args:
        sequences: FASTA header -> sequence mapping to search for motif matches.
        motif_index: The motif_id -> MotifResource index to match against each sequence.

    Returns:
        A dict mapping each UniProt accession to its list of (motif_id, start, end) regex
        match tuples. Accessions with no matches are omitted.
    """

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
    dmi_type: str,
) -> list[RoleBasedMatch]:
    """Predict one direction's role-based interaction matches.
    Intially finds all the regex motif matches in the motif-bearing species,
    then finds all compatible partner proteins based on the partmer species' Pfam domains,
    and returns a RoleBasedMatch for each combination.

    Args:
        motif_sequences: FASTA header -> sequence mapping for the motif-bearing species.
        motif_index: The motif_id -> MotifResource index restricted to this direction's partner domains.
        partner_domains: Pfam ID -> partner protein accessions for this direction's partner species.
        dmi_type: 'forward' if the motif-bearing species is human, 'reverse' if bacterial.

    Returns:
        A list of RoleBasedMatch instances, one per (motif hit, compatible domain, partner protein)
        combination.
    """

    motif_matches = _find_motif_matches(motif_sequences, motif_index)
    rows: list[RoleBasedMatch] = []

    for motif_protein, hits in motif_matches.items():
        for motif_id, start, end in hits:
            resource = motif_index[motif_id]

            for domain_name in resource.domains:
                for partner_protein in partner_domains[domain_name]:
                    rows.append(
                        RoleBasedMatch(
                            dmi_type=dmi_type,
                            motif_protein=motif_protein,
                            motif_id=motif_id,
                            partner_protein=partner_protein,
                            domain_name=domain_name,
                            start=start,
                            end=end,
                            resource=resource.source,
                        ),
                    )

    return rows


def _reverse_compatible_motif_resources(
    motif_resources: dict[str, MotifResource],
) -> dict[str, MotifResource]:
    """Drop CLV_-prefixed (cleavage-site) motif classes, incompatible with reverse-mode DMI.

    Args:
        motif_resources: The full motif_id -> MotifResource index.

    Returns:
        A motif_id -> MotifResource index with CLV_-prefixed motif classes removed.
    """

    return {
        motif_id: resource
        for motif_id, resource in motif_resources.items()
        if not motif_id.startswith(REVERSE_EXCLUDED_MOTIF_PREFIX)
    }


def _predict_direction(
    sequences: dict[str, str] | None,
    partner_domains: dict[str, list[str]] | None,
    motif_resources: dict[str, MotifResource],
    dmi_type: str,
) -> list[RoleBasedMatch]:
    """Validate inputs, then predict one direction's role-based interaction matches.

    Args:
        sequences: FASTA header -> sequence mapping for the motif-bearing species, or None
            if not supplied.
        partner_domains: Pfam ID -> partner protein accessions for this direction's partner
            species, or None if not supplied.
        motif_resources: The motif_id -> MotifResource index to search with (already filtered
            for reverse mode, if applicable).
        dmi_type: 'forward' if the motif-bearing species is human, 'reverse' if bacterial.

    Returns:
        A list of RoleBasedMatch instances predicted for this direction.

    Raises:
        ValueError: If sequences or partner_domains is None, meaning the required input for
            this direction's mode was not supplied.
    """

    if sequences is None or partner_domains is None:
        raise ValueError(f"{dmi_type} mode requires sequences and partner domains.")

    motif_index = _build_motif_index(motif_resources, partner_domains)
    return _predict_directional_interactions(
        sequences,
        motif_index,
        partner_domains,
        dmi_type=dmi_type,
    )


def _predict_rows(
    mode: str,
    human_sequences: dict[str, str] | None,
    bacterial_domains: dict[str, list[str]] | None,
    bacterial_sequences: dict[str, str] | None,
    human_domains: dict[str, list[str]] | None,
    motif_resources: dict[str, MotifResource],
) -> list[RoleBasedMatch]:
    """Predict combined role-based matches for whichever direction(s) mode selects.
    A match wiull be returned if the motif is found in the sequences and the domain found
    in the Pfam Ids for a domain-motif pair. Forward mode is a bacterial domain interacting
    with a human motif, and reverse mode is a human domain interacting with a bacterial motif.

    Args:
        mode: 'forward', 'reverse', or 'both'.
        human_sequences: Human FASTA header -> sequence mapping. Required for forward/both.
        bacterial_domains: Pfam ID -> bacterial UniProt accessions. Required for forward/both.
        bacterial_sequences: Bacterial FASTA header -> sequence mapping. Required for reverse/both.
        human_domains: Pfam ID -> human UniProt accessions. Required for reverse/both.
        motif_resources: The full motif_id -> MotifResource index.

    Returns:
        A list of RoleBasedMatch instances predicted across the selected direction(s).
    """

    rows: list[RoleBasedMatch] = []

    if mode in {"forward", "both"}:
        rows.extend(
            _predict_direction(
                human_sequences,
                bacterial_domains,
                motif_resources,
                dmi_type="forward",
            )
        )

    if mode in {"reverse", "both"}:
        reverse_resources = _reverse_compatible_motif_resources(motif_resources)
        rows.extend(
            _predict_direction(
                bacterial_sequences,
                human_domains,
                reverse_resources,
                dmi_type="reverse",
            )
        )

    return rows


# SECTION 3: Output translation and duplicate merging


def _merge_duplicate_matches(rows: list[RoleBasedMatch]) -> list[RoleBasedMatch]:
    """Merge role-based matches sharing protein pair, domain, and motif position across resources.

    Args:
        rows: The role-based matches to merge.

    Returns:
        A list of RoleBasedMatch instances, one per unique (dmi_type, motif_protein,
        partner_protein, domain_name, start, end) combination, with motif_id and resource
        joined by '|' across all contributing matches.
    """

    groups: dict[tuple, dict[str, list[str]]] = {}

    for match in rows:
        key = (
            match.dmi_type,
            match.motif_protein,
            match.partner_protein,
            match.domain_name,
            match.start,
            match.end,
        )
        if key not in groups:
            groups[key] = {"motif_ids": [], "resources": []}

        group = groups[key]
        if match.motif_id not in group["motif_ids"]:
            group["motif_ids"].append(match.motif_id)
        if match.resource not in group["resources"]:
            group["resources"].append(match.resource)

    merged_rows: list[RoleBasedMatch] = []
    for key, group in groups.items():
        dmi_type, motif_protein, partner_protein, domain_name, start, end = key
        merged_rows.append(
            RoleBasedMatch(
                dmi_type=dmi_type,
                motif_protein=motif_protein,
                motif_id="|".join(group["motif_ids"]),
                partner_protein=partner_protein,
                domain_name=domain_name,
                start=start,
                end=end,
                resource="|".join(group["resources"]),
            ),
        )

    return merged_rows


def _to_output_row(row: RoleBasedMatch) -> tuple:
    """Translate one role-based match into the bacterial/human output schema (forward: motif=human; reverse: motif=bacterial).

    Args:
        row: The role-based match to translate.

    Returns:
        A tuple matching OUTPUT_COLUMNS: (dmi_type, bacterial_uniprot_id, bacterial_annotation,
        human_uniprot_id, human_annotation, start, end, resource).
    """

    if row.dmi_type == "forward":
        bacterial_id, bacterial_annotation = row.partner_protein, row.domain_name
        human_id, human_annotation = row.motif_protein, row.motif_id
    else:
        bacterial_id, bacterial_annotation = row.motif_protein, row.motif_id
        human_id, human_annotation = row.partner_protein, row.domain_name

    return (
        row.dmi_type,
        bacterial_id,
        bacterial_annotation,
        human_id,
        human_annotation,
        row.start,
        row.end,
        row.resource,
    )


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
