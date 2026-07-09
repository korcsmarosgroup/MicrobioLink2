#!/usr/bin/env python

"""Domain-domain interaction utilities for the library-style MicrobioLink API."""

from __future__ import annotations

from dataclasses import asdict
from dataclasses import dataclass
from dataclasses import field
import importlib.resources
from pathlib import Path
from typing import Union

import pandas as pd

from microbiolink_api.dmi import read_protein_domain_table
import microbiolink_api.resources


PathLike = Union[str, Path]


@dataclass(frozen = True)
class DomainDomainInteraction:
    """Represent one host-microbe domain-domain interaction."""

    bacterial_protein: str
    human_protein: str
    bacterial_domain: str
    human_domain: str
    resource: str = 'custom'


@dataclass(frozen = True)
class DDIResourceBundle:
    """Collect the built-in or user-provided DDI resource tables."""

    pfam_pairs: set[tuple[str, str]]
    pair_sources: dict[tuple[str, str], tuple[str, ...]] = field(
        default_factory = dict,
    )
    name: str = 'custom'


def _canonicalize_pfam_pair(
    domain_a: str,
    domain_b: str,
) -> tuple[str, str]:
    """Return a deterministic representation of an undirected Pfam pair."""

    return tuple(sorted((domain_a, domain_b)))


def read_ddi_resource_table(filename: PathLike) -> set[tuple[str, str]]:
    """Read a Pfam-Pfam interaction table from disk."""

    pfam_pairs: set[tuple[str, str]] = set()

    with open(filename, encoding = 'utf-8') as ddi_table:
        for line in ddi_table:
            fields = line.strip().split('\t')

            if len(fields) < 2:
                continue

            pfam_pairs.add(_canonicalize_pfam_pair(fields[0], fields[1]))

    return pfam_pairs


def _build_pair_sources(
    pfam_pairs: set[tuple[str, str]],
    source_name: str,
) -> dict[tuple[str, str], tuple[str, ...]]:
    """Assign one source label to every Pfam pair in a resource bundle."""

    return {
        pfam_pair: (source_name,)
        for pfam_pair in pfam_pairs
    }


def merge_ddi_resource_bundles(
    *resource_bundles: DDIResourceBundle,
) -> DDIResourceBundle:
    """Merge multiple DDI resource bundles."""

    merged_pairs: set[tuple[str, str]] = set()
    merged_sources: dict[tuple[str, str], tuple[str, ...]] = {}
    merged_names: list[str] = []

    for bundle in resource_bundles:
        merged_pairs.update(bundle.pfam_pairs)

        if bundle.name not in merged_names:
            merged_names.append(bundle.name)

        for pfam_pair in bundle.pfam_pairs:
            existing_sources = list(merged_sources.get(pfam_pair, ()))
            bundle_sources = bundle.pair_sources.get(pfam_pair, (bundle.name,))

            for source_name in bundle_sources:
                if source_name not in existing_sources:
                    existing_sources.append(source_name)

            merged_sources[pfam_pair] = tuple(existing_sources)

    return DDIResourceBundle(
        pfam_pairs = merged_pairs,
        pair_sources = merged_sources,
        name = '+'.join(merged_names) if merged_names else 'custom',
    )


def load_default_3did_ddi_resource_bundle() -> DDIResourceBundle:
    """Load the packaged 3did Pfam-Pfam interaction table."""

    resource_path = importlib.resources.files(
        microbiolink_api.resources,
    ).joinpath('pfam_interactions_3did_current.tsv')
    pfam_pairs = read_ddi_resource_table(str(resource_path))

    return DDIResourceBundle(
        pfam_pairs = pfam_pairs,
        pair_sources = _build_pair_sources(pfam_pairs, '3did'),
        name = '3did',
    )


def load_default_domine_hc_ddi_resource_bundle() -> DDIResourceBundle:
    """Load the packaged high-confidence DOMINE Pfam-Pfam interaction table."""

    resource_path = importlib.resources.files(
        microbiolink_api.resources,
    ).joinpath('domine_v2_hc_pfam_pairs.tsv')
    pfam_pairs = read_ddi_resource_table(str(resource_path))

    return DDIResourceBundle(
        pfam_pairs = pfam_pairs,
        pair_sources = _build_pair_sources(pfam_pairs, 'DOMINE_hc'),
        name = 'DOMINE_hc',
    )


def load_default_domine_all_ddi_resource_bundle() -> DDIResourceBundle:
    """Load the packaged full DOMINE Pfam-Pfam interaction table."""

    resource_path = importlib.resources.files(
        microbiolink_api.resources,
    ).joinpath('domine_v2_all_pfam_pairs.tsv')
    pfam_pairs = read_ddi_resource_table(str(resource_path))

    return DDIResourceBundle(
        pfam_pairs = pfam_pairs,
        pair_sources = _build_pair_sources(pfam_pairs, 'DOMINE_all'),
        name = 'DOMINE_all',
    )


def load_default_ddi_resource_bundle() -> DDIResourceBundle:
    """Load the default DDI resource bundle for MicrobioLink 2.1."""

    return merge_ddi_resource_bundles(
        load_default_3did_ddi_resource_bundle(),
        load_default_domine_all_ddi_resource_bundle(),
    )


def _resolve_ddi_resource_bundle(
    ddi_resource_file: PathLike | None,
    resource_bundle: DDIResourceBundle | None,
) -> DDIResourceBundle:
    """Resolve built-in and user-provided DDI resource tables."""

    if ddi_resource_file is not None:
        pfam_pairs = read_ddi_resource_table(ddi_resource_file)
        return DDIResourceBundle(
            pfam_pairs = pfam_pairs,
            pair_sources = _build_pair_sources(pfam_pairs, 'custom'),
            name = 'custom',
        )

    if resource_bundle is not None:
        return resource_bundle

    return load_default_ddi_resource_bundle()


def predict_domain_domain_interactions_from_data(
    bacterial_domains: dict[str, list[str]],
    human_domains: dict[str, list[str]],
    resource_bundle: DDIResourceBundle | None = None,
) -> list[DomainDomainInteraction]:
    """Predict domain-domain interactions from in-memory inputs."""

    resolved_bundle = resource_bundle or load_default_ddi_resource_bundle()
    interactions: list[DomainDomainInteraction] = []
    seen_interactions: set[tuple[str, str, str, str, str]] = set()

    for bacterial_domain, bacterial_proteins in bacterial_domains.items():
        for human_domain, human_proteins in human_domains.items():
            pfam_pair = _canonicalize_pfam_pair(
                bacterial_domain,
                human_domain,
            )

            if pfam_pair not in resolved_bundle.pfam_pairs:
                continue

            resource_name = '|'.join(
                resolved_bundle.pair_sources.get(
                    pfam_pair,
                    (resolved_bundle.name,),
                ),
            )

            for bacterial_protein in bacterial_proteins:
                for human_protein in human_proteins:
                    interaction_key = (
                        bacterial_protein,
                        human_protein,
                        bacterial_domain,
                        human_domain,
                        resource_name,
                    )

                    if interaction_key in seen_interactions:
                        continue

                    seen_interactions.add(interaction_key)
                    interactions.append(
                        DomainDomainInteraction(
                            bacterial_protein = bacterial_protein,
                            human_protein = human_protein,
                            bacterial_domain = bacterial_domain,
                            human_domain = human_domain,
                            resource = resource_name,
                        ),
                    )

    return interactions


def predict_domain_domain_interactions(
    bacterial_domain_file: PathLike,
    human_domain_file: PathLike,
    ddi_resource_file: PathLike | None = None,
    resource_bundle: DDIResourceBundle | None = None,
) -> list[DomainDomainInteraction]:
    """Predict domain-domain interactions from input files."""

    bacterial_domains = read_protein_domain_table(bacterial_domain_file)
    human_domains = read_protein_domain_table(human_domain_file)
    resolved_bundle = _resolve_ddi_resource_bundle(
        ddi_resource_file = ddi_resource_file,
        resource_bundle = resource_bundle,
    )

    return predict_domain_domain_interactions_from_data(
        bacterial_domains = bacterial_domains,
        human_domains = human_domains,
        resource_bundle = resolved_bundle,
    )


def ddi_interactions_to_dataframe(
    interactions: list[DomainDomainInteraction],
) -> pd.DataFrame:
    """Convert DDI interaction records to a data frame."""

    columns = [
        'bacterial_protein',
        'human_protein',
        'bacterial_domain',
        'human_domain',
        'resource',
    ]

    return pd.DataFrame(
        [asdict(interaction) for interaction in interactions],
        columns = columns,
    )


def write_domain_domain_interactions(
    interactions: list[DomainDomainInteraction],
    output_file: PathLike,
    separator: str = ';',
) -> None:
    """Write predicted DDI interactions to disk."""

    interaction_frame = ddi_interactions_to_dataframe(interactions)
    interaction_frame.to_csv(
        output_file,
        sep = separator,
        index = False,
    )
