"""Predicting domain-domain interactions between bacterial and human proteins."""

import functools
import importlib.resources

import pandas as pd

from .. import data

# Public API — Essential functions (see docs/api/ddi.md).
__all__ = [
    'predict_domain_domain_interactions',
]

OUTPUT_COLUMNS = [
    'bacterial_uniprot_id',
    'bacterial_pfam_domain',
    'human_uniprot_id',
    'human_pfam_domain',
    'resource',
]

RESOURCE_FILES = {
    '3did': 'pfam_interactions_3did_current.tsv',
    'DOMINE_hc': 'domine_v2_hc_pfam_pairs.tsv',
}


def _canonical_pair(domain_a: str, domain_b: str) -> tuple[str, str]:
    """Return a deterministic representation of an undirected Pfam pair."""

    first, second = sorted((domain_a, domain_b))
    return first, second


def _read_pfam_pair_table(filename: str) -> set[tuple[str, str]]:
    """Read a headerless two-column Pfam-Pfam interaction TSV into a canonical pair set."""

    resource_path = importlib.resources.files(data).joinpath(filename)
    pairs: set[tuple[str, str]] = set()

    with resource_path.open(encoding='utf-8') as pair_table:
        for line in pair_table:
            fields = line.strip().split('\t')
            if len(fields) < 2:
                continue
            pairs.add(_canonical_pair(fields[0], fields[1]))

    return pairs


@functools.lru_cache(maxsize=None)
def _load_ddi_resource_pairs() -> dict[tuple[str, str], list[str]]:
    """Load the packaged 3did and DOMINE high-confidence Pfam pair resources."""

    resource_pairs: dict[tuple[str, str], list[str]] = {}
    for source_name, filename in RESOURCE_FILES.items():
        for pair in _read_pfam_pair_table(filename):
            resource_pairs.setdefault(pair, []).append(source_name)

    return resource_pairs


def predict_domain_domain_interactions(
    bacterial_domains: dict[str, list[str]],
    human_domains: dict[str, list[str]],
) -> pd.DataFrame:
    """Predict domain-domain interactions between bacterial and human proteins.

    Args:
        bacterial_domains: Mapping of Pfam ID to bacterial UniProt accessions
            carrying that domain (Module 4's per-species output shape).
        human_domains: Mapping of Pfam ID to human UniProt accessions
            carrying that domain.

    Returns:
        A data frame with columns bacterial_uniprot_id, bacterial_pfam_domain,
        human_uniprot_id, human_pfam_domain, and resource ('3did', 'DOMINE_hc',
        or '3did|DOMINE_hc'), one row per known-interacting domain pair
        instance between a bacterial and a human protein.
    """

    resource_pairs = _load_ddi_resource_pairs()
    rows = []
    seen = set()

    for bacterial_pfam, bacterial_proteins in bacterial_domains.items():
        for human_pfam, human_proteins in human_domains.items():
            pair = _canonical_pair(bacterial_pfam, human_pfam)
            sources = resource_pairs.get(pair)
            if sources is None:
                continue

            resource = '|'.join(sources)
            for bacterial_protein in bacterial_proteins:
                for human_protein in human_proteins:
                    key = (bacterial_protein, human_protein, bacterial_pfam, human_pfam)
                    if key in seen:
                        continue
                    seen.add(key)
                    rows.append((bacterial_protein, bacterial_pfam, human_protein, human_pfam, resource))

    return pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
