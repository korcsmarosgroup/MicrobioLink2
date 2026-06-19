#!/usr/bin/env python3

"""Benchmark broad PHISTO bacteria-human interactions with MicrobioLink rules.

This script benchmarks three upstream interaction approaches against a broad
true-positive interaction panel exported from PHISTO:

1. forward DMI using the packaged ELM + 3did DMI resources,
2. reverse DMI using the packaged ELM + 3did DMI resources,
3. DDI using 3did and DOMINE Pfam-Pfam resources.

The script is intentionally standard-library only so it can run in a minimal
environment while still following the same resource logic as the current
MicrobioLink codebase.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from dataclasses import field
from pathlib import Path
import re
import time
from typing import Iterable
from typing import Iterator
from urllib.error import HTTPError
from urllib.error import URLError
from urllib.parse import urlencode
from urllib.request import urlopen


UNIPROT_STREAM_URL = 'https://rest.uniprot.org/uniprotkb/stream'
UNIPROT_SEARCH_URL = 'https://rest.uniprot.org/uniprotkb/search'
UNIPROT_FIELDS = 'accession,id,xref_pfam,sequence'
UNIPROT_BATCH_SIZE = 100
REQUEST_RETRIES = 4
REQUEST_BACKOFF_SECONDS = 2


@dataclass(frozen = True)
class PhistoRow:
    """Represent one raw PHISTO export row."""

    pathogen_name: str
    taxonomy_id: str
    pathogen_accession: str
    pathogen_protein_name: str
    human_accession: str
    human_protein_name: str
    experimental_method: str
    pubmed_id: str


@dataclass
class BenchmarkPair:
    """Aggregate one unique pathogen-human pair from PHISTO."""

    pathogen_accession: str
    human_accession: str
    pathogen_name: str
    pathogen_protein_name: str
    human_protein_name: str
    taxonomy_ids: set[str] = field(default_factory = set)
    methods: set[str] = field(default_factory = set)
    pubmed_ids: set[str] = field(default_factory = set)
    supporting_rows: int = 0


@dataclass(frozen = True)
class UniProtAnnotation:
    """Store the sequence and Pfam annotations resolved for one accession."""

    requested_accession: str
    current_accession: str
    entry_name: str
    pfams: tuple[str, ...]
    sequence: str
    resolution_mode: str


@dataclass(frozen = True)
class DMIResources:
    """Store merged motif regexes and motif-domain relationships."""

    motif_patterns: dict[str, re.Pattern[str]]
    motif_domains: dict[str, tuple[str, ...]]
    motif_sources: dict[str, str]


@dataclass(frozen = True)
class MotifPresence:
    """Summarize motif presence on one protein."""

    count: int
    first_start: int
    first_end: int


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description = 'Benchmark PHISTO bacteria-human interactions with '
        'forward DMI, reverse DMI, and DDI.',
    )
    parser.add_argument(
        '--phisto_csv',
        required = True,
        help = 'Path to the PHISTO bacteria-only CSV export.',
    )
    parser.add_argument(
        '--domine_interaction_file',
        required = True,
        help = 'Path to DOMINE INTERACTION.txt from the v2.0 archive.',
    )
    parser.add_argument(
        '--output_dir',
        required = True,
        help = 'Directory where the benchmark outputs will be written.',
    )
    return parser.parse_args()


def ensure_directory(path: Path) -> None:
    """Create a directory if needed."""

    path.mkdir(parents = True, exist_ok = True)


def read_phisto_rows(path: Path) -> list[PhistoRow]:
    """Read a PHISTO bacteria export."""

    rows: list[PhistoRow] = []

    with open(path, newline = '', encoding = 'latin-1') as infile:
        reader = csv.reader(infile)
        next(reader, None)

        for fields in reader:
            if len(fields) < 8:
                continue

            rows.append(
                PhistoRow(
                    pathogen_name = fields[0].strip(),
                    taxonomy_id = fields[1].strip(),
                    pathogen_accession = fields[2].strip(),
                    pathogen_protein_name = fields[3].strip(),
                    human_accession = fields[4].strip(),
                    human_protein_name = fields[5].strip(),
                    experimental_method = fields[6].strip(),
                    pubmed_id = fields[7].strip(),
                ),
            )

    return rows


def aggregate_benchmark_pairs(rows: Iterable[PhistoRow]) -> dict[tuple[str, str], BenchmarkPair]:
    """Aggregate PHISTO rows into unique pathogen-human pairs."""

    pairs: dict[tuple[str, str], BenchmarkPair] = {}

    for row in rows:
        key = (row.pathogen_accession, row.human_accession)

        if key not in pairs:
            pairs[key] = BenchmarkPair(
                pathogen_accession = row.pathogen_accession,
                human_accession = row.human_accession,
                pathogen_name = row.pathogen_name,
                pathogen_protein_name = row.pathogen_protein_name,
                human_protein_name = row.human_protein_name,
            )

        pair = pairs[key]
        pair.taxonomy_ids.add(row.taxonomy_id)
        pair.methods.add(row.experimental_method)
        pair.pubmed_ids.add(row.pubmed_id)
        pair.supporting_rows += 1

    return pairs


def write_tsv(
    path: Path,
    rows: Iterable[dict[str, object]],
    fieldnames: list[str],
) -> None:
    """Write dictionaries to a TSV file."""

    with open(path, 'w', newline = '', encoding = 'utf-8') as outfile:
        writer = csv.DictWriter(
            outfile,
            fieldnames = fieldnames,
            delimiter = '\t',
            extrasaction = 'ignore',
        )
        writer.writeheader()

        for row in rows:
            writer.writerow(row)


def write_text_lines(path: Path, values: Iterable[str]) -> None:
    """Write one value per line."""

    with open(path, 'w', encoding = 'utf-8') as outfile:
        for value in values:
            outfile.write(f'{value}\n')


def build_uniprot_query(accessions: Iterable[str]) -> str:
    """Build a UniProt accession query."""

    clauses = [f'accession:{accession}' for accession in accessions]
    return '(' + ' OR '.join(clauses) + ')'


def fetch_text(url: str) -> str:
    """Fetch text from a URL with small retry logic."""

    delay = REQUEST_BACKOFF_SECONDS

    for attempt in range(REQUEST_RETRIES):
        try:
            with urlopen(url, timeout = 120) as response:
                return response.read().decode('utf-8')
        except HTTPError as error:
            retryable = error.code in {429, 500, 502, 503, 504}

            if not retryable or attempt == REQUEST_RETRIES - 1:
                raise
        except URLError:
            if attempt == REQUEST_RETRIES - 1:
                raise

        time.sleep(delay)
        delay *= 2

    raise RuntimeError(f'Failed to fetch URL after retries: {url}')


def parse_uniprot_tsv(
    text: str,
    requested_accession: str | None = None,
    resolution_mode: str = 'direct',
) -> list[UniProtAnnotation]:
    """Parse a UniProt TSV response."""

    if not text.strip():
        return []

    reader = csv.DictReader(text.splitlines(), delimiter = '\t')
    annotations: list[UniProtAnnotation] = []

    for row in reader:
        current_accession = row.get('Entry', '').strip()

        if not current_accession:
            continue

        pfam_field = row.get('Pfam', '').strip()
        pfams = tuple(
            pfam.strip()
            for pfam in pfam_field.split(';')
            if pfam.strip()
        )

        annotations.append(
            UniProtAnnotation(
                requested_accession = requested_accession or current_accession,
                current_accession = current_accession,
                entry_name = row.get('Entry Name', '').strip(),
                pfams = pfams,
                sequence = row.get('Sequence', '').strip(),
                resolution_mode = resolution_mode,
            ),
        )

    return annotations


def fetch_uniprot_batch(accessions: list[str]) -> list[UniProtAnnotation]:
    """Fetch one batch of UniProt annotations."""

    query = build_uniprot_query(accessions)
    params = {
        'fields': UNIPROT_FIELDS,
        'format': 'tsv',
        'query': query,
        'size': str(len(accessions)),
    }
    url = UNIPROT_SEARCH_URL + '?' + urlencode(params)
    return parse_uniprot_tsv(fetch_text(url))


def fetch_uniprot_single(accession: str) -> UniProtAnnotation | None:
    """Resolve one accession from UniProt."""

    params = {
        'fields': UNIPROT_FIELDS,
        'format': 'tsv',
        'query': f'accession:{accession}',
        'size': '1',
    }
    url = UNIPROT_SEARCH_URL + '?' + urlencode(params)
    annotations = parse_uniprot_tsv(
        fetch_text(url),
        requested_accession = accession,
        resolution_mode = 'search_fallback',
    )

    if not annotations:
        return None

    return annotations[0]


def resolve_uniprot_annotations(accessions: Iterable[str]) -> dict[str, UniProtAnnotation]:
    """Resolve UniProt sequence and Pfam annotations for requested accessions."""

    accession_list = sorted(set(accessions))
    resolved: dict[str, UniProtAnnotation] = {}

    for start in range(0, len(accession_list), UNIPROT_BATCH_SIZE):
        batch = accession_list[start : start + UNIPROT_BATCH_SIZE]
        batch_annotations = fetch_uniprot_batch(batch)

        for annotation in batch_annotations:
            if annotation.current_accession in batch:
                resolved[annotation.current_accession] = UniProtAnnotation(
                    requested_accession = annotation.current_accession,
                    current_accession = annotation.current_accession,
                    entry_name = annotation.entry_name,
                    pfams = annotation.pfams,
                    sequence = annotation.sequence,
                    resolution_mode = 'direct',
                )

    unresolved = [
        accession
        for accession in accession_list
        if accession not in resolved
    ]

    for accession in unresolved:
        annotation = fetch_uniprot_single(accession)

        if annotation is not None:
            resolved[accession] = annotation

    return resolved


def write_annotation_table(path: Path, annotations: dict[str, UniProtAnnotation]) -> None:
    """Write resolved UniProt annotations to TSV."""

    rows = []

    for accession in sorted(annotations):
        annotation = annotations[accession]
        rows.append(
            {
                'requested_accession': annotation.requested_accession,
                'current_accession': annotation.current_accession,
                'entry_name': annotation.entry_name,
                'pfam_count': len(annotation.pfams),
                'pfams': ';'.join(annotation.pfams),
                'sequence_length': len(annotation.sequence),
                'resolution_mode': annotation.resolution_mode,
            },
        )

    write_tsv(
        path,
        rows,
        [
            'requested_accession',
            'current_accession',
            'entry_name',
            'pfam_count',
            'pfams',
            'sequence_length',
            'resolution_mode',
        ],
    )


def write_fasta(path: Path, annotations: dict[str, UniProtAnnotation]) -> None:
    """Write a FASTA file keyed by the requested accession."""

    with open(path, 'w', encoding = 'utf-8') as outfile:
        for accession in sorted(annotations):
            annotation = annotations[accession]

            if not annotation.sequence:
                continue

            header = f'sp|{annotation.requested_accession}|{annotation.entry_name}'
            outfile.write(f'>{header}\n')
            outfile.write(f'{annotation.sequence}\n')


def write_domain_table(path: Path, annotations: dict[str, UniProtAnnotation]) -> None:
    """Write a simple protein-to-Pfam table."""

    with open(path, 'w', encoding = 'utf-8', newline = '') as outfile:
        writer = csv.writer(outfile, delimiter = '\t')
        writer.writerow(['protein', 'pfam'])

        for accession in sorted(annotations):
            annotation = annotations[accession]

            if not annotation.pfams:
                continue

            writer.writerow([annotation.requested_accession, ';'.join(annotation.pfams)])


def parse_elm_regex_table(path: Path) -> dict[str, str]:
    """Parse an ELM-style motif regex table."""

    motif_regex: dict[str, str] = {}

    with open(path, encoding = 'utf-8') as infile:
        next(infile, None)

        for line in infile:
            if not line or line.startswith('#'):
                continue

            fields = line.replace('"', '').strip().split('\t')

            if len(fields) > 4:
                motif_regex[fields[1]] = fields[4]

    return motif_regex


def parse_motif_domain_table(path: Path) -> dict[str, tuple[str, ...]]:
    """Parse an ELM-style motif-domain table."""

    motif_domains: dict[str, list[str]] = {}

    with open(path, encoding = 'utf-8') as infile:
        next(infile, None)

        for line in infile:
            fields = line.replace('"', '').strip().split('\t')

            if len(fields) <= 1:
                continue

            motif_domains.setdefault(fields[0], []).append(fields[1])

    return {
        motif: tuple(dict.fromkeys(domains))
        for motif, domains in motif_domains.items()
    }


def load_dmi_resources(repo_root: Path) -> DMIResources:
    """Load the merged ELM + 3did DMI resources."""

    resource_dir = repo_root.joinpath('microbiolink_api', 'resources')

    resource_specs = [
        (
            'ELM',
            resource_dir.joinpath('elm_classes.tsv'),
            resource_dir.joinpath('elm_interaction_domains.tsv'),
        ),
        (
            '3did',
            resource_dir.joinpath('3did_dmi_classes.tsv'),
            resource_dir.joinpath('3did_dmi_interaction_domains.tsv'),
        ),
    ]

    motif_patterns: dict[str, re.Pattern[str]] = {}
    motif_domains: dict[str, tuple[str, ...]] = {}
    motif_sources: dict[str, str] = {}

    for source_name, regex_path, motif_domain_path in resource_specs:
        regex_table = parse_elm_regex_table(regex_path)
        domain_table = parse_motif_domain_table(motif_domain_path)

        for motif_name, motif_pattern in regex_table.items():
            existing_pattern = motif_patterns.get(motif_name)

            if existing_pattern is not None and existing_pattern.pattern != motif_pattern:
                raise ValueError(
                    f'Conflicting motif pattern for {motif_name}: '
                    f'{existing_pattern.pattern} vs {motif_pattern}',
                )

            motif_patterns[motif_name] = re.compile(motif_pattern)
            motif_sources[motif_name] = source_name

        for motif_name, domains in domain_table.items():
            merged_domains = [*motif_domains.get(motif_name, ()), *domains]
            motif_domains[motif_name] = tuple(dict.fromkeys(merged_domains))
            motif_sources.setdefault(motif_name, source_name)

    return DMIResources(
        motif_patterns = motif_patterns,
        motif_domains = motif_domains,
        motif_sources = motif_sources,
    )


def build_domain_to_proteins(
    annotations: dict[str, UniProtAnnotation],
) -> dict[str, list[str]]:
    """Convert UniProt annotations into a domain-to-proteins mapping."""

    mapping: dict[str, list[str]] = {}

    for accession in sorted(annotations):
        for pfam_domain in annotations[accession].pfams:
            mapping.setdefault(pfam_domain, []).append(accession)

    return mapping


def scan_motif_presence(
    sequences: dict[str, str],
    motif_patterns: dict[str, re.Pattern[str]],
    index_map: dict[str, int],
) -> dict[str, dict[str, object]]:
    """Scan proteins for motif presence and summarize matches."""

    summary: dict[str, dict[str, object]] = {}

    for protein_id, sequence in sequences.items():
        bit = 1 << index_map[protein_id]

        for motif_name, motif_pattern in motif_patterns.items():
            count = 0
            first_start = -1
            first_end = -1

            for match in motif_pattern.finditer(sequence):
                count += 1

                if first_start == -1:
                    first_start = match.start()
                    first_end = match.end()

            if count == 0:
                continue

            motif_entry = summary.setdefault(
                motif_name,
                {
                    'mask': 0,
                    'proteins': {},
                    'total_matches': 0,
                },
            )
            motif_entry['mask'] |= bit
            motif_entry['total_matches'] += count
            motif_entry['proteins'][protein_id] = MotifPresence(
                count = count,
                first_start = first_start,
                first_end = first_end,
            )

    return summary


def build_true_masks(
    benchmark_pairs: Iterable[tuple[str, str]],
    human_index: dict[str, int],
) -> dict[str, int]:
    """Build human-partner bitmasks for each bacterial protein."""

    masks: dict[str, int] = {}

    for bacterial_protein, human_protein in benchmark_pairs:
        masks.setdefault(bacterial_protein, 0)
        masks[bacterial_protein] |= 1 << human_index[human_protein]

    return masks


def iter_mask_indices(mask: int) -> Iterator[int]:
    """Yield the set bit indices of a mask."""

    while mask:
        least_significant_bit = mask & -mask
        yield least_significant_bit.bit_length() - 1
        mask ^= least_significant_bit


def count_set_bits(mask: int) -> int:
    """Count the number of set bits in an integer mask."""

    return bin(mask).count('1')


def count_analyzable_pairs(
    benchmark_pairs: Iterable[tuple[str, str]],
    allowed_bacterial: set[str],
    allowed_human: set[str],
) -> int:
    """Count benchmark pairs analyzable by a method."""

    return sum(
        1
        for bacterial_protein, human_protein in benchmark_pairs
        if bacterial_protein in allowed_bacterial and human_protein in allowed_human
    )


def benchmark_forward_dmi(
    benchmark_pairs: dict[tuple[str, str], BenchmarkPair],
    human_annotations: dict[str, UniProtAnnotation],
    bacterial_annotations: dict[str, UniProtAnnotation],
    dmi_resources: DMIResources,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Benchmark forward host-motif vs bacterial-domain interactions."""

    human_proteins = sorted({pair.human_accession for pair in benchmark_pairs.values()})
    human_index = {protein: index for index, protein in enumerate(human_proteins)}
    index_to_human = {index: protein for protein, index in human_index.items()}

    true_masks = build_true_masks(benchmark_pairs, human_index)
    predicted_masks: dict[str, int] = {}

    human_sequences = {
        accession: annotation.sequence
        for accession, annotation in human_annotations.items()
        if annotation.sequence
    }
    bacterial_by_domain = build_domain_to_proteins(bacterial_annotations)
    motif_presence = scan_motif_presence(
        sequences = human_sequences,
        motif_patterns = dmi_resources.motif_patterns,
        index_map = human_index,
    )

    raw_prediction_rows = 0
    details: list[dict[str, object]] = []

    for motif_name, motif_entry in motif_presence.items():
        human_mask = int(motif_entry['mask'])
        human_support = motif_entry['proteins']
        total_matches = int(motif_entry['total_matches'])

        for domain_name in dmi_resources.motif_domains.get(motif_name, ()):
            bacterial_partners = bacterial_by_domain.get(domain_name, [])

            if not bacterial_partners:
                continue

            raw_prediction_rows += len(bacterial_partners) * total_matches

            for bacterial_protein in bacterial_partners:
                predicted_masks[bacterial_protein] = (
                    predicted_masks.get(bacterial_protein, 0) | human_mask
                )

                tp_mask = human_mask & true_masks.get(bacterial_protein, 0)

                if tp_mask == 0:
                    continue

                for human_index_value in iter_mask_indices(tp_mask):
                    human_protein = index_to_human[human_index_value]
                    motif_support = human_support[human_protein]
                    details.append(
                        {
                            'direction': 'forward',
                            'bacterial_accession': bacterial_protein,
                            'human_accession': human_protein,
                            'motif': motif_name,
                            'domain': domain_name,
                            'resource': dmi_resources.motif_sources.get(motif_name, 'custom'),
                            'motif_match_count': motif_support.count,
                            'first_start': motif_support.first_start,
                            'first_end': motif_support.first_end,
                        },
                    )

    unique_predicted_pairs = sum(count_set_bits(mask) for mask in predicted_masks.values())
    true_positives_recovered = sum(
        count_set_bits(predicted_masks.get(bacterial_protein, 0) & true_mask)
        for bacterial_protein, true_mask in true_masks.items()
    )

    analyzable_true_pairs = count_analyzable_pairs(
        benchmark_pairs = benchmark_pairs,
        allowed_bacterial = set(bacterial_annotations),
        allowed_human = set(human_sequences),
    )

    summary = {
        'approach': 'forward_dmi',
        'resource_name': 'elm_plus_3did',
        'benchmark_unique_pairs': len(benchmark_pairs),
        'analyzable_true_pairs': analyzable_true_pairs,
        'unique_predicted_pairs': unique_predicted_pairs,
        'raw_prediction_rows': raw_prediction_rows,
        'true_positives_recovered': true_positives_recovered,
        'false_positive_pairs_within_panel_cross_product': (
            unique_predicted_pairs - true_positives_recovered
        ),
        'recall_total': true_positives_recovered / len(benchmark_pairs),
        'recall_on_analyzable_pairs': (
            true_positives_recovered / analyzable_true_pairs
            if analyzable_true_pairs
            else 0.0
        ),
        'precision_unique_pairs': (
            true_positives_recovered / unique_predicted_pairs
            if unique_predicted_pairs
            else 0.0
        ),
    }

    return summary, details


def benchmark_reverse_dmi(
    benchmark_pairs: dict[tuple[str, str], BenchmarkPair],
    human_annotations: dict[str, UniProtAnnotation],
    bacterial_annotations: dict[str, UniProtAnnotation],
    dmi_resources: DMIResources,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Benchmark reverse bacterial-motif vs host-domain interactions."""

    human_proteins = sorted({pair.human_accession for pair in benchmark_pairs.values()})
    human_index = {protein: index for index, protein in enumerate(human_proteins)}
    index_to_human = {index: protein for protein, index in human_index.items()}

    true_masks = build_true_masks(benchmark_pairs, human_index)
    predicted_masks: dict[str, int] = {}

    bacterial_sequences = {
        accession: annotation.sequence
        for accession, annotation in bacterial_annotations.items()
        if annotation.sequence
    }
    bacterial_index = {
        protein: index
        for index, protein in enumerate(sorted(bacterial_sequences))
    }
    bacterial_motif_presence = scan_motif_presence(
        sequences = bacterial_sequences,
        motif_patterns = dmi_resources.motif_patterns,
        index_map = bacterial_index,
    )

    human_by_domain = build_domain_to_proteins(human_annotations)
    human_mask_by_domain = {
        domain_name: sum(1 << human_index[protein] for protein in proteins)
        for domain_name, proteins in human_by_domain.items()
        if proteins
    }

    raw_prediction_rows = 0
    details: list[dict[str, object]] = []

    for motif_name, motif_entry in bacterial_motif_presence.items():
        bacterial_support = motif_entry['proteins']
        bacterial_proteins = sorted(bacterial_support)
        total_matches = int(motif_entry['total_matches'])

        for domain_name in dmi_resources.motif_domains.get(motif_name, ()):
            human_mask = human_mask_by_domain.get(domain_name, 0)
            human_partners = human_by_domain.get(domain_name, [])

            if human_mask == 0 or not human_partners:
                continue

            raw_prediction_rows += len(human_partners) * total_matches

            for bacterial_protein in bacterial_proteins:
                predicted_masks[bacterial_protein] = (
                    predicted_masks.get(bacterial_protein, 0) | human_mask
                )

                tp_mask = human_mask & true_masks.get(bacterial_protein, 0)

                if tp_mask == 0:
                    continue

                motif_support = bacterial_support[bacterial_protein]

                for human_index_value in iter_mask_indices(tp_mask):
                    details.append(
                        {
                            'direction': 'reverse',
                            'bacterial_accession': bacterial_protein,
                            'human_accession': index_to_human[human_index_value],
                            'motif': motif_name,
                            'domain': domain_name,
                            'resource': dmi_resources.motif_sources.get(motif_name, 'custom'),
                            'motif_match_count': motif_support.count,
                            'first_start': motif_support.first_start,
                            'first_end': motif_support.first_end,
                        },
                    )

    unique_predicted_pairs = sum(count_set_bits(mask) for mask in predicted_masks.values())
    true_positives_recovered = sum(
        count_set_bits(predicted_masks.get(bacterial_protein, 0) & true_mask)
        for bacterial_protein, true_mask in true_masks.items()
    )

    analyzable_true_pairs = count_analyzable_pairs(
        benchmark_pairs = benchmark_pairs,
        allowed_bacterial = set(bacterial_sequences),
        allowed_human = {
            accession
            for accession, annotation in human_annotations.items()
            if annotation.pfams
        },
    )

    summary = {
        'approach': 'reverse_dmi',
        'resource_name': 'elm_plus_3did',
        'benchmark_unique_pairs': len(benchmark_pairs),
        'analyzable_true_pairs': analyzable_true_pairs,
        'unique_predicted_pairs': unique_predicted_pairs,
        'raw_prediction_rows': raw_prediction_rows,
        'true_positives_recovered': true_positives_recovered,
        'false_positive_pairs_within_panel_cross_product': (
            unique_predicted_pairs - true_positives_recovered
        ),
        'recall_total': true_positives_recovered / len(benchmark_pairs),
        'recall_on_analyzable_pairs': (
            true_positives_recovered / analyzable_true_pairs
            if analyzable_true_pairs
            else 0.0
        ),
        'precision_unique_pairs': (
            true_positives_recovered / unique_predicted_pairs
            if unique_predicted_pairs
            else 0.0
        ),
    }

    return summary, details


def load_undirected_pair_resource(path: Path) -> set[tuple[str, str]]:
    """Load a two-column undirected Pfam-Pfam resource."""

    pairs: set[tuple[str, str]] = set()

    with open(path, encoding = 'utf-8') as infile:
        for line in infile:
            fields = line.strip().split()

            if len(fields) < 2:
                continue

            pairs.add(tuple(sorted((fields[0], fields[1]))))

    return pairs


def load_domine_resources(path: Path) -> tuple[set[tuple[str, str]], set[tuple[str, str]]]:
    """Load DOMINE all-pair and high-confidence pair sets."""

    all_pairs: set[tuple[str, str]] = set()
    hc_pairs: set[tuple[str, str]] = set()

    with open(path, encoding = 'utf-8') as infile:
        for line in infile:
            fields = line.strip().split('|')

            if len(fields) < 18:
                continue

            pair = tuple(sorted((fields[0], fields[1])))
            all_pairs.add(pair)

            if fields[17] == 'HC':
                hc_pairs.add(pair)

    return all_pairs, hc_pairs


def benchmark_ddi(
    benchmark_pairs: dict[tuple[str, str], BenchmarkPair],
    human_annotations: dict[str, UniProtAnnotation],
    bacterial_annotations: dict[str, UniProtAnnotation],
    resource_name: str,
    resource_pairs: set[tuple[str, str]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Benchmark one DDI resource against the PHISTO protein panel."""

    human_proteins = sorted({pair.human_accession for pair in benchmark_pairs.values()})
    human_index = {protein: index for index, protein in enumerate(human_proteins)}
    index_to_human = {index: protein for protein, index in human_index.items()}

    true_masks = build_true_masks(benchmark_pairs, human_index)
    predicted_masks: dict[str, int] = {}
    details: list[dict[str, object]] = []
    raw_prediction_rows = 0

    bacterial_by_domain = build_domain_to_proteins(bacterial_annotations)
    human_by_domain = build_domain_to_proteins(human_annotations)
    human_mask_by_domain = {
        domain_name: sum(1 << human_index[protein] for protein in proteins)
        for domain_name, proteins in human_by_domain.items()
        if proteins
    }

    for domain_a, domain_b in sorted(resource_pairs):
        direction_specs = [(domain_a, domain_b)]

        if domain_a != domain_b:
            direction_specs.append((domain_b, domain_a))

        for bacterial_domain, human_domain in direction_specs:
            bacterial_partners = bacterial_by_domain.get(bacterial_domain, [])
            human_partners = human_by_domain.get(human_domain, [])
            human_mask = human_mask_by_domain.get(human_domain, 0)

            if not bacterial_partners or not human_partners:
                continue

            raw_prediction_rows += len(bacterial_partners) * len(human_partners)

            for bacterial_protein in bacterial_partners:
                predicted_masks[bacterial_protein] = (
                    predicted_masks.get(bacterial_protein, 0) | human_mask
                )

                tp_mask = human_mask & true_masks.get(bacterial_protein, 0)

                if tp_mask == 0:
                    continue

                for human_index_value in iter_mask_indices(tp_mask):
                    details.append(
                        {
                            'resource_name': resource_name,
                            'bacterial_accession': bacterial_protein,
                            'human_accession': index_to_human[human_index_value],
                            'bacterial_domain': bacterial_domain,
                            'human_domain': human_domain,
                        },
                    )

    unique_predicted_pairs = sum(count_set_bits(mask) for mask in predicted_masks.values())
    true_positives_recovered = sum(
        count_set_bits(predicted_masks.get(bacterial_protein, 0) & true_mask)
        for bacterial_protein, true_mask in true_masks.items()
    )

    analyzable_true_pairs = count_analyzable_pairs(
        benchmark_pairs = benchmark_pairs,
        allowed_bacterial = {
            accession
            for accession, annotation in bacterial_annotations.items()
            if annotation.pfams
        },
        allowed_human = {
            accession
            for accession, annotation in human_annotations.items()
            if annotation.pfams
        },
    )

    summary = {
        'approach': 'ddi',
        'resource_name': resource_name,
        'benchmark_unique_pairs': len(benchmark_pairs),
        'analyzable_true_pairs': analyzable_true_pairs,
        'resource_rows': len(resource_pairs),
        'unique_predicted_pairs': unique_predicted_pairs,
        'raw_prediction_rows': raw_prediction_rows,
        'true_positives_recovered': true_positives_recovered,
        'false_positive_pairs_within_panel_cross_product': (
            unique_predicted_pairs - true_positives_recovered
        ),
        'recall_total': true_positives_recovered / len(benchmark_pairs),
        'recall_on_analyzable_pairs': (
            true_positives_recovered / analyzable_true_pairs
            if analyzable_true_pairs
            else 0.0
        ),
        'precision_unique_pairs': (
            true_positives_recovered / unique_predicted_pairs
            if unique_predicted_pairs
            else 0.0
        ),
    }

    return summary, details


def write_benchmark_pair_table(
    path: Path,
    benchmark_pairs: dict[tuple[str, str], BenchmarkPair],
) -> None:
    """Write the unique PHISTO pair panel."""

    rows = []

    for pair in benchmark_pairs.values():
        rows.append(
            {
                'pathogen_accession': pair.pathogen_accession,
                'human_accession': pair.human_accession,
                'pathogen_name': pair.pathogen_name,
                'pathogen_protein_name': pair.pathogen_protein_name,
                'human_protein_name': pair.human_protein_name,
                'supporting_rows': pair.supporting_rows,
                'taxonomy_ids': ';'.join(sorted(tax_id for tax_id in pair.taxonomy_ids if tax_id)),
                'methods': ';'.join(sorted(method for method in pair.methods if method)),
                'pubmed_ids': ';'.join(sorted(pmid for pmid in pair.pubmed_ids if pmid)),
            },
        )

    write_tsv(
        path,
        rows,
        [
            'pathogen_accession',
            'human_accession',
            'pathogen_name',
            'pathogen_protein_name',
            'human_protein_name',
            'supporting_rows',
            'taxonomy_ids',
            'methods',
            'pubmed_ids',
        ],
    )


def write_frequency_table(
    path: Path,
    values: Iterable[str],
    column_name: str,
) -> None:
    """Write simple value counts to TSV."""

    counts: dict[str, int] = {}

    for value in values:
        counts[value] = counts.get(value, 0) + 1

    rows = [
        {
            column_name: key,
            'count': counts[key],
        }
        for key in sorted(counts, key = lambda item: (-counts[item], item))
    ]

    write_tsv(path, rows, [column_name, 'count'])


def main() -> None:
    """Run the PHISTO benchmark."""

    args = parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    phisto_csv = Path(args.phisto_csv).resolve()
    domine_interaction_file = Path(args.domine_interaction_file).resolve()
    output_dir = Path(args.output_dir).resolve()
    raw_dir = output_dir.joinpath('raw')
    inputs_dir = output_dir.joinpath('inputs')

    ensure_directory(output_dir)
    ensure_directory(raw_dir)
    ensure_directory(inputs_dir)

    phisto_rows = read_phisto_rows(phisto_csv)
    benchmark_pairs = aggregate_benchmark_pairs(phisto_rows)
    bacterial_accessions = sorted(
        {
            pair.pathogen_accession
            for pair in benchmark_pairs.values()
        },
    )
    human_accessions = sorted(
        {
            pair.human_accession
            for pair in benchmark_pairs.values()
        },
    )

    write_benchmark_pair_table(
        output_dir.joinpath('phisto_unique_pairs.tsv'),
        benchmark_pairs,
    )
    write_frequency_table(
        output_dir.joinpath('phisto_method_counts.tsv'),
        (row.experimental_method for row in phisto_rows),
        'experimental_method',
    )
    write_frequency_table(
        output_dir.joinpath('phisto_pmid_counts.tsv'),
        (row.pubmed_id for row in phisto_rows),
        'pubmed_id',
    )
    write_text_lines(output_dir.joinpath('bacterial_accessions.txt'), bacterial_accessions)
    write_text_lines(output_dir.joinpath('human_accessions.txt'), human_accessions)

    bacterial_annotations = resolve_uniprot_annotations(bacterial_accessions)
    human_annotations = resolve_uniprot_annotations(human_accessions)
    unresolved_bacterial = sorted(
        set(bacterial_accessions) - set(bacterial_annotations),
    )
    unresolved_human = sorted(set(human_accessions) - set(human_annotations))

    write_text_lines(
        output_dir.joinpath('unresolved_bacterial_accessions.txt'),
        unresolved_bacterial,
    )
    write_text_lines(
        output_dir.joinpath('unresolved_human_accessions.txt'),
        unresolved_human,
    )

    write_annotation_table(
        output_dir.joinpath('bacterial_uniprot_annotations.tsv'),
        bacterial_annotations,
    )
    write_annotation_table(
        output_dir.joinpath('human_uniprot_annotations.tsv'),
        human_annotations,
    )

    write_fasta(inputs_dir.joinpath('bacterial_sequences.fasta'), bacterial_annotations)
    write_fasta(inputs_dir.joinpath('human_sequences.fasta'), human_annotations)
    write_domain_table(inputs_dir.joinpath('bacterial_domains.tsv'), bacterial_annotations)
    write_domain_table(inputs_dir.joinpath('human_domains.tsv'), human_annotations)

    dmi_resources = load_dmi_resources(repo_root)

    forward_summary, forward_details = benchmark_forward_dmi(
        benchmark_pairs = benchmark_pairs,
        human_annotations = human_annotations,
        bacterial_annotations = bacterial_annotations,
        dmi_resources = dmi_resources,
    )
    reverse_summary, reverse_details = benchmark_reverse_dmi(
        benchmark_pairs = benchmark_pairs,
        human_annotations = human_annotations,
        bacterial_annotations = bacterial_annotations,
        dmi_resources = dmi_resources,
    )

    write_tsv(
        output_dir.joinpath('forward_dmi_true_positive_details.tsv'),
        forward_details,
        [
            'direction',
            'bacterial_accession',
            'human_accession',
            'motif',
            'domain',
            'resource',
            'motif_match_count',
            'first_start',
            'first_end',
        ],
    )
    write_tsv(
        output_dir.joinpath('reverse_dmi_true_positive_details.tsv'),
        reverse_details,
        [
            'direction',
            'bacterial_accession',
            'human_accession',
            'motif',
            'domain',
            'resource',
            'motif_match_count',
            'first_start',
            'first_end',
        ],
    )

    ddi_3did = load_undirected_pair_resource(
        repo_root.joinpath('DDI', 'resources', 'pfam_interactions_3did_current.tsv'),
    )
    domine_all, domine_hc = load_domine_resources(domine_interaction_file)

    ddi_resource_sets = [
        ('3did_current', ddi_3did),
        ('domine_v2_hc', domine_hc),
        ('domine_v2_all', domine_all),
        ('3did_plus_domine_v2_hc', ddi_3did | domine_hc),
        ('3did_plus_domine_v2_all', ddi_3did | domine_all),
    ]

    ddi_summaries: list[dict[str, object]] = []
    ddi_details_with_resource: list[dict[str, object]] = []

    for resource_name, resource_pairs in ddi_resource_sets:
        summary, details = benchmark_ddi(
            benchmark_pairs = benchmark_pairs,
            human_annotations = human_annotations,
            bacterial_annotations = bacterial_annotations,
            resource_name = resource_name,
            resource_pairs = resource_pairs,
        )
        ddi_summaries.append(summary)
        ddi_details_with_resource.extend(details)

    write_tsv(
        output_dir.joinpath('ddi_true_positive_details.tsv'),
        ddi_details_with_resource,
        [
            'resource_name',
            'bacterial_accession',
            'human_accession',
            'bacterial_domain',
            'human_domain',
        ],
    )

    benchmark_summary_rows = [
        {
            'metric': 'phisto_rows',
            'value': len(phisto_rows),
        },
        {
            'metric': 'phisto_unique_pairs',
            'value': len(benchmark_pairs),
        },
        {
            'metric': 'unique_bacterial_accessions',
            'value': len(bacterial_accessions),
        },
        {
            'metric': 'unique_human_accessions',
            'value': len(human_accessions),
        },
        {
            'metric': 'resolved_bacterial_annotations',
            'value': len(bacterial_annotations),
        },
        {
            'metric': 'resolved_human_annotations',
            'value': len(human_annotations),
        },
        {
            'metric': 'bacterial_with_pfam',
            'value': sum(1 for annotation in bacterial_annotations.values() if annotation.pfams),
        },
        {
            'metric': 'human_with_pfam',
            'value': sum(1 for annotation in human_annotations.values() if annotation.pfams),
        },
        {
            'metric': 'unresolved_bacterial_accessions',
            'value': len(unresolved_bacterial),
        },
        {
            'metric': 'unresolved_human_accessions',
            'value': len(unresolved_human),
        },
        {
            'metric': 'bacterial_with_sequence',
            'value': sum(1 for annotation in bacterial_annotations.values() if annotation.sequence),
        },
        {
            'metric': 'human_with_sequence',
            'value': sum(1 for annotation in human_annotations.values() if annotation.sequence),
        },
    ]

    write_tsv(
        output_dir.joinpath('benchmark_overview.tsv'),
        benchmark_summary_rows,
        ['metric', 'value'],
    )

    approach_summaries = [forward_summary, reverse_summary, *ddi_summaries]
    write_tsv(
        output_dir.joinpath('approach_summary.tsv'),
        approach_summaries,
        [
            'approach',
            'resource_name',
            'benchmark_unique_pairs',
            'analyzable_true_pairs',
            'resource_rows',
            'unique_predicted_pairs',
            'raw_prediction_rows',
            'true_positives_recovered',
            'false_positive_pairs_within_panel_cross_product',
            'recall_total',
            'recall_on_analyzable_pairs',
            'precision_unique_pairs',
        ],
    )

    print(
        'PHISTO benchmark complete\n'
        f'- Unique PHISTO pairs: {len(benchmark_pairs)}\n'
        f'- Bacterial accessions: {len(bacterial_accessions)} '
        f'({len(bacterial_annotations)} resolved)\n'
        f'- Human accessions: {len(human_accessions)} '
        f'({len(human_annotations)} resolved)\n'
        f'- Forward DMI true positives: '
        f"{forward_summary['true_positives_recovered']}\n"
        f'- Reverse DMI true positives: '
        f"{reverse_summary['true_positives_recovered']}\n"
        f'- Best DDI resource by recall: '
        f"{max(ddi_summaries, key = lambda item: item['recall_total'])['resource_name']}\n"
        f'- Output directory: {output_dir}'
    )


if __name__ == '__main__':
    main()
