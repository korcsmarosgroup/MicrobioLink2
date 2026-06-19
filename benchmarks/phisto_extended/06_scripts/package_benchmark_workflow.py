#!/usr/bin/env python3

"""
Run the PHISTO benchmark through the packaged MicrobioLink 2.1 API.

The original documented bundle in this benchmark folder preserves the exact
standard-library script and outputs that were produced during the initial
analysis. This helper module adds a more interactive, package-driven route for
recreating the same upstream benchmark logic from a notebook:

1. build the unique PHISTO true-positive panel,
2. export accession lists,
3. optionally redownload FASTA and Pfam annotations with the packaged generic
   UniProt helpers,
4. rerun forward DMI, reverse DMI, and DDI with the packaged resources,
5. write compact benchmark summaries and recovered true-positive detail tables.

The prediction steps are batched for forward and reverse DMI so the full PHISTO
panel can be rerun through `microbiolink_api` without keeping every predicted
interaction row in memory at once.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
from typing import Iterator
from typing import Mapping

import pandas as pd

from microbiolink.download_protein_domains import (
    DEFAULT_UNIPROT_FIELDS,
    UNIPROT_BATCH_SIZE,
    download_protein_list_with_fields,
)
from microbiolink.get_protein_fasta import DEFAULT_BATCH_SIZE, fetch_fasta_sequences
from microbiolink_api import (
    bidirectional_interactions_to_dataframe,
    ddi_interactions_to_dataframe,
    extract_uniprot_id,
    interactions_to_dataframe,
    load_default_3did_ddi_resource_bundle,
    load_default_ddi_resource_bundle,
    load_default_dmi_resource_bundle,
    load_default_domine_all_ddi_resource_bundle,
    load_default_domine_hc_ddi_resource_bundle,
    merge_ddi_resource_bundles,
    predict_domain_domain_interactions_from_data,
    predict_domain_motif_interactions_from_data,
    predict_reverse_domain_motif_interactions_from_data,
    read_fasta_sequences,
    read_protein_domain_table,
)


PHISTO_COLUMNS = [
    'pathogen_name',
    'taxonomy_id',
    'pathogen_accession',
    'pathogen_protein_name',
    'human_accession',
    'human_protein_name',
    'experimental_method',
    'pubmed_id',
]
FORWARD_DMI_RESOURCE_NAME = 'elm_plus_3did'
DEFAULT_FORWARD_BATCH_SIZE = 100
DEFAULT_REVERSE_BATCH_SIZE = 25


@dataclass(frozen = True)
class BenchmarkInputPaths:
    """Collect the canonical file paths inside the documented benchmark bundle."""

    phisto_csv: Path
    bacterial_fasta: Path
    human_fasta: Path
    bacterial_domains: Path
    human_domains: Path
    bacterial_annotations: Path | None = None
    human_annotations: Path | None = None


def resolve_benchmark_input_paths(
    benchmark_root: str | Path,
) -> BenchmarkInputPaths:
    """
    Resolve the main input files used by the documented PHISTO bundle.

    Args:
        benchmark_root: Root directory of the benchmark bundle.

    Returns:
        Structured paths for the raw PHISTO export and prepared inputs.
    """

    root = Path(benchmark_root)

    return BenchmarkInputPaths(
        phisto_csv = root / '01_raw_sources' / 'phisto_bacteria.csv',
        bacterial_fasta = (
            root / '04_benchmark_inputs' / 'bacterial_sequences.fasta'
        ),
        human_fasta = root / '04_benchmark_inputs' / 'human_sequences.fasta',
        bacterial_domains = root / '04_benchmark_inputs' / 'bacterial_domains.tsv',
        human_domains = root / '04_benchmark_inputs' / 'human_domains.tsv',
        bacterial_annotations = (
            root
            / '03_resolved_annotations'
            / 'bacterial_uniprot_annotations.tsv'
        ),
        human_annotations = (
            root
            / '03_resolved_annotations'
            / 'human_uniprot_annotations.tsv'
        ),
    )


def read_phisto_export(
    phisto_csv: str | Path,
) -> pd.DataFrame:
    """
    Read the raw PHISTO bacteria-human export.

    Args:
        phisto_csv: Path to the PHISTO CSV export.

    Returns:
        Data frame with the canonical PHISTO benchmark columns.
    """

    phisto_frame = pd.read_csv(
        phisto_csv,
        encoding = 'latin-1',
        header = 0,
        names = PHISTO_COLUMNS,
    )

    return phisto_frame.fillna('')


def _join_unique(values: Iterable[object]) -> str:
    """Join unique non-empty values while preserving input order."""

    ordered_values = OrderedDict()

    for value in values:
        text_value = str(value).strip()

        if not text_value:
            continue

        ordered_values[text_value] = None

    return ';'.join(ordered_values)


def build_unique_pair_panel(
    phisto_csv: str | Path,
) -> pd.DataFrame:
    """
    Collapse the raw PHISTO export into unique pathogen-human pairs.

    Args:
        phisto_csv: Path to the PHISTO CSV export.

    Returns:
        Deduplicated true-positive panel matching the documented benchmark
        layout.
    """

    phisto_frame = read_phisto_export(phisto_csv)
    phisto_frame = phisto_frame[
        (phisto_frame['pathogen_accession'].astype(str).str.strip() != '')
        & (phisto_frame['human_accession'].astype(str).str.strip() != '')
    ].copy()

    grouped_frame = (
        phisto_frame
        .groupby(
            ['pathogen_accession', 'human_accession'],
            as_index = False,
            sort = False,
        )
        .agg(
            pathogen_name = ('pathogen_name', 'first'),
            pathogen_protein_name = ('pathogen_protein_name', 'first'),
            human_protein_name = ('human_protein_name', 'first'),
            supporting_rows = ('pathogen_accession', 'size'),
            taxonomy_ids = ('taxonomy_id', _join_unique),
            methods = ('experimental_method', _join_unique),
            pubmed_ids = ('pubmed_id', _join_unique),
        )
    )

    ordered_columns = [
        'pathogen_accession',
        'human_accession',
        'pathogen_name',
        'pathogen_protein_name',
        'human_protein_name',
        'supporting_rows',
        'taxonomy_ids',
        'methods',
        'pubmed_ids',
    ]

    return grouped_frame[ordered_columns]


def write_accession_lists(
    unique_pairs: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """
    Write one-column bacterial and human accession tables for reuse.

    Args:
        unique_pairs: Deduplicated PHISTO benchmark panel.
        output_dir: Directory where accession tables will be written.

    Returns:
        Mapping from logical accession list name to output path.
    """

    output_path = Path(output_dir)
    output_path.mkdir(parents = True, exist_ok = True)

    bacterial_accessions = list(
        OrderedDict.fromkeys(unique_pairs['pathogen_accession'].tolist()),
    )
    human_accessions = list(
        OrderedDict.fromkeys(unique_pairs['human_accession'].tolist()),
    )

    bacterial_path = output_path / 'bacterial_accessions.tsv'
    human_path = output_path / 'human_accessions.tsv'

    pd.DataFrame({'accession': bacterial_accessions}).to_csv(
        bacterial_path,
        sep = '\t',
        index = False,
    )
    pd.DataFrame({'accession': human_accessions}).to_csv(
        human_path,
        sep = '\t',
        index = False,
    )

    return {
        'bacterial': bacterial_path,
        'human': human_path,
    }


def _unique_preserve_order(
    values: Iterable[str],
) -> list[str]:
    """Return unique strings while preserving order."""

    return list(OrderedDict.fromkeys(values))


def download_fasta_from_accessions(
    accessions: Iterable[str],
    output_file: str | Path,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> Path:
    """
    Download FASTA sequences for a list of UniProt accessions.

    Args:
        accessions: UniProt accession strings.
        output_file: Path to the FASTA file to be written.
        batch_size: Number of accessions per UniProt request.

    Returns:
        Path to the written FASTA file.
    """

    output_path = Path(output_file)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    ordered_accessions = _unique_preserve_order(
        [
            accession
            for accession in accessions
            if str(accession).strip()
        ],
    )

    with open(output_path, 'w', encoding = 'utf-8') as fasta_file:
        for batch_start in range(0, len(ordered_accessions), batch_size):
            batch = ordered_accessions[batch_start : batch_start + batch_size]
            fasta_file.write(fetch_fasta_sequences(batch))

    return output_path


def download_domain_annotations_from_accessions(
    accessions: Iterable[str],
    output_file: str | Path,
    batch_size: int = UNIPROT_BATCH_SIZE,
    fields: list[str] | None = None,
) -> Path:
    """
    Download UniProt domain annotation tables for a list of accessions.

    Args:
        accessions: UniProt accession strings.
        output_file: Path to the raw UniProt TSV file to be written.
        batch_size: Number of accessions per UniProt request.
        fields: UniProt field list. Defaults to the packaged helper fields.

    Returns:
        Path to the written raw TSV file.
    """

    output_path = Path(output_file)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    ordered_accessions = _unique_preserve_order(
        [
            accession
            for accession in accessions
            if str(accession).strip()
        ],
    )
    selected_fields = fields or list(DEFAULT_UNIPROT_FIELDS)
    header_written = False

    with open(output_path, 'w', encoding = 'utf-8') as domain_file:
        for batch_start in range(0, len(ordered_accessions), batch_size):
            batch = ordered_accessions[batch_start : batch_start + batch_size]
            response_text = download_protein_list_with_fields(
                batch,
                fields = selected_fields,
            )
            lines = response_text.rstrip().splitlines()

            if not lines:
                continue

            if not header_written:
                domain_file.write(lines[0] + '\n')
                header_written = True

            for line in lines[1:]:
                if line.strip():
                    domain_file.write(line + '\n')

    return output_path


def _find_first_matching_column(
    columns: list[str],
    candidate_names: Iterable[str],
) -> str:
    """Return the first case-insensitive matching column name."""

    lower_to_original = {
        column_name.lower(): column_name
        for column_name in columns
    }

    for candidate_name in candidate_names:
        matched_name = lower_to_original.get(candidate_name.lower())

        if matched_name is not None:
            return matched_name

    raise KeyError(
        'Could not find a matching column in the UniProt download table for '
        f'any of: {", ".join(candidate_names)}',
    )


def normalize_uniprot_domain_download(
    raw_domain_table: str | Path,
    output_file: str | Path,
) -> pd.DataFrame:
    """
    Convert a raw UniProt domain download into `protein`/`pfam` format.

    Args:
        raw_domain_table: Raw UniProt TSV produced by the generic helper.
        output_file: Path where the normalized MicrobioLink input table will
            be written.

    Returns:
        The normalized two-column protein-to-Pfam data frame.
    """

    raw_frame = pd.read_csv(raw_domain_table, sep = '\t').fillna('')
    accession_column = _find_first_matching_column(
        list(raw_frame.columns),
        ['Entry', 'Accession', 'accession'],
    )
    pfam_column = _find_first_matching_column(
        list(raw_frame.columns),
        ['Pfam', 'Cross-reference (Pfam)', 'xref_pfam'],
    )

    normalized_records: list[dict[str, str]] = []

    for _, row in raw_frame.iterrows():
        protein_id = str(row[accession_column]).strip()
        pfam_value = str(row[pfam_column]).strip()

        if not protein_id or not pfam_value:
            continue

        pfam_values = _unique_preserve_order(
            [
                pfam_name.strip()
                for pfam_name in pfam_value.split(';')
                if pfam_name.strip()
            ],
        )

        if not pfam_values:
            continue

        normalized_records.append(
            {
                'protein': protein_id,
                'pfam': ';'.join(pfam_values),
            },
        )

    normalized_frame = pd.DataFrame(normalized_records)

    if not normalized_frame.empty:
        normalized_frame = (
            normalized_frame
            .drop_duplicates()
            .sort_values('protein')
            .reset_index(drop = True)
        )
    else:
        normalized_frame = pd.DataFrame(columns = ['protein', 'pfam'])

    output_path = Path(output_file)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    normalized_frame.to_csv(
        output_path,
        sep = '\t',
        index = False,
    )

    return normalized_frame


def download_panel_inputs(
    unique_pairs: pd.DataFrame,
    output_dir: str | Path,
) -> dict[str, Path]:
    """
    Regenerate the PHISTO benchmark FASTA and Pfam inputs via package helpers.

    Args:
        unique_pairs: Deduplicated PHISTO benchmark panel.
        output_dir: Root directory where accession lists, raw downloads, and
            normalized input files will be written.

    Returns:
        Mapping from logical input name to generated file path.
    """

    root = Path(output_dir)
    accessions_dir = root / '01_accessions'
    raw_dir = root / '02_raw_uniprot_downloads'
    inputs_dir = root / '04_benchmark_inputs'

    accession_paths = write_accession_lists(unique_pairs, accessions_dir)

    bacterial_accessions = unique_pairs['pathogen_accession'].tolist()
    human_accessions = unique_pairs['human_accession'].tolist()

    bacterial_fasta = download_fasta_from_accessions(
        bacterial_accessions,
        inputs_dir / 'bacterial_sequences.fasta',
    )
    human_fasta = download_fasta_from_accessions(
        human_accessions,
        inputs_dir / 'human_sequences.fasta',
    )
    bacterial_domains_raw = download_domain_annotations_from_accessions(
        bacterial_accessions,
        raw_dir / 'bacterial_domains_raw.tsv',
    )
    human_domains_raw = download_domain_annotations_from_accessions(
        human_accessions,
        raw_dir / 'human_domains_raw.tsv',
    )
    bacterial_domains = normalize_uniprot_domain_download(
        bacterial_domains_raw,
        inputs_dir / 'bacterial_domains.tsv',
    )
    human_domains = normalize_uniprot_domain_download(
        human_domains_raw,
        inputs_dir / 'human_domains.tsv',
    )

    return {
        'bacterial_accessions': accession_paths['bacterial'],
        'human_accessions': accession_paths['human'],
        'bacterial_fasta': bacterial_fasta,
        'human_fasta': human_fasta,
        'bacterial_domains_raw': bacterial_domains_raw,
        'human_domains_raw': human_domains_raw,
        'bacterial_domains': inputs_dir / 'bacterial_domains.tsv',
        'human_domains': inputs_dir / 'human_domains.tsv',
    }


def _iter_sequence_batches(
    sequences: Mapping[str, str],
    batch_size: int,
) -> Iterator[dict[str, str]]:
    """Yield deterministic batches from a FASTA sequence mapping."""

    items = list(sequences.items())

    for batch_start in range(0, len(items), batch_size):
        yield dict(items[batch_start : batch_start + batch_size])


def _sequence_accessions(
    sequences: Mapping[str, str],
) -> set[str]:
    """Return the UniProt accessions represented by a FASTA mapping."""

    return {
        extract_uniprot_id(header)
        for header in sequences
    }


def _domain_proteins(
    protein_domains: Mapping[str, list[str]],
) -> set[str]:
    """Return the proteins represented by a Pfam-domain mapping."""

    return {
        protein_id
        for proteins in protein_domains.values()
        for protein_id in proteins
    }


def _subset_domain_mapping_by_proteins(
    protein_domains: Mapping[str, list[str]],
    allowed_proteins: set[str],
) -> dict[str, list[str]]:
    """Filter a Pfam-to-protein mapping down to a small protein subset."""

    subset_mapping: dict[str, list[str]] = {}

    for domain_name, proteins in protein_domains.items():
        kept_proteins = [
            protein_id
            for protein_id in proteins
            if protein_id in allowed_proteins
        ]

        if kept_proteins:
            subset_mapping[domain_name] = kept_proteins

    return subset_mapping


def _true_positive_pairs(
    unique_pairs: pd.DataFrame,
) -> set[tuple[str, str]]:
    """Return the benchmark true-positive pair set."""

    return set(
        zip(
            unique_pairs['pathogen_accession'],
            unique_pairs['human_accession'],
        ),
    )


def _count_analyzable_pairs(
    true_positive_pairs: set[tuple[str, str]],
    bacterial_proteins: set[str],
    human_proteins: set[str],
) -> int:
    """Count PHISTO pairs with all structural inputs required by a method."""

    return sum(
        1
        for bacterial_protein, human_protein in true_positive_pairs
        if (
            bacterial_protein in bacterial_proteins
            and human_protein in human_proteins
        )
    )


def _safe_divide(
    numerator: int,
    denominator: int,
) -> float:
    """Return a zero-safe division result."""

    if denominator == 0:
        return 0.0

    return numerator / denominator


def _write_frame(
    data_frame: pd.DataFrame,
    output_file: str | Path,
) -> Path:
    """Write a TSV data frame and return its path."""

    output_path = Path(output_file)
    output_path.parent.mkdir(parents = True, exist_ok = True)
    data_frame.to_csv(output_path, sep = '\t', index = False)
    return output_path


def _detail_frame_from_records(
    detail_records: list[dict[str, object]],
    columns: list[str],
) -> pd.DataFrame:
    """Create a deterministically ordered detail data frame."""

    if not detail_records:
        return pd.DataFrame(columns = columns)

    detail_frame = pd.DataFrame(detail_records, columns = columns)

    return detail_frame.sort_values(columns).reset_index(drop = True)


def run_forward_dmi_benchmark(
    unique_pairs: pd.DataFrame,
    human_fasta_file: str | Path,
    bacterial_domain_file: str | Path,
    batch_size: int = DEFAULT_FORWARD_BATCH_SIZE,
) -> tuple[dict[str, object], pd.DataFrame]:
    """
    Benchmark forward DMI on the PHISTO panel via the packaged API.

    Args:
        unique_pairs: Deduplicated PHISTO benchmark panel.
        human_fasta_file: Host FASTA input file.
        bacterial_domain_file: Bacterial protein-to-Pfam table.
        batch_size: Number of host proteins per prediction batch.

    Returns:
        Summary-row dictionary and recovered true-positive detail frame.
    """

    dmi_bundle = load_default_dmi_resource_bundle()
    human_sequences = read_fasta_sequences(human_fasta_file)
    bacterial_domains = read_protein_domain_table(bacterial_domain_file)
    true_positive_pairs = _true_positive_pairs(unique_pairs)
    analyzable_true_pairs = _count_analyzable_pairs(
        true_positive_pairs,
        bacterial_proteins = _domain_proteins(bacterial_domains),
        human_proteins = _sequence_accessions(human_sequences),
    )

    raw_prediction_rows = 0
    unique_predicted_pairs: set[tuple[str, str]] = set()
    recovered_pairs: set[tuple[str, str]] = set()
    detail_lookup: dict[
        tuple[str, str, str, str, str],
        dict[str, object],
    ] = {}

    for sequence_batch in _iter_sequence_batches(human_sequences, batch_size):
        interactions = predict_domain_motif_interactions_from_data(
            human_sequences = sequence_batch,
            elm_regex = dmi_bundle.elm_regex,
            motif_domains = dmi_bundle.motif_domains,
            bacterial_domains = bacterial_domains,
            motif_sources = dmi_bundle.motif_sources,
        )
        raw_prediction_rows += len(interactions)

        for interaction in interactions:
            pair = (
                interaction.bacterial_protein,
                interaction.human_protein,
            )
            unique_predicted_pairs.add(pair)

            if pair not in true_positive_pairs:
                continue

            recovered_pairs.add(pair)
            detail_key = (
                interaction.bacterial_protein,
                interaction.human_protein,
                interaction.motif,
                interaction.bacterial_domain,
                interaction.resource,
            )
            detail_row = detail_lookup.setdefault(
                detail_key,
                {
                    'direction': 'forward',
                    'bacterial_accession': interaction.bacterial_protein,
                    'human_accession': interaction.human_protein,
                    'motif': interaction.motif,
                    'domain': interaction.bacterial_domain,
                    'resource': interaction.resource,
                    'motif_match_count': 0,
                    'first_start': interaction.start,
                    'first_end': interaction.end,
                },
            )
            detail_row['motif_match_count'] += 1

            if interaction.start < int(detail_row['first_start']):
                detail_row['first_start'] = interaction.start
                detail_row['first_end'] = interaction.end

    true_positives_recovered = len(recovered_pairs)
    unique_predicted_pair_count = len(unique_predicted_pairs)
    false_positive_pairs = len(unique_predicted_pairs - true_positive_pairs)

    summary_row = {
        'approach': 'forward_dmi',
        'resource_name': FORWARD_DMI_RESOURCE_NAME,
        'benchmark_unique_pairs': len(true_positive_pairs),
        'analyzable_true_pairs': analyzable_true_pairs,
        'resource_rows': '',
        'unique_predicted_pairs': unique_predicted_pair_count,
        'raw_prediction_rows': raw_prediction_rows,
        'true_positives_recovered': true_positives_recovered,
        'false_positive_pairs_within_panel_cross_product': (
            false_positive_pairs
        ),
        'recall_total': _safe_divide(
            true_positives_recovered,
            len(true_positive_pairs),
        ),
        'recall_on_analyzable_pairs': _safe_divide(
            true_positives_recovered,
            analyzable_true_pairs,
        ),
        'precision_unique_pairs': _safe_divide(
            true_positives_recovered,
            unique_predicted_pair_count,
        ),
    }

    detail_frame = _detail_frame_from_records(
        list(detail_lookup.values()),
        columns = [
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

    return summary_row, detail_frame


def run_reverse_dmi_benchmark(
    unique_pairs: pd.DataFrame,
    bacterial_fasta_file: str | Path,
    human_domain_file: str | Path,
    batch_size: int = DEFAULT_REVERSE_BATCH_SIZE,
) -> tuple[dict[str, object], pd.DataFrame]:
    """
    Benchmark reverse DMI on the PHISTO panel via the packaged API.

    Args:
        unique_pairs: Deduplicated PHISTO benchmark panel.
        bacterial_fasta_file: Bacterial FASTA input file.
        human_domain_file: Human protein-to-Pfam table.
        batch_size: Number of bacterial proteins per prediction batch.

    Returns:
        Summary-row dictionary and recovered true-positive detail frame.
    """

    dmi_bundle = load_default_dmi_resource_bundle()
    bacterial_sequences = read_fasta_sequences(bacterial_fasta_file)
    human_domains = read_protein_domain_table(human_domain_file)
    true_positive_pairs = _true_positive_pairs(unique_pairs)
    analyzable_true_pairs = _count_analyzable_pairs(
        true_positive_pairs,
        bacterial_proteins = _sequence_accessions(bacterial_sequences),
        human_proteins = _domain_proteins(human_domains),
    )

    raw_prediction_rows = 0
    unique_predicted_pairs: set[tuple[str, str]] = set()
    recovered_pairs: set[tuple[str, str]] = set()
    detail_lookup: dict[
        tuple[str, str, str, str, str],
        dict[str, object],
    ] = {}

    for sequence_batch in _iter_sequence_batches(
        bacterial_sequences,
        batch_size,
    ):
        interactions = predict_reverse_domain_motif_interactions_from_data(
            bacterial_sequences = sequence_batch,
            elm_regex = dmi_bundle.elm_regex,
            motif_domains = dmi_bundle.motif_domains,
            human_domains = human_domains,
            motif_sources = dmi_bundle.motif_sources,
        )
        raw_prediction_rows += len(interactions)

        for interaction in interactions:
            pair = (
                interaction.microbial_protein,
                interaction.host_protein,
            )
            unique_predicted_pairs.add(pair)

            if pair not in true_positive_pairs:
                continue

            recovered_pairs.add(pair)
            detail_key = (
                interaction.microbial_protein,
                interaction.host_protein,
                interaction.motif,
                interaction.domain,
                interaction.resource,
            )
            detail_row = detail_lookup.setdefault(
                detail_key,
                {
                    'direction': 'reverse',
                    'bacterial_accession': interaction.microbial_protein,
                    'human_accession': interaction.host_protein,
                    'motif': interaction.motif,
                    'domain': interaction.domain,
                    'resource': interaction.resource,
                    'motif_match_count': 0,
                    'first_start': interaction.start,
                    'first_end': interaction.end,
                },
            )
            detail_row['motif_match_count'] += 1

            if interaction.start < int(detail_row['first_start']):
                detail_row['first_start'] = interaction.start
                detail_row['first_end'] = interaction.end

    true_positives_recovered = len(recovered_pairs)
    unique_predicted_pair_count = len(unique_predicted_pairs)
    false_positive_pairs = len(unique_predicted_pairs - true_positive_pairs)

    summary_row = {
        'approach': 'reverse_dmi',
        'resource_name': FORWARD_DMI_RESOURCE_NAME,
        'benchmark_unique_pairs': len(true_positive_pairs),
        'analyzable_true_pairs': analyzable_true_pairs,
        'resource_rows': '',
        'unique_predicted_pairs': unique_predicted_pair_count,
        'raw_prediction_rows': raw_prediction_rows,
        'true_positives_recovered': true_positives_recovered,
        'false_positive_pairs_within_panel_cross_product': (
            false_positive_pairs
        ),
        'recall_total': _safe_divide(
            true_positives_recovered,
            len(true_positive_pairs),
        ),
        'recall_on_analyzable_pairs': _safe_divide(
            true_positives_recovered,
            analyzable_true_pairs,
        ),
        'precision_unique_pairs': _safe_divide(
            true_positives_recovered,
            unique_predicted_pair_count,
        ),
    }

    detail_frame = _detail_frame_from_records(
        list(detail_lookup.values()),
        columns = [
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

    return summary_row, detail_frame


def build_default_ddi_resource_bundles():
    """
    Build the DDI resource combinations used in the documented benchmark.

    Returns:
        Ordered mapping from benchmark resource-set label to package bundle.
    """

    three_did_bundle = load_default_3did_ddi_resource_bundle()
    domine_hc_bundle = load_default_domine_hc_ddi_resource_bundle()
    domine_all_bundle = load_default_domine_all_ddi_resource_bundle()

    return OrderedDict(
        {
            '3did_current': three_did_bundle,
            'domine_v2_hc': domine_hc_bundle,
            'domine_v2_all': domine_all_bundle,
            '3did_plus_domine_v2_hc': merge_ddi_resource_bundles(
                three_did_bundle,
                domine_hc_bundle,
            ),
            '3did_plus_domine_v2_all': load_default_ddi_resource_bundle(),
        },
    )


def run_ddi_benchmark(
    unique_pairs: pd.DataFrame,
    bacterial_domain_file: str | Path,
    human_domain_file: str | Path,
    resource_name: str,
    resource_bundle,
) -> tuple[dict[str, object], pd.DataFrame]:
    """
    Benchmark one DDI resource set on the PHISTO panel.

    Args:
        unique_pairs: Deduplicated PHISTO benchmark panel.
        bacterial_domain_file: Bacterial protein-to-Pfam table.
        human_domain_file: Human protein-to-Pfam table.
        resource_name: Benchmark label for the DDI resource set.
        resource_bundle: Packaged or merged DDI resource bundle.

    Returns:
        Summary-row dictionary and recovered true-positive detail frame.
    """

    bacterial_domains = read_protein_domain_table(bacterial_domain_file)
    human_domains = read_protein_domain_table(human_domain_file)
    true_positive_pairs = _true_positive_pairs(unique_pairs)
    analyzable_true_pairs = _count_analyzable_pairs(
        true_positive_pairs,
        bacterial_proteins = _domain_proteins(bacterial_domains),
        human_proteins = _domain_proteins(human_domains),
    )

    interactions = predict_domain_domain_interactions_from_data(
        bacterial_domains = bacterial_domains,
        human_domains = human_domains,
        resource_bundle = resource_bundle,
    )
    raw_prediction_rows = len(interactions)
    unique_predicted_pairs: set[tuple[str, str]] = set()
    recovered_pairs: set[tuple[str, str]] = set()
    detail_records: list[dict[str, object]] = []

    for interaction in interactions:
        pair = (
            interaction.bacterial_protein,
            interaction.human_protein,
        )
        unique_predicted_pairs.add(pair)

        if pair not in true_positive_pairs:
            continue

        recovered_pairs.add(pair)
        detail_records.append(
            {
                'resource_name': resource_name,
                'bacterial_accession': interaction.bacterial_protein,
                'human_accession': interaction.human_protein,
                'bacterial_domain': interaction.bacterial_domain,
                'human_domain': interaction.human_domain,
            },
        )

    true_positives_recovered = len(recovered_pairs)
    unique_predicted_pair_count = len(unique_predicted_pairs)
    false_positive_pairs = len(unique_predicted_pairs - true_positive_pairs)

    summary_row = {
        'approach': 'ddi',
        'resource_name': resource_name,
        'benchmark_unique_pairs': len(true_positive_pairs),
        'analyzable_true_pairs': analyzable_true_pairs,
        'resource_rows': len(resource_bundle.pfam_pairs),
        'unique_predicted_pairs': unique_predicted_pair_count,
        'raw_prediction_rows': raw_prediction_rows,
        'true_positives_recovered': true_positives_recovered,
        'false_positive_pairs_within_panel_cross_product': (
            false_positive_pairs
        ),
        'recall_total': _safe_divide(
            true_positives_recovered,
            len(true_positive_pairs),
        ),
        'recall_on_analyzable_pairs': _safe_divide(
            true_positives_recovered,
            analyzable_true_pairs,
        ),
        'precision_unique_pairs': _safe_divide(
            true_positives_recovered,
            unique_predicted_pair_count,
        ),
    }

    detail_frame = _detail_frame_from_records(
        detail_records,
        columns = [
            'resource_name',
            'bacterial_accession',
            'human_accession',
            'bacterial_domain',
            'human_domain',
        ],
    )

    return summary_row, detail_frame


def build_benchmark_overview(
    unique_pairs: pd.DataFrame,
    phisto_row_count: int,
    bacterial_fasta_file: str | Path,
    human_fasta_file: str | Path,
    bacterial_domain_file: str | Path,
    human_domain_file: str | Path,
    bacterial_annotation_file: str | Path | None = None,
    human_annotation_file: str | Path | None = None,
) -> pd.DataFrame:
    """
    Build an overview table for the benchmark panel and structural inputs.

    Args:
        unique_pairs: Deduplicated PHISTO benchmark panel.
        phisto_row_count: Number of rows in the raw PHISTO export.
        bacterial_fasta_file: Bacterial FASTA file.
        human_fasta_file: Human FASTA file.
        bacterial_domain_file: Bacterial protein-to-Pfam table.
        human_domain_file: Human protein-to-Pfam table.
        bacterial_annotation_file: Optional bacterial annotation TSV.
        human_annotation_file: Optional human annotation TSV.

    Returns:
        Overview table matching the documented benchmark summary schema.
    """

    bacterial_sequences = read_fasta_sequences(bacterial_fasta_file)
    human_sequences = read_fasta_sequences(human_fasta_file)
    bacterial_domains = read_protein_domain_table(bacterial_domain_file)
    human_domains = read_protein_domain_table(human_domain_file)

    resolved_bacterial_annotations = ''
    resolved_human_annotations = ''
    unresolved_bacterial_accessions = ''
    unresolved_human_accessions = ''

    if bacterial_annotation_file is not None and Path(
        bacterial_annotation_file,
    ).exists():
        resolved_bacterial_annotations = len(
            pd.read_csv(bacterial_annotation_file, sep = '\t'),
        )
        unresolved_bacterial_accessions = (
            unique_pairs['pathogen_accession'].nunique()
            - resolved_bacterial_annotations
        )

    if human_annotation_file is not None and Path(
        human_annotation_file,
    ).exists():
        resolved_human_annotations = len(
            pd.read_csv(human_annotation_file, sep = '\t'),
        )
        unresolved_human_accessions = (
            unique_pairs['human_accession'].nunique()
            - resolved_human_annotations
        )

    overview_rows = [
        {'metric': 'phisto_rows', 'value': phisto_row_count},
        {'metric': 'phisto_unique_pairs', 'value': len(unique_pairs)},
        {
            'metric': 'unique_bacterial_accessions',
            'value': unique_pairs['pathogen_accession'].nunique(),
        },
        {
            'metric': 'unique_human_accessions',
            'value': unique_pairs['human_accession'].nunique(),
        },
        {
            'metric': 'resolved_bacterial_annotations',
            'value': resolved_bacterial_annotations,
        },
        {
            'metric': 'resolved_human_annotations',
            'value': resolved_human_annotations,
        },
        {
            'metric': 'bacterial_with_pfam',
            'value': len(_domain_proteins(bacterial_domains)),
        },
        {
            'metric': 'human_with_pfam',
            'value': len(_domain_proteins(human_domains)),
        },
        {
            'metric': 'unresolved_bacterial_accessions',
            'value': unresolved_bacterial_accessions,
        },
        {
            'metric': 'unresolved_human_accessions',
            'value': unresolved_human_accessions,
        },
        {
            'metric': 'bacterial_with_sequence',
            'value': len(_sequence_accessions(bacterial_sequences)),
        },
        {
            'metric': 'human_with_sequence',
            'value': len(_sequence_accessions(human_sequences)),
        },
    ]

    return pd.DataFrame(overview_rows)


def run_broad_phisto_benchmark(
    unique_pairs: pd.DataFrame,
    human_fasta_file: str | Path,
    bacterial_domain_file: str | Path,
    bacterial_fasta_file: str | Path,
    human_domain_file: str | Path,
    output_dir: str | Path,
    phisto_row_count: int | None = None,
    bacterial_annotation_file: str | Path | None = None,
    human_annotation_file: str | Path | None = None,
    forward_batch_size: int = DEFAULT_FORWARD_BATCH_SIZE,
    reverse_batch_size: int = DEFAULT_REVERSE_BATCH_SIZE,
) -> dict[str, pd.DataFrame]:
    """
    Run forward DMI, reverse DMI, and DDI on the broad PHISTO panel.

    Args:
        unique_pairs: Deduplicated PHISTO benchmark panel.
        human_fasta_file: Host FASTA file for forward DMI.
        bacterial_domain_file: Bacterial domain table for forward DMI and DDI.
        bacterial_fasta_file: Bacterial FASTA file for reverse DMI.
        human_domain_file: Human domain table for reverse DMI and DDI.
        output_dir: Root directory where benchmark outputs will be written.
        phisto_row_count: Optional raw PHISTO row count for overview output.
        bacterial_annotation_file: Optional bacterial UniProt annotation file.
        human_annotation_file: Optional human UniProt annotation file.
        forward_batch_size: Number of human proteins per forward DMI batch.
        reverse_batch_size: Number of bacterial proteins per reverse DMI batch.

    Returns:
        Mapping from logical output name to written or computed data frame.
    """

    output_root = Path(output_dir)
    results_dir = output_root / '05_results'
    metadata_dir = output_root / '07_metadata'
    results_dir.mkdir(parents = True, exist_ok = True)
    metadata_dir.mkdir(parents = True, exist_ok = True)

    forward_summary, forward_details = run_forward_dmi_benchmark(
        unique_pairs = unique_pairs,
        human_fasta_file = human_fasta_file,
        bacterial_domain_file = bacterial_domain_file,
        batch_size = forward_batch_size,
    )
    reverse_summary, reverse_details = run_reverse_dmi_benchmark(
        unique_pairs = unique_pairs,
        bacterial_fasta_file = bacterial_fasta_file,
        human_domain_file = human_domain_file,
        batch_size = reverse_batch_size,
    )

    ddi_summaries: list[dict[str, object]] = []
    ddi_detail_frames: list[pd.DataFrame] = []
    best_ddi_details: pd.DataFrame | None = None

    for resource_name, resource_bundle in build_default_ddi_resource_bundles().items():
        ddi_summary, ddi_details = run_ddi_benchmark(
            unique_pairs = unique_pairs,
            bacterial_domain_file = bacterial_domain_file,
            human_domain_file = human_domain_file,
            resource_name = resource_name,
            resource_bundle = resource_bundle,
        )
        ddi_summaries.append(ddi_summary)
        ddi_detail_frames.append(ddi_details)

        if resource_name == '3did_plus_domine_v2_all':
            best_ddi_details = ddi_details

    approach_summary = pd.DataFrame(
        [forward_summary, reverse_summary, *ddi_summaries],
    )
    approach_union_summary = pd.DataFrame(
        [
            {
                'set_name': 'forward_dmi',
                'unique_true_positive_pairs': (
                    forward_summary['true_positives_recovered']
                ),
            },
            {
                'set_name': 'reverse_dmi',
                'unique_true_positive_pairs': (
                    reverse_summary['true_positives_recovered']
                ),
            },
            {
                'set_name': 'ddi_3did_plus_domine_v2_all',
                'unique_true_positive_pairs': next(
                    summary_row['true_positives_recovered']
                    for summary_row in ddi_summaries
                    if summary_row['resource_name']
                    == '3did_plus_domine_v2_all'
                ),
            },
            {
                'set_name': 'union_all_three',
                'unique_true_positive_pairs': len(
                    set(
                        zip(
                            forward_details['bacterial_accession'],
                            forward_details['human_accession'],
                        ),
                    )
                    | set(
                        zip(
                            reverse_details['bacterial_accession'],
                            reverse_details['human_accession'],
                        ),
                    )
                    | set(
                        zip(
                            best_ddi_details['bacterial_accession'],
                            best_ddi_details['human_accession'],
                        ),
                    )
                ),
            },
        ],
    )

    if phisto_row_count is None:
        phisto_row_count = len(unique_pairs)

    benchmark_overview = build_benchmark_overview(
        unique_pairs = unique_pairs,
        phisto_row_count = phisto_row_count,
        bacterial_fasta_file = bacterial_fasta_file,
        human_fasta_file = human_fasta_file,
        bacterial_domain_file = bacterial_domain_file,
        human_domain_file = human_domain_file,
        bacterial_annotation_file = bacterial_annotation_file,
        human_annotation_file = human_annotation_file,
    )

    ddi_details = pd.concat(ddi_detail_frames, ignore_index = True)

    _write_frame(unique_pairs, results_dir / 'phisto_unique_pairs.tsv')
    _write_frame(benchmark_overview, results_dir / 'benchmark_overview.tsv')
    _write_frame(approach_summary, results_dir / 'approach_summary.tsv')
    _write_frame(
        approach_union_summary,
        results_dir / 'approach_union_summary.tsv',
    )
    _write_frame(
        forward_details,
        results_dir / 'forward_dmi_true_positive_details.tsv',
    )
    _write_frame(
        reverse_details,
        results_dir / 'reverse_dmi_true_positive_details.tsv',
    )
    _write_frame(ddi_details, results_dir / 'ddi_true_positive_details.tsv')

    run_metadata = pd.DataFrame(
        [
            {
                'setting': 'forward_batch_size',
                'value': forward_batch_size,
            },
            {
                'setting': 'reverse_batch_size',
                'value': reverse_batch_size,
            },
            {
                'setting': 'dmi_resource_bundle',
                'value': FORWARD_DMI_RESOURCE_NAME,
            },
            {
                'setting': 'ddi_resource_bundles',
                'value': ';'.join(build_default_ddi_resource_bundles().keys()),
            },
        ],
    )
    _write_frame(run_metadata, metadata_dir / 'package_run_settings.tsv')

    return {
        'benchmark_overview': benchmark_overview,
        'approach_summary': approach_summary,
        'approach_union_summary': approach_union_summary,
        'forward_dmi_true_positive_details': forward_details,
        'reverse_dmi_true_positive_details': reverse_details,
        'ddi_true_positive_details': ddi_details,
        'phisto_unique_pairs': unique_pairs,
    }


def run_broad_phisto_benchmark_from_bundle(
    benchmark_root: str | Path,
    output_dir: str | Path,
    forward_batch_size: int = DEFAULT_FORWARD_BATCH_SIZE,
    reverse_batch_size: int = DEFAULT_REVERSE_BATCH_SIZE,
) -> dict[str, pd.DataFrame]:
    """
    Run the packaged benchmark workflow using the documented bundle layout.

    Args:
        benchmark_root: Root directory of the PHISTO benchmark bundle.
        output_dir: Root directory where regenerated outputs will be written.
        forward_batch_size: Number of host proteins per forward DMI batch.
        reverse_batch_size: Number of bacterial proteins per reverse DMI batch.

    Returns:
        Mapping from logical output name to benchmark result data frame.
    """

    paths = resolve_benchmark_input_paths(benchmark_root)
    phisto_frame = read_phisto_export(paths.phisto_csv)
    unique_pairs = build_unique_pair_panel(paths.phisto_csv)

    return run_broad_phisto_benchmark(
        unique_pairs = unique_pairs,
        human_fasta_file = paths.human_fasta,
        bacterial_domain_file = paths.bacterial_domains,
        bacterial_fasta_file = paths.bacterial_fasta,
        human_domain_file = paths.human_domains,
        output_dir = output_dir,
        phisto_row_count = len(phisto_frame),
        bacterial_annotation_file = paths.bacterial_annotations,
        human_annotation_file = paths.human_annotations,
        forward_batch_size = forward_batch_size,
        reverse_batch_size = reverse_batch_size,
    )


def load_packaged_prediction_tables(
    benchmark_root: str | Path,
) -> dict[str, pd.DataFrame]:
    """
    Load the already documented benchmark result tables for comparison.

    Args:
        benchmark_root: Root directory of the PHISTO benchmark bundle.

    Returns:
        Mapping from table name to data frame.
    """

    root = Path(benchmark_root)
    results_dir = root / '05_results'

    return {
        'benchmark_overview': pd.read_csv(
            results_dir / 'benchmark_overview.tsv',
            sep = '\t',
        ),
        'approach_summary': pd.read_csv(
            results_dir / 'approach_summary.tsv',
            sep = '\t',
        ),
        'approach_union_summary': pd.read_csv(
            results_dir / 'approach_union_summary.tsv',
            sep = '\t',
        ),
        'forward_dmi_true_positive_details': pd.read_csv(
            results_dir / 'forward_dmi_true_positive_details.tsv',
            sep = '\t',
        ),
        'reverse_dmi_true_positive_details': pd.read_csv(
            results_dir / 'reverse_dmi_true_positive_details.tsv',
            sep = '\t',
        ),
        'ddi_true_positive_details': pd.read_csv(
            results_dir / 'ddi_true_positive_details.tsv',
            sep = '\t',
        ),
        'phisto_unique_pairs': pd.read_csv(
            results_dir / 'phisto_unique_pairs.tsv',
            sep = '\t',
        ),
    }


def preview_prediction_shapes(
    benchmark_root: str | Path,
) -> dict[str, pd.DataFrame]:
    """
    Run the three package predictors once and return their tabular schemas.

    This helper is intentionally lightweight. It reuses the prepared benchmark
    inputs and returns small head tables that are useful inside the notebook
    when introducing the forward DMI, reverse DMI, and DDI output schemas.

    Args:
        benchmark_root: Root directory of the PHISTO benchmark bundle.

    Returns:
        Mapping from predictor name to a five-row preview data frame.
    """

    paths = resolve_benchmark_input_paths(benchmark_root)
    dmi_bundle = load_default_dmi_resource_bundle()
    ddi_bundle = load_default_ddi_resource_bundle()
    human_sequences = read_fasta_sequences(paths.human_fasta)
    bacterial_sequences = read_fasta_sequences(paths.bacterial_fasta)
    bacterial_domains = read_protein_domain_table(paths.bacterial_domains)
    human_domains = read_protein_domain_table(paths.human_domains)
    bacterial_subset = {
        extract_uniprot_id(header)
        for header in list(bacterial_sequences)[:25]
    }
    human_subset = {
        extract_uniprot_id(header)
        for header in list(human_sequences)[:25]
    }

    forward_interactions = predict_domain_motif_interactions_from_data(
        human_sequences = dict(list(human_sequences.items())[:25]),
        elm_regex = dmi_bundle.elm_regex,
        motif_domains = dmi_bundle.motif_domains,
        bacterial_domains = bacterial_domains,
        motif_sources = dmi_bundle.motif_sources,
    )
    reverse_interactions = predict_reverse_domain_motif_interactions_from_data(
        bacterial_sequences = dict(list(bacterial_sequences.items())[:25]),
        elm_regex = dmi_bundle.elm_regex,
        motif_domains = dmi_bundle.motif_domains,
        human_domains = _subset_domain_mapping_by_proteins(
            human_domains,
            human_subset,
        ),
        motif_sources = dmi_bundle.motif_sources,
    )
    ddi_interactions = predict_domain_domain_interactions_from_data(
        bacterial_domains = _subset_domain_mapping_by_proteins(
            bacterial_domains,
            bacterial_subset,
        ),
        human_domains = _subset_domain_mapping_by_proteins(
            human_domains,
            human_subset,
        ),
        resource_bundle = ddi_bundle,
    )

    return {
        'forward_dmi': interactions_to_dataframe(forward_interactions).head(),
        'reverse_dmi': bidirectional_interactions_to_dataframe(
            reverse_interactions,
        ).head(),
        'ddi': ddi_interactions_to_dataframe(ddi_interactions).head(),
    }
