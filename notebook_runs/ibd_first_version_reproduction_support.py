from __future__ import annotations

from dataclasses import dataclass
import importlib.util
from io import StringIO
from pathlib import Path
import re
import shlex
import subprocess
import sys
from typing import Iterable
from typing import Optional
from typing import Union

import numpy as np
import pandas as pd
import requests


UNIPROT_ACCESSION_PATTERN = re.compile(
    (
        r'^(?:'
        r'[OPQ][0-9][A-Z0-9]{3}[0-9]'
        r'|'
        r'[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
        r')'
        r'(?:-\d+)?$'
    ),
)
UNIPROT_ACCESSION_CAPTURE_PATTERN = re.compile(
    (
        r'\b(?:'
        r'[OPQ][0-9][A-Z0-9]{3}[0-9]'
        r'|'
        r'[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}'
        r')'
        r'(?:-\d+)?\b'
    ),
)
PFAM_CAPTURE_PATTERN = re.compile(r'PF\d{5}')
EC_CAPTURE_PATTERN = re.compile(r'\b\d+\.\d+\.\d+\.\d+\b')
KEGG_KO_CAPTURE_PATTERN = re.compile(r'\bK\d{5}\b')
KEGG_GENE_CAPTURE_PATTERN = re.compile(r'\b[a-z]{3,4}:[A-Za-z0-9_.-]+\b')
NCBI_PROTEIN_ACCESSION_CAPTURE_PATTERN = re.compile(
    r'\b(?:[A-Z]{2}_[0-9]+\.[0-9]+|[A-Z]{3,5}[0-9]+\.[0-9]+)\b',
)
OMA_GENE_NAME_CAPTURE_PATTERN = re.compile(r'\bGN=([^ ]+)')
TIEDIE_HUMAN_COLUMN_CANDIDATES = [
    '# Human Protein',
    'Human Protein',
    'human_protein',
]
TIEDIE_BACTERIAL_COLUMN_CANDIDATES = [
    'Bacterial protein',
    'Bacteria Protein',
    'bacterial_protein',
]
TIEDIE_SIGN_COLUMN_CANDIDATES = [
    'sign',
    'Sign',
]
DEFAULT_FASTA_BATCH_SIZE = 100
DEFAULT_ANNOTATION_BATCH_SIZE = 250
DEFAULT_UNIPROT_FIELDS = [
    'accession',
    'id',
    'gene_names',
    'xref_pfam',
]
DEFAULT_MMC3_UNIPROT_JSON_BATCH_SIZE = 100
DEFAULT_MMC3_KEGG_BATCH_SIZE = 100
DEFAULT_LOCATION_FILTERS = [
    'plasma_membrane_transmembrane',
    'plasma_membrane_peripheral',
    'secreted',
]
MMC3_DEFAULT_SET_NAMES = [
    'protein abundances',
    'KO-level protein abundances',
    'EC-level protein abundances',
    'KO transcription',
    'EC transcription',
    'metagenomic KO profiles',
    'Metagenomic EC profiles',
]
MMC3_KO_SET_NAMES = [
    'KO-level protein abundances',
    'KO transcription',
    'metagenomic KO profiles',
]
MMC3_EC_SET_NAMES = [
    'EC-level protein abundances',
    'EC transcription',
    'Metagenomic EC profiles',
]
MMC3_PROTEIN_SET_NAMES = [
    'protein abundances',
]
MMC3_GENERIC_DESCRIPTION_TEXT = {
    'no name',
    'hypothetical protein',
    'putative enzyme',
    'putative protein',
    'unknown protein',
}
PREFERRED_PIPELINE_RELATIVE_PATH = Path('worktrees') / 'MicrobioLink-2.1-beta'
PathLike = Union[str, Path]
ColumnName = Union[str, int]
SheetName = Union[int, str]


@dataclass
class Mmc3OrthogroupFilterResult:
    """Hold the tables produced while filtering OMA members with mmc3."""

    microbial_proteins: pd.DataFrame
    member_annotations: pd.DataFrame
    mmc3_feature_tokens: pd.DataFrame
    match_evidence: pd.DataFrame
    matched_groups: pd.DataFrame


def unique_preserve_order(values: Iterable[str]) -> list[str]:
    """Return unique non-empty values while preserving input order."""

    seen: set[str] = set()
    ordered: list[str] = []

    for value in values:
        if value is None:
            continue

        normalized = str(value).strip()
        if not normalized or normalized in seen:
            continue

        seen.add(normalized)
        ordered.append(normalized)

    return ordered


def normalize_text_for_matching(value: object) -> str:
    """Normalize free text for exact, punctuation-insensitive matching."""

    text = str(value).strip().lower()
    text = re.sub(r'[^a-z0-9]+', ' ', text)
    return re.sub(r'\s+', ' ', text).strip()


def split_semicolon_values(value: object) -> list[str]:
    """Split a UniProt-style semicolon-delimited field into tokens."""

    tokens = [
        token.strip()
        for token in str(value).split(';')
    ]
    return [
        token
        for token in tokens
        if token and token.lower() != 'nan'
    ]


def partition_accessions_by_uniprot_pattern(
    accessions: Iterable[str],
) -> tuple[list[str], list[str]]:
    """Split accessions into UniProt-like and non-UniProt-like groups."""

    ordered_accessions = unique_preserve_order(
        str(accession).strip()
        for accession in accessions
    )
    uniprot_like = [
        accession
        for accession in ordered_accessions
        if UNIPROT_ACCESSION_PATTERN.match(accession)
    ]
    other = [
        accession
        for accession in ordered_accessions
        if accession not in set(uniprot_like)
    ]
    return uniprot_like, other


def _extract_generic_fasta_accession(header: str) -> str:
    """Extract the leading accession token from a generic FASTA header."""

    stripped_header = header.strip()

    if '|' in stripped_header:
        fields = stripped_header.split('|')

        if len(fields) >= 2 and fields[1].strip():
            return fields[1].strip()

    return stripped_header.split(None, 1)[0].strip()


def write_selected_fasta_from_reference_fastas(
    selected_accessions: Iterable[str],
    source_fasta_paths: Iterable[PathLike],
    output_path: PathLike,
    strict: bool = True,
    normalize_headers: bool = False,
) -> tuple[Path, list[str]]:
    """Write a FASTA subset by pulling matching accessions from source FASTAs."""

    ordered_accessions = unique_preserve_order(selected_accessions)
    accession_set = set(ordered_accessions)
    output_file = Path(output_path)
    output_file.parent.mkdir(parents = True, exist_ok = True)
    found_records: dict[str, tuple[str, str]] = {}

    for source_fasta_path in source_fasta_paths:
        source_path = Path(source_fasta_path)

        if not source_path.exists():
            continue

        current_header: Optional[str] = None
        current_fragments: list[str] = []

        with open(source_path, encoding = 'utf-8') as fasta_handle:
            for raw_line in fasta_handle:
                line = raw_line.strip()

                if not line:
                    continue

                if line.startswith('>'):
                    if current_header is not None:
                        accession = _extract_generic_fasta_accession(current_header)

                        if accession in accession_set and accession not in found_records:
                            found_records[accession] = (
                                current_header,
                                ''.join(current_fragments),
                            )

                    current_header = line[1:]
                    current_fragments = []
                    continue

                current_fragments.append(line)

            if current_header is not None:
                accession = _extract_generic_fasta_accession(current_header)

                if accession in accession_set and accession not in found_records:
                    found_records[accession] = (
                        current_header,
                        ''.join(current_fragments),
                    )

    missing_accessions = [
        accession
        for accession in ordered_accessions
        if accession not in found_records
    ]

    if missing_accessions and strict:
        raise FileNotFoundError(
            'Could not find all requested accessions in the reference FASTA '
            'sources. Missing accessions:\n'
            + '\n'.join(missing_accessions),
        )

    with open(output_file, 'w', encoding = 'utf-8') as fasta_handle:
        for accession in ordered_accessions:
            if accession not in found_records:
                continue

            header, sequence = found_records[accession]
            normalized_header = (
                f'microbe|{accession}|{accession}'
                if normalize_headers
                else header
            )
            fasta_handle.write(f'>{normalized_header}\n')
            fasta_handle.write(f'{sequence}\n')

    return output_file, missing_accessions


def find_repo_root(start: Optional[Path] = None) -> Path:
    """Find the top-level repository directory from any child path."""

    current = (start or Path.cwd()).resolve()

    for candidate in [current, *current.parents]:
        if (candidate / 'AGENTS.md').exists() and (candidate / 'pyproject.toml').exists():
            return candidate

    raise FileNotFoundError(
        'Could not locate the MicrobioLink repository root from the current '
        'working directory.',
    )


def choose_pipeline_root(
    repo_root: Path,
    preferred_relative_path: Path = PREFERRED_PIPELINE_RELATIVE_PATH,
) -> Path:
    """Pick the MicrobioLink checkout that exposes DDI in this workspace."""

    candidates = [
        repo_root / preferred_relative_path,
        repo_root,
    ]

    for candidate in candidates:
        if (
            candidate.exists()
            and candidate.joinpath('microbiolink').is_dir()
            and candidate.joinpath('microbiolink_api').is_dir()
            and candidate.joinpath('microbiolink_api', 'ddi.py').exists()
        ):
            return candidate

    raise FileNotFoundError(
        'Could not find a MicrobioLink checkout with packaged DDI support. '
        'Expected either the repository root or '
        f'{preferred_relative_path.as_posix()}.',
    )


def choose_monte_carlo_root(
    repo_root: Path,
    preferred_relative_path: Path = PREFERRED_PIPELINE_RELATIVE_PATH,
) -> Path:
    """Pick the checkout that exposes the Monte Carlo motif filter."""

    candidates = [
        repo_root,
        repo_root / preferred_relative_path,
    ]

    for candidate in candidates:
        if candidate.joinpath(
            'microbiolink',
            'motif_monte_carlo_filter.py',
        ).exists():
            return candidate

    raise FileNotFoundError(
        'Could not find a MicrobioLink checkout with '
        '`microbiolink/motif_monte_carlo_filter.py`.',
    )


def activate_pipeline_root(pipeline_root: Path) -> Path:
    """Place the chosen MicrobioLink checkout at the front of sys.path."""

    normalized = str(pipeline_root.resolve())

    if normalized in sys.path:
        sys.path.remove(normalized)

    sys.path.insert(0, normalized)
    return pipeline_root


def load_module_from_path(
    module_name: str,
    module_path: PathLike,
):
    """Load a Python module directly from a filesystem path."""

    resolved_path = Path(module_path).resolve()
    spec = importlib.util.spec_from_file_location(
        module_name,
        resolved_path,
    )

    if spec is None or spec.loader is None:
        raise ImportError(
            f'Could not create an import spec for {resolved_path}.',
        )

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_run_monte_carlo_filter(
    monte_carlo_root: PathLike,
):
    """Load `run_monte_carlo_filter` from the checkout that provides it."""

    module = load_module_from_path(
        'microbiolink_motif_monte_carlo_filter_external',
        Path(monte_carlo_root).joinpath(
            'microbiolink',
            'motif_monte_carlo_filter.py',
        ),
    )
    return module.run_monte_carlo_filter


def inspect_aiupred_runtime() -> dict[str, object]:
    """Report whether AIUPred can run in the current Python environment."""

    runtime = {
        'iupred_available': False,
        'torch_available': False,
        'aiupred_ready': False,
        'message': '',
    }

    if importlib.util.find_spec('iupred') is None:
        runtime['message'] = (
            'The `iupred` package is not installed in this environment.'
        )
        return runtime

    runtime['iupred_available'] = True

    if importlib.util.find_spec('torch') is None:
        runtime['message'] = (
            'AIUPred is installed, but PyTorch is missing: '
            'No module named `torch`.'
        )
        return runtime

    runtime['torch_available'] = True
    runtime['aiupred_ready'] = True
    runtime['message'] = 'AIUPred and PyTorch are available.'
    return runtime


def ensure_directory(path: PathLike) -> Path:
    """Create a directory if it does not exist and return its Path."""

    directory = Path(path)
    directory.mkdir(parents = True, exist_ok = True)
    return directory


def read_table(
    input_path: PathLike,
    sheet_name: SheetName = 0,
) -> pd.DataFrame:
    """Read a CSV, TSV, TXT, or Excel table with light auto-detection."""

    path = Path(input_path)
    suffix = path.suffix.lower()

    if suffix in {'.xlsx', '.xls', '.xlsm'}:
        return pd.read_excel(path, sheet_name = sheet_name)

    if suffix == '.tsv':
        return pd.read_csv(path, sep = '\t')

    return pd.read_csv(path, sep = None, engine = 'python')


def resolve_column_name(
    frame: pd.DataFrame,
    column: ColumnName,
) -> str:
    """Resolve a column identifier supplied as a name or one-based index."""

    if isinstance(column, int):
        column_index = column - 1
        if column_index < 0 or column_index >= len(frame.columns):
            raise ValueError(
                f'Column index {column} is out of range for {list(frame.columns)}.',
            )
        return str(frame.columns[column_index])

    if column not in frame.columns:
        raise ValueError(
            f'Column `{column}` not found in {list(frame.columns)}.',
        )

    return str(column)


def normalize_species_label(value: object) -> str:
    """Normalize a species label to a filesystem-friendly token."""

    text = str(value).strip()
    text = re.sub(r'[^A-Za-z0-9]+', '_', text)
    text = re.sub(r'_+', '_', text)
    return text.strip('_')


def parse_accession_list(value: object) -> list[str]:
    """Extract UniProt accessions from a scalar or list-like cell value."""

    if isinstance(value, (list, tuple, set)):
        tokens = [str(token).strip() for token in value]
        return [
            token
            for token in unique_preserve_order(tokens)
            if UNIPROT_ACCESSION_PATTERN.match(token)
        ]

    if value is None or pd.isna(value):
        return []

    text = str(value)
    matches = UNIPROT_ACCESSION_CAPTURE_PATTERN.findall(text)

    if matches:
        return unique_preserve_order(matches)

    tokens = re.split(r'[\s,;|/]+', text)
    return [
        token
        for token in unique_preserve_order(tokens)
        if UNIPROT_ACCESSION_PATTERN.match(token)
    ]


def load_accession_table(
    input_path: PathLike,
    accession_column: ColumnName,
    sheet_name: SheetName = 0,
    gene_column: Optional[ColumnName] = None,
) -> pd.DataFrame:
    """Standardize a table containing UniProt accessions."""

    frame = read_table(input_path, sheet_name = sheet_name)
    return load_accession_frame(
        frame,
        accession_column = accession_column,
        gene_column = gene_column,
    )


def load_accession_frame(
    frame: pd.DataFrame,
    accession_column: ColumnName,
    gene_column: Optional[ColumnName] = None,
) -> pd.DataFrame:
    """Standardize a data frame containing UniProt accessions."""

    resolved_accession = resolve_column_name(frame, accession_column)
    selected_columns = [resolved_accession]

    if gene_column is not None:
        selected_columns.append(resolve_column_name(frame, gene_column))

    result = frame[selected_columns].copy()
    result['uniprot'] = result[resolved_accession].map(parse_accession_list)
    result = result.explode('uniprot')
    result = result.dropna(subset = ['uniprot'])
    result['uniprot'] = result['uniprot'].astype(str).str.strip()
    result = result[result['uniprot'].str.match(UNIPROT_ACCESSION_PATTERN)]

    rename_map: dict[str, str] = {}
    if gene_column is not None:
        rename_map[selected_columns[1]] = 'genesymbol'

    result = result.rename(columns = rename_map)
    kept_columns = ['uniprot']
    if 'genesymbol' in result.columns:
        kept_columns.append('genesymbol')

    return (
        result[kept_columns]
        .drop_duplicates()
        .sort_values(kept_columns)
        .reset_index(drop = True)
    )


def filter_rows_by_text_keywords(
    frame: pd.DataFrame,
    filter_columns: list[ColumnName],
    include_keywords: Optional[list[str]] = None,
    exclude_keywords: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Filter a frame by keyword presence across one or more text columns."""

    resolved_columns = [
        resolve_column_name(frame, column_name)
        for column_name in filter_columns
    ]

    row_text = pd.Series('', index = frame.index, dtype = str)
    for column_name in resolved_columns:
        row_text = row_text + ' ' + frame[column_name].fillna('').astype(str)

    row_text = row_text.str.lower()
    mask = pd.Series(True, index = frame.index)

    normalized_include = [
        keyword.strip().lower()
        for keyword in include_keywords or []
        if str(keyword).strip()
    ]
    if normalized_include:
        include_pattern = '|'.join(
            re.escape(keyword)
            for keyword in normalized_include
        )
        mask &= row_text.str.contains(include_pattern, regex = True)

    normalized_exclude = [
        keyword.strip().lower()
        for keyword in exclude_keywords or []
        if str(keyword).strip()
    ]
    if normalized_exclude:
        exclude_pattern = '|'.join(
            re.escape(keyword)
            for keyword in normalized_exclude
        )
        mask &= ~row_text.str.contains(exclude_pattern, regex = True)

    return frame.loc[mask].copy()


def load_hpa_candidate_accession_table(
    input_path: PathLike,
    accession_column: ColumnName = 'Uniprot',
    gene_column: ColumnName = 'Gene',
    sheet_name: SheetName = 0,
    filter_columns: Optional[list[ColumnName]] = None,
    include_keywords: Optional[list[str]] = None,
    exclude_keywords: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Load a wide Human Protein Atlas export and normalize accessions."""

    frame = read_table(input_path, sheet_name = sheet_name)

    if filter_columns:
        frame = filter_rows_by_text_keywords(
            frame,
            filter_columns = filter_columns,
            include_keywords = include_keywords,
            exclude_keywords = exclude_keywords,
        )

    return load_accession_frame(
        frame,
        accession_column = accession_column,
        gene_column = gene_column,
    )


def prepare_microbial_protein_table(
    supplementary_table_path: PathLike,
    species_column: Optional[ColumnName],
    oma_group_column: ColumnName,
    microbial_protein_column: Optional[ColumnName] = None,
    sheet_name: SheetName = 0,
    membership_table_path: Optional[PathLike] = None,
    membership_group_column: Optional[ColumnName] = None,
    membership_protein_column: Optional[ColumnName] = None,
    membership_species_column: Optional[ColumnName] = None,
) -> pd.DataFrame:
    """Build the three-column microbial protein table for the notebook."""

    supplementary_table = read_table(
        supplementary_table_path,
        sheet_name = sheet_name,
    )
    resolved_group = resolve_column_name(supplementary_table, oma_group_column)
    resolved_species: Optional[str] = None

    if species_column is not None:
        resolved_species = resolve_column_name(
            supplementary_table,
            species_column,
        )

    if microbial_protein_column is not None:
        if resolved_species is None:
            raise ValueError(
                'Provide `species_column` when the supplementary table already '
                'contains individual microbial protein accessions.',
            )

        resolved_protein = resolve_column_name(
            supplementary_table,
            microbial_protein_column,
        )
        protein_table = supplementary_table[
            [resolved_species, resolved_group, resolved_protein]
        ].copy()
        protein_table['microbial_protein'] = protein_table[
            resolved_protein
        ].map(parse_accession_list)
        species_source_column = resolved_species
    elif membership_table_path is not None:
        membership_table = read_table(membership_table_path)
        if membership_group_column is None or membership_protein_column is None:
            raise ValueError(
                'Provide both `membership_group_column` and '
                '`membership_protein_column` when using `membership_table_path`.',
            )

        resolved_membership_group = resolve_column_name(
            membership_table,
            membership_group_column,
        )
        resolved_membership_protein = resolve_column_name(
            membership_table,
            membership_protein_column,
        )
        membership_columns = [
            resolved_membership_group,
            resolved_membership_protein,
        ]

        resolved_membership_species: Optional[str] = None
        if membership_species_column is not None:
            resolved_membership_species = resolve_column_name(
                membership_table,
                membership_species_column,
            )
            membership_columns.append(resolved_membership_species)

        membership_subset = membership_table[membership_columns].copy()
        membership_subset['microbial_protein'] = membership_subset[
            resolved_membership_protein
        ].map(parse_accession_list)
        membership_subset = membership_subset.explode('microbial_protein')
        membership_subset = membership_subset.dropna(subset = ['microbial_protein'])
        membership_subset['microbial_protein'] = (
            membership_subset['microbial_protein'].astype(str).str.strip()
        )
        membership_subset = membership_subset[
            membership_subset['microbial_protein'].str.match(
                UNIPROT_ACCESSION_PATTERN,
            )
        ]

        if resolved_species is not None:
            supplementary_subset = supplementary_table[
                [resolved_species, resolved_group]
            ].drop_duplicates()
        else:
            supplementary_subset = supplementary_table[
                [resolved_group]
            ].drop_duplicates()

        protein_table = supplementary_subset.merge(
            membership_subset,
            left_on = resolved_group,
            right_on = resolved_membership_group,
            how = 'inner',
        )

        if resolved_species is None:
            if resolved_membership_species is None:
                raise ValueError(
                    'Provide `membership_species_column` when the '
                    'supplementary orthogroup table does not contain species.',
                )
            species_source_column = resolved_membership_species
        else:
            species_source_column = resolved_species

        if resolved_species is not None and resolved_membership_species is not None:
            protein_table[resolved_species] = protein_table[resolved_species].fillna(
                protein_table[resolved_membership_species],
            )
    else:
        raise ValueError(
            'Provide either `microbial_protein_column` or a membership table.',
        )

    protein_table = protein_table.explode('microbial_protein')
    protein_table = protein_table.dropna(subset = ['microbial_protein'])
    protein_table['microbial_protein'] = (
        protein_table['microbial_protein'].astype(str).str.strip()
    )
    protein_table = protein_table[
        protein_table['microbial_protein'].str.match(UNIPROT_ACCESSION_PATTERN)
    ].copy()
    protein_table['species'] = protein_table[species_source_column].map(
        normalize_species_label,
    )
    protein_table['oma_group'] = protein_table[resolved_group].astype(str).str.strip()

    return (
        protein_table[['species', 'oma_group', 'microbial_protein']]
        .drop_duplicates()
        .sort_values(['species', 'oma_group', 'microbial_protein'])
        .reset_index(drop = True)
    )


def _extract_mmc3_feature_description(feature: str) -> str:
    """Extract the descriptive part of an mmc3 feature label."""

    if ':' not in feature:
        return ''

    description = feature.split(':', 1)[1]
    description = description.split('[EC:', 1)[0]
    return description.strip()


def parse_mmc3_feature_table(
    input_path: PathLike,
    set_column: ColumnName = 'Set',
    feature_column: ColumnName = 'Feature',
    effect_column: Optional[ColumnName] = 'Dysbiosis Coefficient (CD)',
    allowed_sets: Optional[list[str]] = None,
    cd_increased_only: bool = True,
    sheet_name: SheetName = 0,
) -> pd.DataFrame:
    """Convert mmc3 features into an exact-match token table."""

    frame = read_table(input_path, sheet_name = sheet_name)
    resolved_set = resolve_column_name(frame, set_column)
    resolved_feature = resolve_column_name(frame, feature_column)

    if effect_column is not None:
        resolved_effect = resolve_column_name(frame, effect_column)
        frame[resolved_effect] = pd.to_numeric(
            frame[resolved_effect],
            errors = 'coerce',
        )

        if cd_increased_only:
            frame = frame[frame[resolved_effect] > 0].copy()
    else:
        resolved_effect = None

    selected_sets = allowed_sets or MMC3_DEFAULT_SET_NAMES
    frame = frame[frame[resolved_set].astype(str).isin(selected_sets)].copy()
    rows: list[dict[str, object]] = []

    for _, row in frame.iterrows():
        set_name = str(row[resolved_set]).strip()
        feature = str(row[resolved_feature]).strip()
        description = _extract_mmc3_feature_description(feature)
        normalized_description = normalize_text_for_matching(description)
        emitted: set[tuple[str, str]] = set()

        def emit(
            match_type: str,
            match_token: str,
            token_source: str,
        ) -> None:
            normalized_token = str(match_token).strip()

            if not normalized_token:
                return

            key = (match_type, normalized_token)

            if key in emitted:
                return

            emitted.add(key)
            rows.append(
                {
                    'set_name': set_name,
                    'feature': feature,
                    'match_type': match_type,
                    'match_token': normalized_token,
                    'token_source': token_source,
                    'feature_description': description,
                    'normalized_feature_description': normalized_description,
                    'dysbiosis_coefficient_cd': (
                        row[resolved_effect]
                        if resolved_effect is not None
                        else np.nan
                    ),
                },
            )

        if set_name in MMC3_PROTEIN_SET_NAMES:
            for accession in unique_preserve_order(
                NCBI_PROTEIN_ACCESSION_CAPTURE_PATTERN.findall(feature),
            ):
                emit(
                    match_type = 'protein_accession',
                    match_token = accession,
                    token_source = 'protein_accession',
                )

        if set_name in MMC3_KO_SET_NAMES:
            for ko_id in unique_preserve_order(
                KEGG_KO_CAPTURE_PATTERN.findall(feature),
            ):
                emit(
                    match_type = 'ko_id',
                    match_token = ko_id,
                    token_source = 'ko_id',
                )

            if (
                normalized_description
                and normalized_description not in MMC3_GENERIC_DESCRIPTION_TEXT
            ):
                emit(
                    match_type = 'protein_name',
                    match_token = normalized_description,
                    token_source = 'normalized_description',
                )

        if set_name in MMC3_KO_SET_NAMES or set_name in MMC3_EC_SET_NAMES:
            for ec_number in unique_preserve_order(
                EC_CAPTURE_PATTERN.findall(feature),
            ):
                emit(
                    match_type = 'ec_number',
                    match_token = ec_number,
                    token_source = 'ec_number',
                )

    if not rows:
        return pd.DataFrame(
            columns = [
                'set_name',
                'feature',
                'match_type',
                'match_token',
                'token_source',
                'feature_description',
                'normalized_feature_description',
                'dysbiosis_coefficient_cd',
            ],
        )

    return (
        pd.DataFrame(rows)
        .drop_duplicates()
        .sort_values(['set_name', 'feature', 'match_type', 'match_token'])
        .reset_index(drop = True)
    )


def _parse_oma_member_label(
    microbial_protein: str,
    oma_member_label: str,
) -> dict[str, str]:
    """Parse a local OMA member label into names that can be matched."""

    label = str(oma_member_label).strip()
    accession = str(microbial_protein).strip()
    gene_match = OMA_GENE_NAME_CAPTURE_PATTERN.search(label)
    gene_name = gene_match.group(1).strip() if gene_match else ''
    protein_name = ''
    organism_name = ''

    if ' OS=' in label:
        prefix, suffix = label.split(' OS=', 1)
        organism_name = suffix.split(' OX=', 1)[0].strip()

        if '|' in prefix:
            trailing = prefix.split('|', 2)[-1]

            if ' ' in trailing:
                protein_name = trailing.split(' ', 1)[1].strip()
        else:
            protein_name = prefix.strip()
    else:
        parts = label.split(' ', 1)

        if len(parts) == 2:
            protein_name = parts[1].strip()

        if protein_name.endswith(']') and '[' in protein_name:
            protein_name, bracketed_species = protein_name.rsplit('[', 1)
            protein_name = protein_name.strip()
            organism_name = bracketed_species.rstrip(']').strip()

    return {
        'local_accession': accession,
        'local_gene_name': gene_name,
        'local_protein_name': protein_name,
        'local_protein_name_normalized': normalize_text_for_matching(protein_name),
        'local_organism_name': organism_name,
    }


def download_uniprot_json_annotation_table(
    uniprot_ids: Iterable[str],
    output_path: PathLike,
    batch_size: int = DEFAULT_MMC3_UNIPROT_JSON_BATCH_SIZE,
    timeout: int = 120,
) -> pd.DataFrame:
    """Download UniProt JSON records and flatten the fields needed for mmc3."""

    ordered_ids = unique_preserve_order(
        accession
        for accession in uniprot_ids
        if UNIPROT_ACCESSION_PATTERN.match(str(accession).strip())
    )
    output_file = Path(output_path)
    output_file.parent.mkdir(parents = True, exist_ok = True)
    empty_columns = [
        'microbial_protein',
        'uniprot_protein_name',
        'uniprot_gene_names',
        'uniprot_organism_name',
        'ec_numbers',
        'refseq_accessions',
        'embl_protein_accessions',
        'kegg_gene_ids',
    ]

    if not ordered_ids:
        empty_result = pd.DataFrame(columns = empty_columns)
        empty_result.to_csv(output_file, sep = '\t', index = False)
        return empty_result

    rows: list[dict[str, str]] = []

    for start in range(0, len(ordered_ids), batch_size):
        batch = ordered_ids[start : start + batch_size]
        response = requests.get(
            'https://rest.uniprot.org/uniprotkb/stream',
            params = {
                'format': 'json',
                'query': _build_uniprot_query(batch),
            },
            timeout = timeout,
        )
        response.raise_for_status()
        payload = response.json()

        for record in payload.get('results', []):
            accession = str(record.get('primaryAccession', '')).strip()
            protein_description = record.get('proteinDescription', {})
            recommended_name = protein_description.get('recommendedName', {})
            full_name = recommended_name.get('fullName', {}).get('value', '')
            ec_numbers = [
                str(item.get('value', '')).strip()
                for item in recommended_name.get('ecNumbers', [])
                if str(item.get('value', '')).strip()
            ]
            gene_names: list[str] = []

            for gene_entry in record.get('genes', []):
                for key in ['geneName', 'orderedLocusNames', 'synonyms']:
                    gene_value = gene_entry.get(key)

                    if isinstance(gene_value, dict):
                        candidate = str(gene_value.get('value', '')).strip()
                        if candidate:
                            gene_names.append(candidate)
                    elif isinstance(gene_value, list):
                        for item in gene_value:
                            candidate = str(item.get('value', '')).strip()
                            if candidate:
                                gene_names.append(candidate)

            refseq_ids: list[str] = []
            embl_protein_ids: list[str] = []
            kegg_gene_ids: list[str] = []

            for cross_reference in record.get('uniProtKBCrossReferences', []):
                database_name = cross_reference.get('database')

                if database_name == 'RefSeq':
                    refseq_ids.append(str(cross_reference.get('id', '')).strip())

                if database_name == 'EMBL':
                    for property_entry in cross_reference.get('properties', []):
                        if property_entry.get('key') == 'ProteinId':
                            candidate = str(property_entry.get('value', '')).strip()
                            if candidate:
                                embl_protein_ids.append(candidate)

                if database_name == 'KEGG':
                    candidate = str(cross_reference.get('id', '')).strip()
                    if candidate:
                        kegg_gene_ids.append(candidate)

            rows.append(
                {
                    'microbial_protein': accession,
                    'uniprot_protein_name': full_name,
                    'uniprot_gene_names': ';'.join(
                        unique_preserve_order(gene_names),
                    ),
                    'uniprot_organism_name': str(
                        record.get('organism', {}).get('scientificName', ''),
                    ).strip(),
                    'ec_numbers': ';'.join(unique_preserve_order(ec_numbers)),
                    'refseq_accessions': ';'.join(
                        unique_preserve_order(refseq_ids),
                    ),
                    'embl_protein_accessions': ';'.join(
                        unique_preserve_order(embl_protein_ids),
                    ),
                    'kegg_gene_ids': ';'.join(
                        unique_preserve_order(kegg_gene_ids),
                    ),
                },
            )

    if rows:
        result = (
            pd.DataFrame(rows)
            .drop_duplicates()
            .sort_values('microbial_protein')
            .reset_index(drop = True)
        )
    else:
        result = pd.DataFrame(columns = empty_columns)

    result.to_csv(output_file, sep = '\t', index = False)
    return result


def download_kegg_gene_to_ko_table(
    kegg_gene_ids: Iterable[str],
    output_path: PathLike,
    batch_size: int = DEFAULT_MMC3_KEGG_BATCH_SIZE,
    timeout: int = 120,
) -> pd.DataFrame:
    """Map KEGG gene identifiers to KO identifiers in batches."""

    ordered_gene_ids = unique_preserve_order(
        gene_id
        for gene_id in kegg_gene_ids
        if KEGG_GENE_CAPTURE_PATTERN.match(str(gene_id).strip())
    )
    output_file = Path(output_path)
    output_file.parent.mkdir(parents = True, exist_ok = True)
    rows: list[dict[str, str]] = []

    for start in range(0, len(ordered_gene_ids), batch_size):
        batch = ordered_gene_ids[start : start + batch_size]
        batch_query = '+'.join(batch)
        response = requests.get(
            f'https://rest.kegg.jp/link/ko/{batch_query}',
            timeout = timeout,
        )
        response.raise_for_status()

        for line in response.text.splitlines():
            if not line.strip():
                continue

            left, right = line.split('\t', 1)
            rows.append(
                {
                    'kegg_gene_id': left.strip(),
                    'ko_id': right.replace('ko:', '').strip(),
                },
            )

    if rows:
        result = (
            pd.DataFrame(rows)
            .drop_duplicates()
            .sort_values(['kegg_gene_id', 'ko_id'])
            .reset_index(drop = True)
        )
    else:
        result = pd.DataFrame(columns = ['kegg_gene_id', 'ko_id'])

    result.to_csv(output_file, sep = '\t', index = False)
    return result


def build_mmc3_member_annotation_table(
    membership_frame: pd.DataFrame,
    uniprot_annotation_frame: Optional[pd.DataFrame] = None,
    kegg_ko_frame: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """Combine local OMA labels with optional UniProt and KEGG annotations."""

    annotation_frame = membership_frame.copy()
    parsed_local = pd.DataFrame(
        [
            _parse_oma_member_label(
                microbial_protein = row['microbial_protein'],
                oma_member_label = row.get('oma_member_label', ''),
            )
            for _, row in annotation_frame.iterrows()
        ],
    )
    annotation_frame = pd.concat(
        [
            annotation_frame.reset_index(drop = True),
            parsed_local.reset_index(drop = True),
        ],
        axis = 1,
    )

    if uniprot_annotation_frame is not None and not uniprot_annotation_frame.empty:
        uniprot_subset = uniprot_annotation_frame.copy()
        uniprot_subset['microbial_protein'] = (
            uniprot_subset['microbial_protein'].astype(str).str.strip()
        )
        annotation_frame = annotation_frame.merge(
            uniprot_subset,
            on = 'microbial_protein',
            how = 'left',
        )
    else:
        for column_name in [
            'uniprot_protein_name',
            'uniprot_gene_names',
            'uniprot_organism_name',
            'ec_numbers',
            'refseq_accessions',
            'embl_protein_accessions',
            'kegg_gene_ids',
        ]:
            annotation_frame[column_name] = ''

    if kegg_ko_frame is not None and not kegg_ko_frame.empty:
        kegg_map = (
            kegg_ko_frame.groupby('kegg_gene_id')['ko_id']
            .apply(lambda values: ';'.join(unique_preserve_order(values)))
            .to_dict()
        )
        annotation_frame['ko_ids'] = annotation_frame['kegg_gene_ids'].map(
            lambda value: ';'.join(
                unique_preserve_order(
                    ko_id
                    for gene_id in split_semicolon_values(value)
                    for ko_id in split_semicolon_values(kegg_map.get(gene_id, ''))
                ),
            ),
        )
    else:
        annotation_frame['ko_ids'] = ''

    annotation_frame['protein_name'] = annotation_frame['uniprot_protein_name'].fillna(
        '',
    )
    annotation_frame['protein_name'] = annotation_frame['protein_name'].where(
        annotation_frame['protein_name'].astype(str).str.strip() != '',
        annotation_frame['local_protein_name'],
    )
    annotation_frame['gene_names'] = annotation_frame['uniprot_gene_names'].fillna('')
    annotation_frame['gene_names'] = annotation_frame['gene_names'].where(
        annotation_frame['gene_names'].astype(str).str.strip() != '',
        annotation_frame['local_gene_name'],
    )
    annotation_frame['organism_name'] = annotation_frame['uniprot_organism_name'].fillna(
        '',
    )
    annotation_frame['organism_name'] = annotation_frame['organism_name'].where(
        annotation_frame['organism_name'].astype(str).str.strip() != '',
        annotation_frame['local_organism_name'],
    )
    annotation_frame['protein_name_normalized'] = annotation_frame['protein_name'].map(
        normalize_text_for_matching,
    )

    return annotation_frame


def build_mmc3_member_token_table(
    annotation_frame: pd.DataFrame,
) -> pd.DataFrame:
    """Expand annotated OMA members into a token table for exact matching."""

    rows: list[dict[str, str]] = []

    for _, row in annotation_frame.iterrows():
        shared_fields = {
            'species': str(row['species']).strip(),
            'oma_group': str(row['oma_group']).strip(),
            'microbial_protein': str(row['microbial_protein']).strip(),
            'protein_name': str(row.get('protein_name', '')).strip(),
            'gene_names': str(row.get('gene_names', '')).strip(),
        }
        emitted: set[tuple[str, str]] = set()

        def emit(
            match_type: str,
            match_token: str,
            token_source: str,
        ) -> None:
            normalized_token = str(match_token).strip()

            if not normalized_token:
                return

            key = (match_type, normalized_token)

            if key in emitted:
                return

            emitted.add(key)
            rows.append(
                {
                    **shared_fields,
                    'match_type': match_type,
                    'match_token': normalized_token,
                    'token_source': token_source,
                },
            )

        emit(
            match_type = 'protein_accession',
            match_token = row['microbial_protein'],
            token_source = 'microbial_protein',
        )

        for value in split_semicolon_values(row.get('refseq_accessions', '')):
            emit(
                match_type = 'protein_accession',
                match_token = value,
                token_source = 'refseq_accession',
            )

        for value in split_semicolon_values(row.get('embl_protein_accessions', '')):
            emit(
                match_type = 'protein_accession',
                match_token = value,
                token_source = 'embl_protein_accession',
            )

        for value in split_semicolon_values(row.get('ec_numbers', '')):
            emit(
                match_type = 'ec_number',
                match_token = value,
                token_source = 'ec_number',
            )

        for value in split_semicolon_values(row.get('ko_ids', '')):
            emit(
                match_type = 'ko_id',
                match_token = value,
                token_source = 'ko_id',
            )

        normalized_protein_name = normalize_text_for_matching(
            row.get('protein_name', ''),
        )

        if normalized_protein_name:
            emit(
                match_type = 'protein_name',
                match_token = normalized_protein_name,
                token_source = 'protein_name',
            )

    if not rows:
        return pd.DataFrame(
            columns = [
                'species',
                'oma_group',
                'microbial_protein',
                'protein_name',
                'gene_names',
                'match_type',
                'match_token',
                'token_source',
            ],
        )

    return (
        pd.DataFrame(rows)
        .drop_duplicates()
        .sort_values(['oma_group', 'microbial_protein', 'match_type', 'match_token'])
        .reset_index(drop = True)
    )


def filter_oma_members_with_mmc3(
    membership_frame: pd.DataFrame,
    mmc3_feature_tokens: pd.DataFrame,
    uniprot_annotation_frame: Optional[pd.DataFrame] = None,
    kegg_ko_frame: Optional[pd.DataFrame] = None,
) -> Mmc3OrthogroupFilterResult:
    """Filter OMA members to groups supported by CD-increased mmc3 features."""

    member_annotations = build_mmc3_member_annotation_table(
        membership_frame = membership_frame,
        uniprot_annotation_frame = uniprot_annotation_frame,
        kegg_ko_frame = kegg_ko_frame,
    )
    member_tokens = build_mmc3_member_token_table(member_annotations)

    match_evidence = member_tokens.merge(
        mmc3_feature_tokens,
        on = ['match_type', 'match_token'],
        how = 'inner',
        suffixes = ('_member', '_mmc3'),
    )
    match_evidence = (
        match_evidence.drop_duplicates()
        .sort_values(
            [
                'oma_group',
                'microbial_protein',
                'set_name',
                'feature',
                'match_type',
                'match_token',
            ],
        )
        .reset_index(drop = True)
    )

    if match_evidence.empty:
        matched_groups = pd.DataFrame(
            columns = [
                'oma_group',
                'matched_member_count',
                'matched_species_count',
                'matched_feature_count',
                'matched_set_count',
                'matched_sets',
            ],
        )
        microbial_proteins = membership_frame.iloc[0:0][
            ['species', 'oma_group', 'microbial_protein']
        ].copy()
    else:
        matched_groups = (
            match_evidence.groupby('oma_group')
            .agg(
                matched_member_count = ('microbial_protein', 'nunique'),
                matched_species_count = ('species', 'nunique'),
                matched_feature_count = ('feature', 'nunique'),
                matched_set_count = ('set_name', 'nunique'),
                matched_sets = (
                    'set_name',
                    lambda values: ';'.join(unique_preserve_order(values)),
                ),
            )
            .reset_index()
            .sort_values(
                [
                    'matched_feature_count',
                    'matched_member_count',
                    'oma_group',
                ],
                ascending = [False, False, True],
            )
            .reset_index(drop = True)
        )
        microbial_proteins = (
            membership_frame[
                membership_frame['oma_group'].isin(matched_groups['oma_group'])
            ][['species', 'oma_group', 'microbial_protein']]
            .drop_duplicates()
            .sort_values(['species', 'oma_group', 'microbial_protein'])
            .reset_index(drop = True)
        )

    return Mmc3OrthogroupFilterResult(
        microbial_proteins = microbial_proteins,
        member_annotations = member_annotations,
        mmc3_feature_tokens = mmc3_feature_tokens,
        match_evidence = match_evidence,
        matched_groups = matched_groups,
    )


def prepare_microbial_protein_table_from_mmc3(
    membership_table_path: PathLike,
    mmc3_path: PathLike,
    membership_group_column: ColumnName = 'oma_group',
    membership_protein_column: ColumnName = 'microbial_protein',
    membership_species_column: ColumnName = 'species',
    membership_label_column: Optional[ColumnName] = 'oma_member_label',
    mmc3_set_column: ColumnName = 'Set',
    mmc3_feature_column: ColumnName = 'Feature',
    mmc3_effect_column: Optional[ColumnName] = 'Dysbiosis Coefficient (CD)',
    mmc3_allowed_sets: Optional[list[str]] = None,
    cd_increased_only: bool = True,
    uniprot_annotation_path: Optional[PathLike] = None,
    kegg_ko_mapping_path: Optional[PathLike] = None,
    download_uniprot_annotations: bool = False,
    download_kegg_ko_mapping: bool = False,
) -> Mmc3OrthogroupFilterResult:
    """Build the microbial protein table by filtering OMA members with mmc3."""

    membership_table = read_table(membership_table_path)
    resolved_group = resolve_column_name(
        membership_table,
        membership_group_column,
    )
    resolved_protein = resolve_column_name(
        membership_table,
        membership_protein_column,
    )
    resolved_species = resolve_column_name(
        membership_table,
        membership_species_column,
    )
    selected_columns = [
        resolved_species,
        resolved_group,
        resolved_protein,
    ]

    if membership_label_column is not None:
        resolved_label = resolve_column_name(
            membership_table,
            membership_label_column,
        )
        selected_columns.append(resolved_label)
    else:
        resolved_label = None

    membership_frame = membership_table[selected_columns].copy().rename(
        columns = {
            resolved_species: 'species',
            resolved_group: 'oma_group',
            resolved_protein: 'microbial_protein',
        },
    )

    if resolved_label is not None:
        membership_frame = membership_frame.rename(
            columns = {resolved_label: 'oma_member_label'},
        )
    else:
        membership_frame['oma_member_label'] = membership_frame['microbial_protein']

    membership_frame['species'] = membership_frame['species'].map(
        normalize_species_label,
    )
    membership_frame['oma_group'] = membership_frame['oma_group'].astype(str).str.strip()
    membership_frame['microbial_protein'] = (
        membership_frame['microbial_protein'].astype(str).str.strip()
    )
    membership_frame = membership_frame.drop_duplicates().reset_index(drop = True)

    mmc3_feature_tokens = parse_mmc3_feature_table(
        input_path = mmc3_path,
        set_column = mmc3_set_column,
        feature_column = mmc3_feature_column,
        effect_column = mmc3_effect_column,
        allowed_sets = mmc3_allowed_sets,
        cd_increased_only = cd_increased_only,
    )

    uniprot_annotation_frame: Optional[pd.DataFrame] = None
    uniprot_ids = membership_frame.loc[
        membership_frame['microbial_protein'].str.match(UNIPROT_ACCESSION_PATTERN),
        'microbial_protein',
    ].tolist()

    if uniprot_annotation_path is not None and Path(uniprot_annotation_path).exists():
        uniprot_annotation_frame = read_table(uniprot_annotation_path)
    elif download_uniprot_annotations and uniprot_annotation_path is not None:
        uniprot_annotation_frame = download_uniprot_json_annotation_table(
            uniprot_ids = uniprot_ids,
            output_path = uniprot_annotation_path,
        )

    kegg_ko_frame: Optional[pd.DataFrame] = None

    if kegg_ko_mapping_path is not None and Path(kegg_ko_mapping_path).exists():
        kegg_ko_frame = read_table(kegg_ko_mapping_path)
    elif (
        download_kegg_ko_mapping
        and kegg_ko_mapping_path is not None
        and uniprot_annotation_frame is not None
    ):
        kegg_gene_ids = unique_preserve_order(
            gene_id
            for value in uniprot_annotation_frame.get('kegg_gene_ids', [])
            for gene_id in split_semicolon_values(value)
        )
        kegg_ko_frame = download_kegg_gene_to_ko_table(
            kegg_gene_ids = kegg_gene_ids,
            output_path = kegg_ko_mapping_path,
        )

    return filter_oma_members_with_mmc3(
        membership_frame = membership_frame,
        mmc3_feature_tokens = mmc3_feature_tokens,
        uniprot_annotation_frame = uniprot_annotation_frame,
        kegg_ko_frame = kegg_ko_frame,
    )


def query_current_human_candidate_accessions(
    location_filters: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Query the current OmniPath intercell resource for host candidates."""

    import omnipath as op

    resolved_filters = location_filters or DEFAULT_LOCATION_FILTERS
    frame = op.requests.Intercell.get(
        parent = resolved_filters,
        scope = ['generic', 'specific'],
        source = ['resource_specific', 'composite'],
        entity_type = 'protein',
    )

    if 'uniprot' not in frame.columns:
        raise ValueError(
            'The OmniPath intercell table does not expose a `uniprot` column.',
        )

    optional_columns = [
        column_name
        for column_name in [
            'uniprot',
            'genesymbol',
            'parent',
            'category',
            'database',
            'source',
        ]
        if column_name in frame.columns
    ]

    result = frame[optional_columns].copy()
    result['uniprot'] = result['uniprot'].astype(str).str.strip()
    result = result[result['uniprot'].str.match(UNIPROT_ACCESSION_PATTERN)]

    return (
        result
        .drop_duplicates()
        .sort_values(optional_columns)
        .reset_index(drop = True)
    )


def write_identifier_table(
    values: Iterable[str],
    output_path: PathLike,
    column_name: str = 'uniprot',
) -> Path:
    """Write one identifier per row as a TSV table."""

    output_file = Path(output_path)
    output_file.parent.mkdir(parents = True, exist_ok = True)
    identifiers = unique_preserve_order(values)

    pd.DataFrame({column_name: identifiers}).to_csv(
        output_file,
        sep = '\t',
        index = False,
    )
    return output_file


def _build_uniprot_query(uniprot_ids: list[str]) -> str:
    """Build a UniProt OR query for an accession batch."""

    return 'accession:(' + ' OR '.join(uniprot_ids) + ')'


def _merge_split_uniprot_stream_responses(
    format_name: str,
    response_chunks: list[str],
) -> str:
    """Merge recursively split UniProt responses back into one payload."""

    non_empty_chunks = [
        chunk
        for chunk in response_chunks
        if str(chunk).strip()
    ]

    if not non_empty_chunks:
        return ''

    if format_name != 'tsv':
        return ''.join(non_empty_chunks)

    merged_lines: list[str] = []

    for index, chunk in enumerate(non_empty_chunks):
        lines = chunk.splitlines()

        if index > 0 and lines:
            lines = lines[1:]

        merged_lines.extend(lines)

    if not merged_lines:
        return ''

    return '\n'.join(merged_lines) + '\n'


def _request_uniprot_stream(
    uniprot_ids: list[str],
    format_name: str,
    fields: Optional[list[str]] = None,
    timeout: int = 120,
) -> str:
    """Request a FASTA or TSV stream from UniProt."""

    params: dict[str, str] = {
        'format': format_name,
        'query': _build_uniprot_query(uniprot_ids),
    }

    if fields is not None:
        params['fields'] = ','.join(fields)

    try:
        response = requests.get(
            'https://rest.uniprot.org/uniprotkb/stream',
            params = params,
            timeout = timeout,
        )
        response.raise_for_status()
        return response.text
    except requests.HTTPError as error:
        status_code = (
            error.response.status_code
            if error.response is not None
            else None
        )

        if status_code in {400, 414} and len(uniprot_ids) > 1:
            midpoint = len(uniprot_ids) // 2
            left_text = _request_uniprot_stream(
                uniprot_ids[:midpoint],
                format_name = format_name,
                fields = fields,
                timeout = timeout,
            )
            right_text = _request_uniprot_stream(
                uniprot_ids[midpoint:],
                format_name = format_name,
                fields = fields,
                timeout = timeout,
            )
            return _merge_split_uniprot_stream_responses(
                format_name = format_name,
                response_chunks = [left_text, right_text],
            )

        raise


def download_uniprot_fasta(
    uniprot_ids: Iterable[str],
    output_path: PathLike,
    batch_size: int = DEFAULT_FASTA_BATCH_SIZE,
) -> Path:
    """Download a FASTA file for a list of UniProt accessions."""

    output_file = Path(output_path)
    output_file.parent.mkdir(parents = True, exist_ok = True)
    ordered_ids = unique_preserve_order(uniprot_ids)

    with open(output_file, 'w', encoding = 'utf-8') as fasta_file:
        for start in range(0, len(ordered_ids), batch_size):
            batch = ordered_ids[start : start + batch_size]
            fasta_file.write(
                _request_uniprot_stream(
                    batch,
                    format_name = 'fasta',
                ),
            )

    return output_file


def download_uniprot_annotation_table(
    uniprot_ids: Iterable[str],
    output_path: PathLike,
    fields: Optional[list[str]] = None,
    batch_size: int = DEFAULT_ANNOTATION_BATCH_SIZE,
) -> pd.DataFrame:
    """Download a UniProt TSV annotation table for a list of accessions."""

    output_file = Path(output_path)
    output_file.parent.mkdir(parents = True, exist_ok = True)
    ordered_ids = unique_preserve_order(uniprot_ids)
    frames: list[pd.DataFrame] = []

    for start in range(0, len(ordered_ids), batch_size):
        batch = ordered_ids[start : start + batch_size]
        response_text = _request_uniprot_stream(
            batch,
            format_name = 'tsv',
            fields = fields or DEFAULT_UNIPROT_FIELDS,
        )
        frame = pd.read_csv(StringIO(response_text), sep = '\t')
        frames.append(frame)

    if frames:
        result = pd.concat(frames, ignore_index = True).drop_duplicates()
    else:
        result = pd.DataFrame()

    result.to_csv(output_file, sep = '\t', index = False)
    return result


def _first_present_column(
    frame: pd.DataFrame,
    candidates: list[str],
) -> str:
    """Return the first column name present in the frame."""

    for candidate in candidates:
        if candidate in frame.columns:
            return candidate

    raise ValueError(
        f'Could not find any of {candidates} in {list(frame.columns)}.',
    )


def _first_present_optional_column(
    frame: pd.DataFrame,
    candidates: list[str],
) -> Optional[str]:
    """Return the first matching column present in the frame, if any."""

    for candidate in candidates:
        if candidate in frame.columns:
            return candidate

    return None


def uniprot_annotation_table_to_domain_table(
    annotation_frame: pd.DataFrame,
) -> pd.DataFrame:
    """Convert a raw UniProt TSV response to the MicrobioLink domain format."""

    accession_column = _first_present_column(
        annotation_frame,
        ['Entry', 'accession', 'Accession'],
    )
    pfam_column = _first_present_column(
        annotation_frame,
        ['Pfam', 'Cross-reference (Pfam)', 'xref_pfam'],
    )

    rows: list[dict[str, str]] = []

    for _, row in annotation_frame.iterrows():
        accession = str(row[accession_column]).strip()
        pfams = unique_preserve_order(
            PFAM_CAPTURE_PATTERN.findall(str(row[pfam_column])),
        )

        if not accession or not pfams:
            continue

        rows.append(
            {
                'protein': accession,
                'pfam': ';'.join(pfams),
            },
        )

    result = pd.DataFrame(rows)

    if result.empty:
        return pd.DataFrame(columns = ['protein', 'pfam'])

    return (
        result
        .drop_duplicates()
        .sort_values(['protein', 'pfam'])
        .reset_index(drop = True)
    )


def download_uniprot_domain_table(
    uniprot_ids: Iterable[str],
    raw_output_path: PathLike,
    domain_output_path: PathLike,
    fields: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Download and normalize UniProt Pfam annotations."""

    annotation_frame = download_uniprot_annotation_table(
        uniprot_ids = uniprot_ids,
        output_path = raw_output_path,
        fields = fields,
    )
    domain_frame = uniprot_annotation_table_to_domain_table(annotation_frame)
    domain_frame.to_csv(domain_output_path, sep = '\t', index = False)
    return domain_frame


def standardize_deg_table(
    input_path: PathLike,
    gene_column: ColumnName,
    value_column: ColumnName,
    pvalue_column: ColumnName,
    sheet_name: SheetName = 0,
) -> pd.DataFrame:
    """Select and rename the gene, fold-change, and p-value columns."""

    frame = read_table(input_path, sheet_name = sheet_name)
    resolved_gene = resolve_column_name(frame, gene_column)
    resolved_value = resolve_column_name(frame, value_column)
    resolved_pvalue = resolve_column_name(frame, pvalue_column)

    result = frame[
        [resolved_gene, resolved_value, resolved_pvalue]
    ].copy().rename(
        columns = {
            resolved_gene: 'gene_symbol',
            resolved_value: 'log2FC',
            resolved_pvalue: 'padj',
        },
    )
    result['gene_symbol'] = result['gene_symbol'].astype(str).str.strip()
    result['log2FC'] = pd.to_numeric(result['log2FC'], errors = 'coerce')
    result['padj'] = pd.to_numeric(result['padj'], errors = 'coerce')

    return (
        result
        .dropna(subset = ['gene_symbol', 'log2FC', 'padj'])
        .drop_duplicates()
        .reset_index(drop = True)
    )


def normalize_tiedie_hmi_table(
    interaction_frame: pd.DataFrame,
) -> pd.DataFrame:
    """Rename an interaction table to the TieDIE-compatible HMI schema."""

    human_column = _first_present_column(
        interaction_frame,
        TIEDIE_HUMAN_COLUMN_CANDIDATES,
    )
    bacterial_column = _first_present_column(
        interaction_frame,
        TIEDIE_BACTERIAL_COLUMN_CANDIDATES,
    )
    sign_column = _first_present_optional_column(
        interaction_frame,
        TIEDIE_SIGN_COLUMN_CANDIDATES,
    )

    selected_columns = [human_column, bacterial_column]
    rename_map = {
        human_column: '# Human Protein',
        bacterial_column: 'Bacterial protein',
    }

    if sign_column is not None:
        selected_columns.append(sign_column)
        rename_map[sign_column] = 'sign'

    normalized = interaction_frame[selected_columns].copy().rename(
        columns = rename_map,
    )
    normalized['# Human Protein'] = (
        normalized['# Human Protein'].astype(str).str.strip()
    )
    normalized['Bacterial protein'] = (
        normalized['Bacterial protein'].astype(str).str.strip()
    )

    if 'sign' in normalized.columns:
        normalized['sign'] = normalized['sign'].astype(str).str.strip()

    return normalized.drop_duplicates().reset_index(drop = True)


def write_tiedie_hmi_table(
    interaction_frame: pd.DataFrame,
    output_path: PathLike,
) -> Path:
    """Write a TieDIE-compatible HMI table to disk."""

    output_file = Path(output_path)
    output_file.parent.mkdir(parents = True, exist_ok = True)
    normalize_tiedie_hmi_table(interaction_frame).to_csv(
        output_file,
        sep = '\t',
        index = False,
    )
    return output_file


def load_hmi_table_for_tiedie(
    file_path: PathLike,
    repo_root: Optional[Path] = None,
) -> pd.DataFrame:
    """Load and harmonize an HMI table for TieDIE processing."""

    resolved_repo_root = repo_root or find_repo_root()
    helper_module = load_module_from_path(
        'microbiolink_tiedie_hmi_external',
        resolved_repo_root / 'microbiolink' / '_tiedie_hmi.py',
    )
    return helper_module.load_hmi_table(str(file_path))


def _consensus_direction_label(value: object) -> str:
    """Convert OmniPath consensus values to TieDIE edge labels."""

    if pd.isna(value):
        return 'unknown'

    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {
            '1',
            'true',
            'activates',
            'activates>',
            'stimulates',
            'stimulates>',
        }:
            return 'activates>'
        if normalized in {
            '0',
            'false',
            'inhibits',
            'inhibits>',
        }:
            return 'inhibits>'
        return 'unknown'

    if value is True or value == 1:
        return 'activates>'

    if value is False or value == 0:
        return 'inhibits>'

    return 'unknown'


def read_endpoint_gene_table(
    endpoint_file: PathLike,
    separator: str = '\t',
    pvalue_column: Optional[int] = 3,
) -> pd.DataFrame:
    """Read a DEG table and rename its first column to `target_genesymbol`."""

    endpoint_frame = pd.read_csv(endpoint_file, sep = separator).copy()
    endpoint_frame = endpoint_frame.rename(
        columns = {endpoint_frame.columns[0]: 'target_genesymbol'},
    )
    endpoint_frame['target_genesymbol'] = (
        endpoint_frame['target_genesymbol'].astype(str).str.strip()
    )

    if pvalue_column is not None:
        resolved_pvalue_column = endpoint_frame.columns[pvalue_column - 1]
        endpoint_frame[resolved_pvalue_column] = pd.to_numeric(
            endpoint_frame[resolved_pvalue_column],
            errors = 'coerce',
        )
        endpoint_frame = endpoint_frame[
            endpoint_frame[resolved_pvalue_column] < 0.05
        ]

    return endpoint_frame.drop_duplicates().reset_index(drop = True)


def build_tiedie_inputs_from_omnipath_network(
    endpoint_file: PathLike,
    hmi_prediction_file: PathLike,
    output_dir: PathLike,
    repo_root: Optional[Path] = None,
    endpoint_separator: str = '\t',
    endpoint_pvalue_column: Optional[int] = 3,
    endpoint_value_column: int = 2,
    pathway_input_filename: str = 'pathway.sif',
    upstream_input_filename: str = 'upstream.input',
    downstream_input_filename: str = 'downstream.input',
) -> dict[str, Path]:
    """Build TieDIE inputs from the full directed OmniPath network and DEGs."""

    import omnipath as op

    resolved_output_dir = ensure_directory(output_dir)
    ppis = op.interactions.OmniPath.get(genesymbols = 1).copy()
    tf_tg = op.interactions.Transcriptional.get(
        databases = 'CollecTRI',
        genesymbols = 1,
    ).copy()
    endpoint_genes = read_endpoint_gene_table(
        endpoint_file,
        separator = endpoint_separator,
        pvalue_column = endpoint_pvalue_column,
    )

    contextualised_tf_tg = tf_tg[
        tf_tg['target_genesymbol'].isin(endpoint_genes['target_genesymbol'])
    ].drop_duplicates().copy()

    if 'source_genesymbol' in contextualised_tf_tg.columns:
        contextualised_tf_tg['source'] = contextualised_tf_tg['source_genesymbol']

    if 'target_genesymbol' in contextualised_tf_tg.columns:
        contextualised_tf_tg['target'] = contextualised_tf_tg['target_genesymbol']

    regulators = contextualised_tf_tg.pivot_table(
        columns = ['source'],
        aggfunc = 'size',
    )
    regulators = pd.DataFrame(regulators).reset_index().rename(
        columns = {'source': 'TF_name', 0: 'num_degs'},
    )

    contextualised_tf_tg_path = resolved_output_dir / (
        'contextualised_regulator-target_network.txt'
    )
    regulators_path = resolved_output_dir / (
        'contextualised_regulators_of_targets.txt'
    )
    contextualised_tf_tg.to_csv(
        contextualised_tf_tg_path,
        sep = '\t',
        index = False,
    )
    regulators.to_csv(
        regulators_path,
        sep = '\t',
        index = False,
    )

    source_column = (
        'source_genesymbol'
        if 'source_genesymbol' in ppis.columns
        else 'source'
    )
    target_column = (
        'target_genesymbol'
        if 'target_genesymbol' in ppis.columns
        else 'target'
    )
    pathway_frame = ppis[[source_column, target_column]].copy().rename(
        columns = {
            source_column: 'source',
            target_column: 'target',
        },
    )
    pathway_frame['direction'] = ppis['consensus_stimulation'].map(
        _consensus_direction_label,
    )
    pathway_frame = pathway_frame.dropna(subset = ['source', 'target'])
    pathway_frame = pathway_frame[['source', 'direction', 'target']].drop_duplicates()
    pathway_path = resolved_output_dir / pathway_input_filename
    pathway_frame.to_csv(
        pathway_path,
        sep = '\t',
        index = False,
        header = False,
    )

    hmi_frame = load_hmi_table_for_tiedie(
        hmi_prediction_file,
        repo_root = repo_root,
    )
    if 'sign' in hmi_frame.columns:
        filtered_hmi = hmi_frame[hmi_frame['sign'].isin(['-', '+'])].copy()
        upstream_frame = filtered_hmi[
            ['# Human Protein', 'sign']
        ].rename(columns = {'sign': 'direction'})
        upstream_frame = (
            upstream_frame
            .groupby(['# Human Protein', 'direction'])
            .size()
            .reset_index(name = 'n')
        )
    else:
        upstream_frame = (
            hmi_frame[['# Human Protein']]
            .groupby('# Human Protein')
            .size()
            .reset_index(name = 'n')
        )
        upstream_frame['direction'] = '-'

    upstream_path = resolved_output_dir / upstream_input_filename
    upstream_frame[['# Human Protein', 'n', 'direction']].to_csv(
        upstream_path,
        sep = '\t',
        index = False,
        header = False,
    )

    expression_column = endpoint_genes.columns[endpoint_value_column - 1]
    endpoint_numeric = endpoint_genes.copy()
    endpoint_numeric[expression_column] = pd.to_numeric(
        endpoint_numeric[expression_column],
        errors = 'coerce',
    )
    downstream_join = contextualised_tf_tg.merge(
        endpoint_numeric,
        on = 'target_genesymbol',
        how = 'inner',
    )
    downstream_subset = downstream_join[
        ['source', 'target', 'consensus_stimulation', expression_column]
    ].copy()
    downstream_subset['exp_sign'] = np.where(
        downstream_subset['consensus_stimulation'] == True,
        downstream_subset[expression_column],
        -1 * downstream_subset[expression_column],
    )
    grouped = downstream_subset.groupby('source')
    downstream_frame = grouped.count()
    downstream_frame['sumof'] = grouped.exp_sign.sum()
    downstream_frame['final_val'] = (
        downstream_frame['sumof'] / downstream_frame['exp_sign']
    )
    downstream_frame['sign'] = np.where(
        downstream_frame['final_val'] >= 0,
        '+',
        '-',
    )
    downstream_frame = downstream_frame.drop(
        ['exp_sign', 'sumof', 'target'],
        axis = 1,
    )
    downstream_export = downstream_frame.reset_index()[['source', 'final_val', 'sign']]
    downstream_path = resolved_output_dir / downstream_input_filename
    downstream_export.to_csv(
        downstream_path,
        sep = '\t',
        index = False,
        header = False,
    )

    return {
        'pathway': pathway_path,
        'upstream': upstream_path,
        'downstream': downstream_path,
        'contextualised_tf_tg': contextualised_tf_tg_path,
        'regulators': regulators_path,
    }


def build_reverse_monte_carlo_input_table(
    reverse_interactions: pd.DataFrame,
) -> pd.DataFrame:
    """Reshape reverse DMI output for the generic Monte Carlo filter."""

    required_columns = [
        'microbial_protein',
        'motif',
        'start',
        'end',
    ]
    missing_columns = [
        column_name
        for column_name in required_columns
        if column_name not in reverse_interactions.columns
    ]

    if missing_columns:
        raise ValueError(
            f'Reverse DMI frame is missing columns: {missing_columns}.',
        )

    return (
        reverse_interactions[required_columns]
        .rename(columns = {'microbial_protein': 'human_protein'})
        .drop_duplicates()
        .sort_values(['human_protein', 'motif', 'start', 'end'])
        .reset_index(drop = True)
    )


def merge_reverse_monte_carlo_results(
    reverse_interactions: pd.DataFrame,
    reverse_monte_carlo_results: pd.DataFrame,
) -> pd.DataFrame:
    """Merge reverse-DMI Monte Carlo annotations back onto the full table."""

    renamed_results = reverse_monte_carlo_results.rename(
        columns = {'human_protein': 'microbial_protein'},
    )
    merge_columns = ['microbial_protein', 'motif', 'start', 'end']

    return reverse_interactions.merge(
        renamed_results.drop_duplicates(subset = merge_columns),
        on = merge_columns,
        how = 'left',
    )


def normalize_forward_interactions_for_tiedie(
    interaction_frame: pd.DataFrame,
    evidence_label: str = 'forward_dmi',
) -> pd.DataFrame:
    """Normalize forward-DMI style tables for TieDIE input building."""

    human_column = _first_present_column(
        interaction_frame,
        ['human_protein', '# Human Protein', 'Human Protein'],
    )
    bacterial_column = _first_present_column(
        interaction_frame,
        ['bacterial_protein', 'Bacterial protein', 'Bacteria Protein'],
    )

    result = interaction_frame[[human_column, bacterial_column]].copy().rename(
        columns = {
            human_column: 'human_protein',
            bacterial_column: 'bacterial_protein',
        },
    )
    result['evidence_type'] = evidence_label

    return result.drop_duplicates().reset_index(drop = True)


def normalize_reverse_interactions_for_tiedie(
    interaction_frame: pd.DataFrame,
    evidence_label: str = 'reverse_dmi',
) -> pd.DataFrame:
    """Normalize reverse-DMI output for TieDIE input building."""

    result = interaction_frame[
        ['host_protein', 'microbial_protein']
    ].copy().rename(
        columns = {
            'host_protein': 'human_protein',
            'microbial_protein': 'bacterial_protein',
        },
    )
    result['evidence_type'] = evidence_label

    return result.drop_duplicates().reset_index(drop = True)


def normalize_ddi_interactions_for_tiedie(
    interaction_frame: pd.DataFrame,
    evidence_label: str = 'ddi',
) -> pd.DataFrame:
    """Normalize DDI output for TieDIE input building."""

    result = interaction_frame[
        ['human_protein', 'bacterial_protein']
    ].copy()
    result['evidence_type'] = evidence_label

    return result.drop_duplicates().reset_index(drop = True)


def combine_tiedie_upstream_tables(
    upstream_tables: list[pd.DataFrame],
    deduplicate_pairs: bool = True,
) -> pd.DataFrame:
    """Combine forward, reverse, and DDI upstream seed tables."""

    combined = pd.concat(upstream_tables, ignore_index = True)
    combined = combined.dropna(subset = ['human_protein', 'bacterial_protein'])
    combined['human_protein'] = combined['human_protein'].astype(str).str.strip()
    combined['bacterial_protein'] = (
        combined['bacterial_protein'].astype(str).str.strip()
    )
    combined['evidence_type'] = combined['evidence_type'].astype(str).str.strip()

    if not deduplicate_pairs:
        return combined.drop_duplicates().reset_index(drop = True)

    aggregated = combined.groupby(
        ['human_protein', 'bacterial_protein'],
        as_index = False,
    ).agg(
        {
            'evidence_type': lambda values: '|'.join(
                unique_preserve_order(values),
            ),
        },
    )

    return aggregated.sort_values(
        ['human_protein', 'bacterial_protein'],
    ).reset_index(drop = True)


def run_command(
    command: list[str],
    cwd: Optional[PathLike] = None,
    env: Optional[dict[str, str]] = None,
) -> subprocess.CompletedProcess[str]:
    """Run a subprocess command and echo it in notebook-friendly form."""

    print('$', shlex.join(command))
    return subprocess.run(
        command,
        cwd = cwd,
        env = env,
        check = True,
        text = True,
        capture_output = False,
    )
