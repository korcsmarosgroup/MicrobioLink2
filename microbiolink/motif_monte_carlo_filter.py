#!/usr/bin/env python3

"""Monte Carlo motif filter based on shuffled disordered-region placements."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from microbiolink.DMI import extract_uniprot_id
from microbiolink.DMI import read_fasta_sequences


PathLike = str | Path

INTERACTION_COLUMN_CANDIDATES = {
    'human_protein': ['humanprotein'],
    'motif': ['motif'],
    'start': ['start'],
    'end': ['end'],
}

STRUCTURE_COLUMN_CANDIDATES = {
    'human_protein': ['humanprotein', 'protein', 'uniprot'],
    'position': ['position', 'residueposition', 'pos'],
    'disorder_score': ['disorderscore', 'iupredscore', 'disorder'],
    'binding_score': ['bindingscore', 'anchorscore', 'binding', 'anchor'],
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(
        description=(
            'Filter motif-mediated host-microbe interactions with a Monte Carlo '
            'test that shuffles disordered-region placements across each host '
            'protein while preserving region lengths.'
        ),
    )
    parser.add_argument(
        '--interaction_file',
        required = True,
        help = (
            'Path to a DMI interaction table. The script accepts both the '
            'semicolon-separated legacy DMI output and structurally annotated '
            'tables as long as they contain human protein, motif, start, and '
            'end columns.'
        ),
    )
    parser.add_argument(
        '--output_file',
        required = True,
        help = 'Path to the full output table with Monte Carlo statistics.',
    )
    parser.add_argument(
        '--filtered_output_file',
        help = 'Optional path to the subset of interactions passing the filter.',
    )
    parser.add_argument(
        '--structure_scores',
        help = (
            'Optional table of residue-level structural scores. Expected '
            'columns include protein, position, disorder score, and optionally '
            'binding score.'
        ),
    )
    parser.add_argument(
        '--fasta_file',
        help = (
            'Optional human FASTA file. Required when `--structure_scores` is '
            'not provided, because the script then predicts structure tracks '
            'with AIUPred.'
        ),
    )
    parser.add_argument(
        '--coordinate_system',
        choices = ['zero_based_half_open', 'one_based_closed'],
        default = 'zero_based_half_open',
        help = (
            'Interpretation of motif coordinates in the interaction table. '
            'Default: zero-based inclusive start with exclusive end.'
        ),
    )
    parser.add_argument(
        '--iterations',
        type = int,
        default = 1000,
        help = 'Number of Monte Carlo shuffles per motif. Default: 1000.',
    )
    parser.add_argument(
        '--alpha',
        type = float,
        default = 0.05,
        help = 'Significance threshold for keeping motifs. Default: 0.05.',
    )
    parser.add_argument(
        '--disorder_threshold',
        type = float,
        default = 0.60,
        help = 'Residue disorder-score threshold. Default: 0.60.',
    )
    parser.add_argument(
        '--binding_threshold',
        type = float,
        default = 0.60,
        help = (
            'Residue binding-score threshold used when `--require_binding` is '
            'set. Default: 0.60.'
        ),
    )
    parser.add_argument(
        '--min_support_fraction',
        type = float,
        default = 1.0,
        help = (
            'Minimum fraction of motif residues that must overlap shuffled '
            'supportive regions for a motif to pass. Default: 1.0.'
        ),
    )
    parser.add_argument(
        '--require_binding',
        action = 'store_true',
        help = (
            'Require both disorder and binding support when defining '
            'supportive structural regions. By default, only disorder is used.'
        ),
    )
    parser.add_argument(
        '--seed',
        type = int,
        default = 0,
        help = 'Random seed for reproducible Monte Carlo sampling. Default: 0.',
    )
    parser.add_argument(
        '--gpu',
        default = 0,
        help = 'GPU index for AIUPred when predicting structure tracks.',
    )
    parser.add_argument(
        '--force_cpu',
        action = 'store_true',
        help = 'Force AIUPred to use the CPU when predicting structure tracks.',
    )
    return parser.parse_args()


def _normalize_column_name(column_name: str) -> str:
    """Normalize a column name for flexible matching."""

    return re.sub(r'[^a-z0-9]+', '', column_name.lower().lstrip('#'))


def _read_delimited_table(filename: PathLike) -> pd.DataFrame:
    """Read a delimited table with automatic separator detection."""

    return pd.read_csv(
        filename,
        sep = None,
        engine = 'python',
    )


def _resolve_columns(
    frame: pd.DataFrame,
    candidates: dict[str, list[str]],
) -> dict[str, str]:
    """Resolve required logical columns to physical data-frame columns."""

    normalized_columns = {
        _normalize_column_name(column_name): column_name
        for column_name in frame.columns
    }
    resolved: dict[str, str] = {}

    for logical_name, aliases in candidates.items():
        for alias in aliases:
            if alias in normalized_columns:
                resolved[logical_name] = normalized_columns[alias]
                break
        else:
            if logical_name == 'binding_score':
                continue
            raise ValueError(
                f'Could not resolve a `{logical_name}` column in {list(frame.columns)}.',
            )

    return resolved


def read_interaction_table(filename: PathLike) -> pd.DataFrame:
    """Read an interaction table and normalize core column names."""

    frame = _read_delimited_table(filename)
    resolved = _resolve_columns(frame, INTERACTION_COLUMN_CANDIDATES)

    normalized = frame.copy()
    for logical_name, column_name in resolved.items():
        normalized[logical_name] = normalized[column_name]

    return normalized


def read_structure_score_table(filename: PathLike) -> dict[str, dict[str, np.ndarray | None]]:
    """Read residue-level structural scores from a table."""

    frame = _read_delimited_table(filename)
    resolved = _resolve_columns(frame, STRUCTURE_COLUMN_CANDIDATES)
    binding_column = resolved.get('binding_score')
    profiles: dict[str, dict[str, np.ndarray | None]] = {}

    grouped = frame.groupby(resolved['human_protein'], sort = False)
    for protein_name, group in grouped:
        max_position = int(pd.to_numeric(group[resolved['position']]).max())
        disorder = np.full(max_position, np.nan, dtype = float)
        binding = (
            np.full(max_position, np.nan, dtype = float)
            if binding_column is not None
            else None
        )

        for _, row in group.iterrows():
            position = int(row[resolved['position']]) - 1
            if position < 0:
                continue

            disorder[position] = float(row[resolved['disorder_score']])
            if binding is not None:
                binding[position] = float(row[binding_column])

        profiles[str(protein_name)] = {
            'disorder': disorder,
            'binding': binding,
        }

    return profiles


def build_sequence_dictionary(fasta_file: PathLike) -> dict[str, str]:
    """Build a UniProt-indexed sequence dictionary from a FASTA file."""

    sequences = read_fasta_sequences(fasta_file)
    return {
        extract_uniprot_id(header): sequence
        for header, sequence in sequences.items()
    }


def predict_structure_profiles_from_sequences(
    sequences: dict[str, str],
    force_cpu: bool,
    gpu: str | int,
) -> dict[str, dict[str, np.ndarray]]:
    """Predict residue-level disorder and binding tracks with AIUPred."""

    try:
        from iupred import init_aiupred_models
        from iupred import predict_aiupred_binding
        from iupred import predict_aiupred_disorder
    except ImportError as error:
        raise ImportError(
            'AIUPred support is not available. Install the `idr` extra or '
            'provide `--structure_scores` instead.',
        ) from error

    disorder_encoder, disorder_decoder, device = init_aiupred_models(
        'disorder',
        force_cpu = force_cpu,
        gpu_num = gpu,
    )
    binding_encoder, binding_decoder, _ = init_aiupred_models(
        'binding',
        force_cpu = force_cpu,
        gpu_num = gpu,
    )

    profiles: dict[str, dict[str, np.ndarray]] = {}
    for protein_name, sequence in sequences.items():
        disorder_profile = np.asarray(
            predict_aiupred_disorder(
                sequence,
                disorder_encoder,
                disorder_decoder,
                device,
                smoothing = True,
            ),
            dtype = float,
        )
        binding_profile = np.asarray(
            predict_aiupred_binding(
                sequence,
                binding_encoder,
                binding_decoder,
                device,
                smoothing = True,
                binding = True,
            ),
            dtype = float,
        )
        profiles[protein_name] = {
            'disorder': disorder_profile,
            'binding': binding_profile,
        }

    return profiles


def load_structure_profiles(
    interaction_frame: pd.DataFrame,
    structure_scores: PathLike | None,
    fasta_file: PathLike | None,
    force_cpu: bool,
    gpu: str | int,
) -> dict[str, dict[str, np.ndarray | None]]:
    """Load or predict per-residue structural profiles."""

    if structure_scores is not None:
        return read_structure_score_table(structure_scores)

    if fasta_file is None:
        raise ValueError(
            'Provide either `--structure_scores` or `--fasta_file`.',
        )

    proteins = {
        str(protein_name)
        for protein_name in interaction_frame['human_protein'].astype(str)
    }
    sequence_dictionary = build_sequence_dictionary(fasta_file)
    selected_sequences = {
        protein_name: sequence_dictionary[protein_name]
        for protein_name in proteins
        if protein_name in sequence_dictionary
    }
    return predict_structure_profiles_from_sequences(
        selected_sequences,
        force_cpu = force_cpu,
        gpu = gpu,
    )


def _convert_coordinates(
    start: Any,
    end: Any,
    coordinate_system: str,
) -> tuple[int, int]:
    """Convert motif coordinates to zero-based half-open indexing."""

    start_int = int(start)
    end_int = int(end)

    if coordinate_system == 'zero_based_half_open':
        return start_int, end_int

    return start_int - 1, end_int


def _support_mask(
    disorder_profile: np.ndarray,
    binding_profile: np.ndarray | None,
    disorder_threshold: float,
    binding_threshold: float,
    require_binding: bool,
) -> np.ndarray:
    """Build a boolean mask of structurally supportive residues."""

    mask = np.isfinite(disorder_profile) & (disorder_profile >= disorder_threshold)
    if require_binding:
        if binding_profile is None:
            raise ValueError(
                'Binding support was requested, but no binding-score column was available.',
            )
        mask &= np.isfinite(binding_profile) & (binding_profile >= binding_threshold)

    return mask


def _segment_lengths(mask: np.ndarray) -> list[int]:
    """Return contiguous True-segment lengths from a boolean mask."""

    lengths: list[int] = []
    current_length = 0

    for value in mask.tolist():
        if value:
            current_length += 1
            continue

        if current_length:
            lengths.append(current_length)
            current_length = 0

    if current_length:
        lengths.append(current_length)

    return lengths


def _random_composition(
    total: int,
    parts: int,
    rng: np.random.Generator,
) -> list[int]:
    """Sample a random composition of `total` into `parts` nonnegative values."""

    if parts == 1:
        return [total]

    bars = np.sort(
        rng.choice(
            total + parts - 1,
            size = parts - 1,
            replace = False,
        ),
    )

    values: list[int] = []
    previous = -1
    for bar in np.concatenate([bars, np.array([total + parts - 1])]):
        values.append(int(bar - previous - 1))
        previous = int(bar)

    return values


def shuffle_segment_mask(
    sequence_length: int,
    segment_lengths: list[int],
    rng: np.random.Generator,
) -> np.ndarray:
    """Shuffle segment placements while preserving segment lengths."""

    if not segment_lengths:
        return np.zeros(sequence_length, dtype = bool)

    total_segment_length = sum(segment_lengths)
    if total_segment_length > sequence_length:
        raise ValueError('Segment lengths exceed sequence length.')

    shuffled_lengths = np.array(segment_lengths, dtype = int)
    rng.shuffle(shuffled_lengths)

    gap_lengths = _random_composition(
        sequence_length - total_segment_length,
        len(shuffled_lengths) + 1,
        rng,
    )

    mask = np.zeros(sequence_length, dtype = bool)
    cursor = gap_lengths[0]

    for index, segment_length in enumerate(shuffled_lengths.tolist()):
        mask[cursor:cursor + segment_length] = True
        cursor += segment_length + gap_lengths[index + 1]

    return mask


def _nanmean_or_nan(values: np.ndarray) -> float:
    """Return the nanmean of an array or NaN when no finite values exist."""

    finite_values = values[np.isfinite(values)]
    if finite_values.size == 0:
        return math.nan

    return float(np.mean(finite_values))


def compute_motif_statistics(
    disorder_profile: np.ndarray,
    binding_profile: np.ndarray | None,
    support_mask: np.ndarray,
    start: int,
    end: int,
) -> dict[str, float]:
    """Compute observed motif statistics for one motif window."""

    motif_mask = support_mask[start:end]
    motif_length = end - start

    return {
        'sequence_length': float(disorder_profile.size),
        'motif_length': float(motif_length),
        'supportive_residue_count': float(np.sum(support_mask)),
        'supportive_residue_fraction': float(np.mean(support_mask)) if support_mask.size else math.nan,
        'observed_support_fraction': float(np.mean(motif_mask)) if motif_mask.size else math.nan,
        'observed_avg_disorder': _nanmean_or_nan(disorder_profile[start:end]),
        'observed_avg_binding': (
            _nanmean_or_nan(binding_profile[start:end])
            if binding_profile is not None
            else math.nan
        ),
    }


def monte_carlo_pvalue(
    sequence_length: int,
    segment_lengths: list[int],
    start: int,
    end: int,
    observed_support_fraction: float,
    iterations: int,
    rng: np.random.Generator,
) -> tuple[int, float]:
    """Estimate a motif p-value by shuffled supportive-region placements."""

    hits = 0
    for _ in range(iterations):
        shuffled_mask = shuffle_segment_mask(
            sequence_length,
            segment_lengths,
            rng,
        )
        shuffled_support_fraction = float(np.mean(shuffled_mask[start:end]))
        if shuffled_support_fraction >= observed_support_fraction:
            hits += 1

    pvalue = (hits + 1) / (iterations + 1)
    return hits, float(pvalue)


def annotate_interactions_with_monte_carlo(
    interaction_frame: pd.DataFrame,
    structure_profiles: dict[str, dict[str, np.ndarray | None]],
    coordinate_system: str,
    iterations: int,
    alpha: float,
    disorder_threshold: float,
    binding_threshold: float,
    min_support_fraction: float,
    require_binding: bool,
    seed: int,
) -> pd.DataFrame:
    """Annotate an interaction table with Monte Carlo motif statistics."""

    annotated_rows: list[dict[str, Any]] = []
    rng = np.random.default_rng(seed)

    for _, row in interaction_frame.iterrows():
        annotated_row = row.to_dict()
        protein_name = str(row['human_protein'])

        if protein_name not in structure_profiles:
            annotated_row.update(
                {
                    'filter_status': 'missing_structure',
                    'sequence_length': math.nan,
                    'motif_length': math.nan,
                    'supportive_residue_count': math.nan,
                    'supportive_residue_fraction': math.nan,
                    'observed_support_fraction': math.nan,
                    'observed_avg_disorder': math.nan,
                    'observed_avg_binding': math.nan,
                    'monte_carlo_hits': math.nan,
                    'monte_carlo_pvalue': math.nan,
                    'passes_monte_carlo': False,
                },
            )
            annotated_rows.append(annotated_row)
            continue

        disorder_profile = structure_profiles[protein_name]['disorder']
        binding_profile = structure_profiles[protein_name]['binding']
        start, end = _convert_coordinates(
            row['start'],
            row['end'],
            coordinate_system = coordinate_system,
        )

        if start < 0 or end <= start or end > len(disorder_profile):
            annotated_row.update(
                {
                    'filter_status': 'invalid_coordinates',
                    'sequence_length': float(len(disorder_profile)),
                    'motif_length': math.nan,
                    'supportive_residue_count': math.nan,
                    'supportive_residue_fraction': math.nan,
                    'observed_support_fraction': math.nan,
                    'observed_avg_disorder': math.nan,
                    'observed_avg_binding': math.nan,
                    'monte_carlo_hits': math.nan,
                    'monte_carlo_pvalue': math.nan,
                    'passes_monte_carlo': False,
                },
            )
            annotated_rows.append(annotated_row)
            continue

        support_mask = _support_mask(
            disorder_profile = disorder_profile,
            binding_profile = binding_profile,
            disorder_threshold = disorder_threshold,
            binding_threshold = binding_threshold,
            require_binding = require_binding,
        )
        segment_lengths = _segment_lengths(support_mask)
        observed_stats = compute_motif_statistics(
            disorder_profile = disorder_profile,
            binding_profile = binding_profile,
            support_mask = support_mask,
            start = start,
            end = end,
        )
        hits, pvalue = monte_carlo_pvalue(
            sequence_length = len(disorder_profile),
            segment_lengths = segment_lengths,
            start = start,
            end = end,
            observed_support_fraction = observed_stats['observed_support_fraction'],
            iterations = iterations,
            rng = rng,
        )
        passes = (
            observed_stats['observed_support_fraction'] >= min_support_fraction
            and pvalue <= alpha
        )

        annotated_row.update(observed_stats)
        annotated_row.update(
            {
                'filter_status': 'ok',
                'monte_carlo_hits': hits,
                'monte_carlo_pvalue': pvalue,
                'passes_monte_carlo': bool(passes),
            },
        )
        annotated_rows.append(annotated_row)

    return pd.DataFrame(annotated_rows)


def run_monte_carlo_filter(
    interaction_file: PathLike,
    output_file: PathLike,
    filtered_output_file: PathLike | None = None,
    structure_scores: PathLike | None = None,
    fasta_file: PathLike | None = None,
    coordinate_system: str = 'zero_based_half_open',
    iterations: int = 1000,
    alpha: float = 0.05,
    disorder_threshold: float = 0.60,
    binding_threshold: float = 0.60,
    min_support_fraction: float = 1.0,
    require_binding: bool = False,
    seed: int = 0,
    force_cpu: bool = False,
    gpu: str | int = 0,
) -> pd.DataFrame:
    """Run the standalone Monte Carlo motif filter."""

    interaction_frame = read_interaction_table(interaction_file)
    structure_profiles = load_structure_profiles(
        interaction_frame = interaction_frame,
        structure_scores = structure_scores,
        fasta_file = fasta_file,
        force_cpu = force_cpu,
        gpu = gpu,
    )
    annotated = annotate_interactions_with_monte_carlo(
        interaction_frame = interaction_frame,
        structure_profiles = structure_profiles,
        coordinate_system = coordinate_system,
        iterations = iterations,
        alpha = alpha,
        disorder_threshold = disorder_threshold,
        binding_threshold = binding_threshold,
        min_support_fraction = min_support_fraction,
        require_binding = require_binding,
        seed = seed,
    )
    annotated.to_csv(output_file, sep = '\t', index = False)

    if filtered_output_file is not None:
        annotated.loc[annotated['passes_monte_carlo']].to_csv(
            filtered_output_file,
            sep = '\t',
            index = False,
        )

    return annotated


def main() -> None:
    """Run the command-line interface."""

    args = parse_args()
    run_monte_carlo_filter(
        interaction_file = args.interaction_file,
        output_file = args.output_file,
        filtered_output_file = args.filtered_output_file,
        structure_scores = args.structure_scores,
        fasta_file = args.fasta_file,
        coordinate_system = args.coordinate_system,
        iterations = args.iterations,
        alpha = args.alpha,
        disorder_threshold = args.disorder_threshold,
        binding_threshold = args.binding_threshold,
        min_support_fraction = args.min_support_fraction,
        require_binding = args.require_binding,
        seed = args.seed,
        force_cpu = args.force_cpu,
        gpu = args.gpu,
    )


if __name__ == '__main__':
    main()
