"""Monte Carlo over-representation filtering of DMIs against a pooled disordered-region composition."""

import functools
import re
import zlib
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from ..utils import dmi_reader, fasta
from . import dmi, idr_filter

MONTE_CARLO_COLUMNS = [
    "monte_carlo_hits",
    "monte_carlo_pvalue",
    "monte_carlo_qvalue",
    "passes_monte_carlo",
]
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"  # canonical order for the composition count vector
AMINO_ACID_ARRAY = np.array(list(AMINO_ACIDS))


class _ScoringParams(NamedTuple):
    """The disorder/sampling parameters shared across the scoring helpers."""

    method: str
    disorder_cutoff: float
    iterations: int
    seed: int
    force_cpu: bool
    gpu_num: int


def _count_matches(regex: str, sequence: str) -> int:
    """Count non-overlapping regex matches in a sequence.

    Uses the same re.finditer matching Module 6 applies, so observed and sampled counts are
    measured identically.

    Args:
        regex: Motif class regex pattern.
        sequence: Amino acid sequence to search.

    Returns:
        The number of non-overlapping matches of regex in sequence.
    """

    return sum(1 for _ in re.finditer(regex, sequence))


def _disordered_residue_counts(sequence: str, params: _ScoringParams) -> np.ndarray:
    """Tally a sequence's disordered residues into an AMINO_ACIDS count vector.

    Profiles the sequence (idr_filter.cached_profile), masks residues at the disorder cutoff, and
    counts the masked residues per canonical amino acid (non-canonical residues are ignored).

    Args:
        sequence: Amino acid sequence to profile and tally.
        params: The shared disorder/sampling parameters.

    Returns:
        An AMINO_ACIDS-length int count vector of the sequence's disordered residues.
    """

    disorder_profile, _ = idr_filter.cached_profile(
        params.method, sequence, params.force_cpu, params.gpu_num
    )
    mask = idr_filter.disordered_mask(disorder_profile, params.disorder_cutoff)

    counts = np.zeros(len(AMINO_ACIDS), dtype=np.int64)
    for index, amino_acid in enumerate(AMINO_ACIDS):
        counts[index] = sum(
            residue == amino_acid
            for residue, is_disordered in zip(sequence, mask)
            if is_disordered
        )

    return counts


def _build_species_pools(
    interaction_table: pd.DataFrame,
    human_lookup: dict[str, str],
    bacterial_lookup: dict[str, str],
    params: _ScoringParams,
) -> dict[str, dict]:
    """Precompute per-species pooled disordered-residue counts and per-protein counts.

    Profiles each unique motif protein referenced by 'forward' (human) and 'reverse' (bacterial)
    rows once, accumulating a per-species global count vector plus a {protein_id: count_vector} map
    for leave-one-out.

    Args:
        interaction_table: Module 7's IDR-filtered table.
        human_lookup: UniProt accession -> sequence for human proteins.
        bacterial_lookup: UniProt accession -> sequence for bacterial proteins.
        params: The shared disorder/sampling parameters.

    Returns:
        A dict keyed by dmi_type ('forward'/'reverse'), each value a dict with 'global' (the pooled
        count vector) and 'per_protein' ({protein_id: count_vector}).
    """

    lookups = {"forward": human_lookup, "reverse": bacterial_lookup}
    pools: dict[str, dict] = {direction: {"per_protein": {}} for direction in lookups}

    for direction, lookup in lookups.items():
        directional = interaction_table[interaction_table["dmi_type"] == direction]
        protein_ids = {
            dmi_reader.motif_side(row)[0] for _, row in directional.iterrows()
        }

        for protein_id in protein_ids:
            sequence = lookup.get(protein_id)
            if sequence is not None:
                pools[direction]["per_protein"][protein_id] = (
                    _disordered_residue_counts(sequence, params)
                )

    for direction in lookups:
        per_protein = pools[direction]["per_protein"].values()
        pools[direction]["global"] = (
            sum(per_protein)
            if per_protein
            else np.zeros(len(AMINO_ACIDS), dtype=np.int64)
        )

    return pools


def _stable_hash(text: str) -> int:
    """Return a deterministic, run-stable integer hash of text (zlib.crc32), for RNG seeding.

    Args:
        text: String to hash.

    Returns:
        The CRC-32 of text's UTF-8 bytes.
    """

    return zlib.crc32(text.encode("utf-8"))


@functools.lru_cache(maxsize=None)
def _cached_null_counts(
    regex: str,
    region_length: int,
    loo_counts: tuple[int, ...],
    iterations: int,
    seed: int,
) -> tuple[int, ...]:
    """Null regex match-count distribution for one (regex, region-length, composition) combination.

    Draws `iterations` synthetic length-`region_length` sequences i.i.d. from the frequencies
    implied by loo_counts (over AMINO_ACIDS) and counts regex matches in each. The RNG is derived
    from seed and stable hashes of regex, loo_counts, and region_length, so the draw is reproducible
    and independent per combination; Module 6's fan-out reuses the cached result.

    Args:
        regex: Motif class regex pattern.
        region_length: Length D of the disordered region being matched against.
        loo_counts: Leave-one-out AMINO_ACIDS count vector (the pooled background composition).
        iterations: Number of synthetic sequences to sample.
        seed: Base random seed.

    Returns:
        The per-sample regex match counts, one int per iteration.
    """

    total = sum(loo_counts)
    if total == 0:
        frequencies = np.full(len(AMINO_ACIDS), 1 / len(AMINO_ACIDS))
    else:
        frequencies = np.array(loo_counts, dtype=np.float64) / total

    rng = np.random.default_rng(
        [seed, _stable_hash(regex), _stable_hash(str(loo_counts)), region_length]
    )
    samples = rng.choice(
        AMINO_ACID_ARRAY, size=(iterations, region_length), p=frequencies
    )

    return tuple(_count_matches(regex, "".join(sample)) for sample in samples)


def _locate_region(
    row: pd.Series,
    sequence_lookup: dict[str, str],
    params: _ScoringParams,
) -> tuple[str, int, int, str] | None:
    """Locate the disordered region containing a row's motif start.

    Args:
        row: One row of Module 7's IDR-filtered table.
        sequence_lookup: UniProt accession -> sequence for the row's motif-bearing species.
        params: The shared disorder/sampling parameters.

    Returns:
        (protein_id, region_start, region_end, region_sequence) for the maximal disordered run
        containing the motif start, or None if the sequence is missing or the motif is not in a
        disordered region.
    """

    protein_id, start, _ = dmi_reader.motif_side(row)
    sequence = sequence_lookup.get(protein_id)
    if sequence is None:
        return None

    disorder_profile, _ = idr_filter.cached_profile(
        params.method, sequence, params.force_cpu, params.gpu_num
    )
    mask = idr_filter.disordered_mask(disorder_profile, params.disorder_cutoff)
    region = idr_filter.disordered_region(mask, start)
    if region is None:
        return None

    region_start, region_end = region
    return protein_id, region_start, region_end, sequence[region_start:region_end]


def _row_instances(
    row: pd.Series,
    sequence_lookup: dict[str, str],
    species_pool: dict,
    regexes: dict[str, str],
    params: _ScoringParams,
) -> list[tuple[tuple, int, float]] | None:
    """Resolve one row into its testable motif instances and their raw MC p-values.

    Locates the disordered region, forms the leave-one-out background, and scores each resolvable
    class regex (observed region count vs. the cached null) into a raw p-value.

    Args:
        row: One row of Module 7's IDR-filtered table.
        sequence_lookup: UniProt accession -> sequence for the row's motif-bearing species.
        species_pool: The 'global'/'per_protein' pool for the row's species.
        regexes: motif_id -> regex mapping.
        params: The shared disorder/sampling parameters.

    Returns:
        One (instance_key, hits, pvalue) tuple per resolvable class, instance_key being
        (regex, protein_id, region_start, region_end). None if the row is not testable.
    """

    located = _locate_region(row, sequence_lookup, params)
    if located is None:
        return None

    protein_id, region_start, region_end, region_seq = located
    resolved = [
        regexes[motif_id]
        for motif_id in dmi_reader.motif_class_ids(row)
        if motif_id in regexes
    ]
    if not resolved:
        return None

    loo_counts = tuple(
        int(count)
        for count in species_pool["global"] - species_pool["per_protein"][protein_id]
    )
    bounds = (protein_id, region_start, region_end)

    return [
        _score_class(regex, region_seq, bounds, loo_counts, params)
        for regex in resolved
    ]


def _score_class(
    regex: str,
    region_seq: str,
    bounds: tuple[str, int, int],
    loo_counts: tuple[int, ...],
    params: _ScoringParams,
) -> tuple[tuple, int, float]:
    """Score one class regex against the cached null for its region.

    Args:
        regex: Motif class regex pattern.
        region_seq: The disordered region subsequence.
        bounds: (protein_id, region_start, region_end) identifying the region.
        loo_counts: Leave-one-out AMINO_ACIDS count vector (background composition).
        params: The shared disorder/sampling parameters.

    Returns:
        (instance_key, hits, pvalue): hits is the number of null samples matching at least as often
        as observed; instance_key is (regex, protein_id, region_start, region_end).
    """

    protein_id, region_start, region_end = bounds
    observed = _count_matches(regex, region_seq)
    null = _cached_null_counts(
        regex, region_end - region_start, loo_counts, params.iterations, params.seed
    )
    hits = sum(count >= observed for count in null)
    pvalue = (hits + 1) / (params.iterations + 1)

    return (regex, protein_id, region_start, region_end), hits, pvalue


def _score_all_instances(
    interaction_table: pd.DataFrame,
    human_lookup: dict[str, str],
    bacterial_lookup: dict[str, str],
    params: _ScoringParams,
) -> tuple[dict[tuple, tuple[int, float]], list[list[tuple] | None]]:
    """Phase 1: build the species pools, then score every row's instances, de-duplicating them.

    Picks the human vs bacterial lookup/pool by dmi_type and calls _row_instances; fan-out
    duplicates of one instance collapse to a single scored entry.

    Args:
        interaction_table: Module 7's IDR-filtered table.
        human_lookup, bacterial_lookup: Per-species accession -> sequence maps.
        params: The shared disorder/sampling parameters.

    Returns:
        (instance_stats, row_keys): instance_stats maps each unique instance_key -> (hits, pvalue),
        computed once; row_keys is aligned to the table, each element the list of instance_keys for
        that row (or None if the row was not testable).
    """

    regexes = dmi.load_motif_regexes()
    pools = _build_species_pools(
        interaction_table, human_lookup, bacterial_lookup, params
    )

    instance_stats: dict[tuple, tuple[int, float]] = {}
    row_keys: list[list[tuple] | None] = []

    for _, row in interaction_table.iterrows():
        if row["dmi_type"] == "forward":
            lookup, pool = human_lookup, pools["forward"]
        else:
            lookup, pool = bacterial_lookup, pools["reverse"]

        instances = _row_instances(row, lookup, pool, regexes, params)
        if instances is None:
            row_keys.append(None)
            continue

        for key, hits, pvalue in instances:
            instance_stats.setdefault(key, (hits, pvalue))
        row_keys.append([key for key, _, _ in instances])

    return instance_stats, row_keys


def _benjamini_hochberg(
    instance_stats: dict[tuple, tuple[int, float]],
) -> dict[tuple, float]:
    """Phase 2: BH-adjust the unique-instance p-values into q-values.

    Runs scipy.stats.false_discovery_control(pvals, method='bh') over the unique-instance p-values
    (one entry per distinct test, so m is the count of distinct instances).

    Args:
        instance_stats: unique instance_key -> (hits, pvalue) from _score_all_instances.

    Returns:
        A dict mapping each instance_key to its BH-adjusted q-value. Empty if there are no
        instances.
    """

    if not instance_stats:
        return {}

    keys = list(instance_stats)
    pvals = [instance_stats[key][1] for key in keys]
    qvals = false_discovery_control(pvals, method="bh")

    return dict(zip(keys, qvals))


def _filter_rows(
    interaction_table: pd.DataFrame,
    row_keys: list[list[tuple] | None],
    instance_stats: dict[tuple, tuple[int, float]],
    instance_qvalues: dict[tuple, float],
    alpha: float,
) -> list[dict]:
    """Phase 3: broadcast q-values back to rows and keep the passing ones.

    For each testable row, selects the class instance with the smallest q-value (== smallest
    p-value, BH being order-preserving) and emits an annotated record iff that q-value <= alpha.

    Args:
        interaction_table: Module 7's IDR-filtered table.
        row_keys: Per-row instance_key lists (None for untestable rows).
        instance_stats: unique instance_key -> (hits, pvalue).
        instance_qvalues: unique instance_key -> q-value.
        alpha: Target false discovery rate; a row passes iff its best q-value <= alpha.

    Returns:
        A list of passing-row dicts, each row.to_dict() updated with the four Monte Carlo columns.
    """

    records = []
    for (_, row), keys in zip(interaction_table.iterrows(), row_keys):
        if keys is None:
            continue

        best_key = min(keys, key=lambda key: instance_qvalues[key])
        hits, pvalue = instance_stats[best_key]
        qvalue = instance_qvalues[best_key]
        if qvalue > alpha:
            continue

        records.append(
            {
                **row.to_dict(),
                "monte_carlo_hits": hits,
                "monte_carlo_pvalue": pvalue,
                "monte_carlo_qvalue": qvalue,
                "passes_monte_carlo": True,
            }
        )

    return records


def _run_filter(
    interaction_table: pd.DataFrame,
    human_sequences: dict[str, str] | None,
    bacterial_sequences: dict[str, str] | None,
    params: _ScoringParams,
    alpha: float,
) -> pd.DataFrame:
    """Run the three-phase scoring after inputs have been validated.

    Args:
        interaction_table: Module 7's IDR-filtered table.
        human_sequences: Human FASTA header -> sequence, or None.
        bacterial_sequences: Bacterial FASTA header -> sequence, or None.
        params: The shared disorder/sampling parameters.
        alpha: Target false discovery rate.

    Returns:
        interaction_table restricted to passing rows, with the four MONTE_CARLO_COLUMNS appended.
    """

    human_lookup = fasta.build_sequence_lookup(human_sequences)
    bacterial_lookup = fasta.build_sequence_lookup(bacterial_sequences)

    instance_stats, row_keys = _score_all_instances(
        interaction_table, human_lookup, bacterial_lookup, params
    )
    instance_qvalues = _benjamini_hochberg(instance_stats)
    records = _filter_rows(
        interaction_table, row_keys, instance_stats, instance_qvalues, alpha
    )

    output_columns = [*interaction_table.columns, *MONTE_CARLO_COLUMNS]
    return pd.DataFrame(records, columns=output_columns)


def filter_by_monte_carlo(
    interaction_table: pd.DataFrame,
    human_sequences: dict[str, str] | None,
    bacterial_sequences: dict[str, str] | None,
    method: str,
    disorder_cutoff: float,
    iterations: int = 1000,
    alpha: float = 0.05,
    seed: int = 0,
    force_cpu: bool = False,
    gpu_num: int = 0,
) -> pd.DataFrame:
    """Filter Module 7's DMI table by Monte Carlo over-representation against pooled disorder.

    Each unique motif instance (regex, protein, disordered region) is counted in its region and
    re-counted over `iterations` length-matched samples from the species' pooled disordered
    composition (leave-one-out); raw p-values are Benjamini-Hochberg corrected across all instances
    and a row passes iff its best class's q-value <= alpha (the target false discovery rate).

    Args:
        interaction_table: Module 7's IDR-filtered table (extra columns are passed through).
        human_sequences: Human FASTA header -> sequence (required for 'forward' rows).
        bacterial_sequences: Bacterial FASTA header -> sequence (required for 'reverse' rows).
        method: 'iupred' or 'aiupred' disorder predictor.
        disorder_cutoff: Per-residue score at/above which a residue is disordered.
        iterations: Samples per instance; the p-value floor is 1/(iterations + 1), so raise it with
            the number of instances tested (see the module docstring).
        alpha: Target false discovery rate; a row passes iff its BH q-value <= alpha.
        seed: Random seed for reproducible sampling.
        force_cpu: Force CPU inference (aiupred only; ignored for iupred).
        gpu_num: GPU index (aiupred only; ignored for iupred).

    Returns:
        interaction_table restricted to passing rows, with the four MONTE_CARLO_COLUMNS appended.

    Raises:
        ValueError: If method is invalid, or a present direction's sequences are missing.
    """

    idr_filter.validate_method(method)
    dmi_reader.validate_direction_sequences(
        interaction_table, human_sequences, bacterial_sequences
    )
    params = _ScoringParams(
        method, disorder_cutoff, iterations, seed, force_cpu, gpu_num
    )

    return _run_filter(
        interaction_table, human_sequences, bacterial_sequences, params, alpha
    )
