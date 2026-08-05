"""Filtering domain-motif interactions by intrinsic disorder and binding likelihood."""

import functools

import numpy as np
import pandas as pd

from ..utils import dmi_reader, fasta
from . import dmi

OUTPUT_COLUMNS = [
    *dmi.OUTPUT_COLUMNS,
    "disordered_score",
    "binding_score",
    "combined_score",
]


def iupred_profile(sequence: str) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-residue IUPred2 disorder and ANCHOR2 binding scores.

    Args:
        sequence: Amino acid sequence to score.

    Returns:
        A (disorder_scores, binding_scores) pair of NumPy arrays, each one score per residue
        in sequence. binding_scores is computed via ANCHOR2 against the long-mode IUPred2
        disorder profile, not the short-mode profile disorder_scores itself uses.
    """

    from iupred import anchor2, iupred

    disorder_short, _ = iupred(sequence, mode="short")
    disorder_long, _ = iupred(sequence, mode="long")
    binding = anchor2(sequence, disorder_long)

    return np.asarray(disorder_short), np.asarray(binding)


def aiupred_profile(
    sequence: str, force_cpu: bool, gpu_num: int
) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-residue AIUPred disorder and binding scores.

    Args:
        sequence: Amino acid sequence to score.
        force_cpu: Force CPU inference even if a GPU is available.
        gpu_num: Index of the GPU to use when force_cpu is False.

    Returns:
        A (disorder_scores, binding_scores) pair of NumPy arrays, each one score per residue
        in sequence.
    """

    from iupred import aiupred_binding, aiupred_disorder

    disorder = aiupred_disorder(sequence, force_cpu=force_cpu, gpu_num=gpu_num)
    binding = aiupred_binding(sequence, force_cpu=force_cpu, gpu_num=gpu_num)

    return np.asarray(disorder), np.asarray(binding)


@functools.lru_cache(maxsize=None)
def cached_profile(
    method: str,
    sequence: str,
    force_cpu: bool,
    gpu_num: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute and cache a sequence's disorder/binding profile, keyed by method and device.

    Module 6's output cross-joins one motif hit against every partner protein sharing a
    compatible domain, so the same motif protein's sequence is commonly scored by many DMI
    rows. Caching here means each unique (method, sequence, force_cpu, gpu_num) combination is
    only ever scored once, regardless of how many rows reference it.

    Args:
        method: 'iupred' or 'aiupred'.
        sequence: Amino acid sequence to score.
        force_cpu: Force CPU inference (aiupred only; accepted but ignored for iupred).
        gpu_num: GPU index (aiupred only; accepted but ignored for iupred).

    Returns:
        A (disorder_scores, binding_scores) pair of NumPy arrays, each one score per residue
        in sequence.

    Raises:
        ValueError: If method is not 'iupred' or 'aiupred'.
    """

    if method == "iupred":
        return iupred_profile(sequence)

    if method == "aiupred":
        return aiupred_profile(sequence, force_cpu, gpu_num)

    raise ValueError(f"method must be 'iupred' or 'aiupred', got {method!r}")


def validate_method(method: str) -> None:
    """Validate that method is a supported disorder predictor.

    Args:
        method: The disorder predictor name to check.

    Raises:
        ValueError: If method is not 'iupred' or 'aiupred'.
    """

    if method not in {"iupred", "aiupred"}:
        raise ValueError(f"method must be 'iupred' or 'aiupred', got {method!r}")


def disordered_mask(disorder_profile: np.ndarray, cutoff: float) -> np.ndarray:
    """Mark which residues are disordered at a per-residue disorder cutoff.

    Args:
        disorder_profile: Per-residue disorder scores.
        cutoff: Score at or above which a residue counts as disordered.

    Returns:
        A boolean array, True where disorder_profile >= cutoff.
    """

    return disorder_profile >= cutoff


def disordered_region(mask: np.ndarray, position: int) -> tuple[int, int] | None:
    """Find the maximal run of disordered residues containing a position.

    Args:
        mask: Boolean per-residue disorder mask (from disordered_mask).
        position: Residue index the region must contain.

    Returns:
        The half-open (start, end) bounds of the maximal True run containing position, or None
        if mask[position] is False.
    """

    if not mask[position]:
        return None

    start = position
    while start > 0 and mask[start - 1]:
        start -= 1

    end = position + 1
    while end < len(mask) and mask[end]:
        end += 1

    return start, end


def _score_row(
    row: pd.Series,
    sequence_lookup: dict[str, str],
    method: str,
    disorder_cutoff: float,
    binding_cutoff: float,
    force_cpu: bool,
    gpu_num: int,
) -> tuple[bool, float, float, float] | None:
    """Score one DMI row's motif window against the per-residue disorder/binding gate.

    Args:
        row: One row of Module 6's DMI output table.
        sequence_lookup: UniProt accession -> sequence mapping for the motif-bearing species
            (built by _build_sequence_lookup).
        method: 'iupred' or 'aiupred'.
        disorder_cutoff: Minimum per-residue disorder score required across the whole window.
        binding_cutoff: Minimum per-residue binding score required across the whole window.
        force_cpu: Passed through to _cached_profile.
        gpu_num: Passed through to _cached_profile.

    Returns:
        (passes, disordered_score, binding_score, combined_score), where passes is True iff
        every residue in the motif window clears both cutoffs, and the three scores are the
        window's mean disorder, mean binding, and their average. None if the motif protein's
        sequence is not present in sequence_lookup.
    """

    motif_id, start, end = dmi_reader.motif_side(row)
    sequence = sequence_lookup.get(motif_id)
    if sequence is None:
        return None

    disorder_profile, binding_profile = cached_profile(
        method, sequence, force_cpu, gpu_num
    )
    disorder_window = disorder_profile[start:end]
    binding_window = binding_profile[start:end]

    passes = bool(
        (disorder_window >= disorder_cutoff).all()
        and (binding_window >= binding_cutoff).all(),
    )
    disordered_score = float(disorder_window.mean())
    binding_score = float(binding_window.mean())
    combined_score = (disordered_score + binding_score) / 2

    return passes, disordered_score, binding_score, combined_score


def _to_output_row(
    row: pd.Series,
    disordered_score: float,
    binding_score: float,
    combined_score: float,
) -> tuple:
    """Translate one passing DMI row plus its scores into an OUTPUT_COLUMNS-shaped tuple.

    Args:
        row: One row of Module 6's DMI output table.
        disordered_score: The motif window's mean disorder score.
        binding_score: The motif window's mean binding score.
        combined_score: The average of disordered_score and binding_score.

    Returns:
        A tuple matching OUTPUT_COLUMNS.
    """

    return (
        row["dmi_type"],
        row["bacterial_uniprot_id"],
        row["bacterial_annotation"],
        row["human_uniprot_id"],
        row["human_annotation"],
        row["start"],
        row["end"],
        row["resource"],
        disordered_score,
        binding_score,
        combined_score,
    )


def _score_all_rows(
    dmi_table: pd.DataFrame,
    human_lookup: dict[str, str],
    bacterial_lookup: dict[str, str],
    method: str,
    disorder_cutoff: float,
    binding_cutoff: float,
    force_cpu: bool,
    gpu_num: int,
) -> list[tuple]:
    """Score every row and collect the ones that pass the per-residue disorder/binding gate.

    Args:
        dmi_table: Module 6's output table.
        human_lookup: UniProt accession -> sequence mapping for human proteins.
        bacterial_lookup: UniProt accession -> sequence mapping for bacterial proteins.
        method: 'iupred' or 'aiupred'.
        disorder_cutoff: Minimum per-residue disorder score required across the window.
        binding_cutoff: Minimum per-residue binding score required across the window.
        force_cpu, gpu_num: Passed through to _cached_profile.

    Returns:
        A list of OUTPUT_COLUMNS-shaped tuples, one per row that passed the gate.
    """

    kept_rows = []
    for _, row in dmi_table.iterrows():
        sequence_lookup = (
            human_lookup if row["dmi_type"] == "forward" else bacterial_lookup
        )
        scored = _score_row(
            row,
            sequence_lookup,
            method,
            disorder_cutoff,
            binding_cutoff,
            force_cpu,
            gpu_num,
        )

        if scored is None:
            continue

        passes, disordered_score, binding_score, combined_score = scored
        if passes:
            kept_rows.append(
                _to_output_row(row, disordered_score, binding_score, combined_score)
            )

    return kept_rows


def filter_by_disorder(
    dmi_table: pd.DataFrame,
    human_sequences: dict[str, str] | None,
    bacterial_sequences: dict[str, str] | None,
    method: str,
    disorder_cutoff: float,
    binding_cutoff: float,
    force_cpu: bool = False,
    gpu_num: int = 0,
) -> pd.DataFrame:
    """Filter Module 6's DMI table to interactions in a disordered, binding-prone region.

    Args:
        dmi_table: Module 6's output table.
        human_sequences: Human FASTA header -> sequence mapping (for 'forward' rows).
        bacterial_sequences: Bacterial FASTA header -> sequence mapping (for 'reverse' rows).
        method: 'iupred' or 'aiupred'.
        disorder_cutoff: Minimum per-residue disorder score required across the motif window.
        binding_cutoff: Minimum per-residue binding score required across the motif window.
        force_cpu: Force CPU inference (aiupred only; ignored for iupred).
        gpu_num: GPU index to use (aiupred only; ignored for iupred).

    Returns:
        dmi_table restricted to rows passing the gate, with disordered_score, binding_score,
        and combined_score columns added.

    Raises:
        ValueError: If method is invalid, or a present DMI direction lacks its sequences.
    """

    validate_method(method)
    dmi_reader.validate_direction_sequences(
        dmi_table, human_sequences, bacterial_sequences
    )

    human_lookup = fasta.build_sequence_lookup(human_sequences)
    bacterial_lookup = fasta.build_sequence_lookup(bacterial_sequences)

    kept_rows = _score_all_rows(
        dmi_table,
        human_lookup,
        bacterial_lookup,
        method,
        disorder_cutoff,
        binding_cutoff,
        force_cpu,
        gpu_num,
    )

    return pd.DataFrame(kept_rows, columns=OUTPUT_COLUMNS)
