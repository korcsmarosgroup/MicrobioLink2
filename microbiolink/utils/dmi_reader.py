"""Read-side helpers for interpreting a DMI interaction table produced by workflow/dmi.py.

workflow/dmi.py *builds* the table (columns per dmi.OUTPUT_COLUMNS, routed by dmi_type in
dmi._to_output_row); the functions here *read* a finished row of that table back - inverting the same
forward/reverse column routing - so consuming modules (idr_filter, monte_carlo) share one
interpretation of the schema. They live in utils, not in dmi.py, because only consumers call them; the
producer never reads its own output rows.
"""

import pandas as pd


def motif_side(row: pd.Series) -> tuple[str, int, int]:
    """Resolve which protein and position carries the motif for one DMI row.

    Args:
        row: One row of a DMI table produced by workflow/dmi.py.

    Returns:
        (motif_uniprot_id, start, end) for whichever side is the motif-bearing side:
        human_uniprot_id/start/end if row['dmi_type'] == 'forward', otherwise the bacterial
        equivalent.
    """

    if row["dmi_type"] == "forward":
        return row["human_uniprot_id"], row["start"], row["end"]

    return row["bacterial_uniprot_id"], row["start"], row["end"]


def motif_class_ids(row: pd.Series) -> list[str]:
    """Resolve the motif class id(s) carried on one DMI row's motif-bearing side.

    Module 6 merges co-located motif classes into one annotation cell joined by '|', so a single
    row may name several classes.

    Args:
        row: One row of a DMI table produced by workflow/dmi.py.

    Returns:
        The motif class ids on the motif-bearing side (human_annotation if
        row['dmi_type'] == 'forward', otherwise bacterial_annotation), split on '|'.
    """

    if row["dmi_type"] == "forward":
        annotation = row["human_annotation"]
    else:
        annotation = row["bacterial_annotation"]

    return annotation.split("|")


def validate_direction_sequences(
    table: pd.DataFrame,
    human_sequences: dict[str, str] | None,
    bacterial_sequences: dict[str, str] | None,
) -> None:
    """Validate that sequences are supplied for whichever directions the table contains.

    Args:
        table: A DMI table produced by workflow/dmi.py.
        human_sequences: Human FASTA header -> sequence mapping, or None.
        bacterial_sequences: Bacterial FASTA header -> sequence mapping, or None.

    Raises:
        ValueError: If table contains 'forward' rows with human_sequences=None (or 'reverse'
            rows with bacterial_sequences=None).
    """

    has_forward = bool((table["dmi_type"] == "forward").any())
    has_reverse = bool((table["dmi_type"] == "reverse").any())

    if has_forward and human_sequences is None:
        raise ValueError(
            "table contains 'forward' rows but human_sequences was not supplied."
        )

    if has_reverse and bacterial_sequences is None:
        raise ValueError(
            "table contains 'reverse' rows but bacterial_sequences was not supplied."
        )
