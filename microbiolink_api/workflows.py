#!/usr/bin/env python

"""High-level workflows for the library-style MicrobioLink API."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from microbiolink_api.dmi import DMIResourceBundle
from microbiolink_api.dmi import DomainMotifInteraction
from microbiolink_api.dmi import interactions_to_dataframe
from microbiolink_api.dmi import load_default_dmi_resource_bundle
from microbiolink_api.dmi import predict_domain_motif_interactions_from_data
from microbiolink_api.dmi import read_bacterial_domain_table
from microbiolink_api.dmi import read_elm_regex_table
from microbiolink_api.dmi import read_fasta_sequences
from microbiolink_api.dmi import read_motif_domain_table
from microbiolink_api.dmi import select_sequences_by_uniprot_ids
from microbiolink_api.dmi import write_domain_motif_interactions
from microbiolink_api.dmi import write_fasta_sequences
from microbiolink_api.expression import filter_count_matrix_file
from microbiolink_api.exceptions import InputFormatError
from microbiolink_api.microbiome import bacterial_domain_dataframe_to_mapping
from microbiolink_api.microbiome import fetch_bacterial_domain_table_from_file
from microbiolink_api.microbiome import filter_bacterial_domain_table_by_location


PathLike = str | Path


@dataclass
class DMIWorkflowResult:
    """Collect the outputs of the expression-to-DMI workflow."""

    filtered_counts: pd.DataFrame
    bacterial_domain_table: pd.DataFrame
    selected_sequences: dict[str, str]
    interactions: list[DomainMotifInteraction]
    interaction_frame: pd.DataFrame


def run_dmi_workflow(
    count_matrix_file: PathLike,
    human_fasta_file: PathLike,
    bacterial_domain_file: PathLike | None = None,
    bacterial_id_file: PathLike | None = None,
    bacterial_id_type: str | None = None,
    bacterial_id_separator: str = '\t',
    bacterial_id_column: int = 1,
    bacterial_location_filters: list[str] | None = None,
    subset_human_fasta_by_expression: bool = False,
    elm_regex_file: PathLike | None = None,
    motif_domain_file: PathLike | None = None,
    zscore_threshold: float = -3,
    resource_bundle: DMIResourceBundle | None = None,
    bacterial_domain_output: PathLike | None = None,
    filtered_counts_output: PathLike | None = None,
    selected_fasta_output: PathLike | None = None,
    dmi_output: PathLike | None = None,
) -> DMIWorkflowResult:
    """Run the user-facing gene-count-to-DMI workflow.

    This workflow filters a host expression matrix, selects the subset of human
    proteins that remain expressed after filtering, and predicts domain-motif
    interactions for that expressed host protein set.

    Args:
        count_matrix_file: Path to the host gene or protein count matrix.
        human_fasta_file: Path to the human protein FASTA file.
        bacterial_domain_file: Optional precomputed bacterial domain table.
        bacterial_id_file: Optional microbiome-side identifier file containing
            UniProt protein IDs or proteome IDs.
        bacterial_id_type: Identifier type for `bacterial_id_file`. Use
            `'Uniprot'` for proteins or `'UP'` for proteomes.
        bacterial_id_separator: Field separator for `bacterial_id_file`.
        bacterial_id_column: One-based identifier column in
            `bacterial_id_file`.
        bacterial_location_filters: Optional case-insensitive location filters
            applied to UniProt subcellular-location annotations for bacterial
            proteins.
        subset_human_fasta_by_expression: Whether to subset the provided human
            FASTA by the filtered count-matrix row identifiers. Leave this as
            `False` when the count matrix uses gene symbols but the FASTA uses
            UniProt headers.
        elm_regex_file: Optional override for the ELM regex table.
        motif_domain_file: Optional override for the motif-domain table.
        zscore_threshold: Minimum z-score to retain a value.
        resource_bundle: Optional in-memory override for the built-in DMI
            resource tables.
        bacterial_domain_output: Optional output path for a fetched bacterial
            domain table.
        filtered_counts_output: Optional path for the filtered matrix.
        selected_fasta_output: Optional path for the filtered FASTA subset.
        dmi_output: Optional path for the DMI result table.

    Returns:
        A structured workflow result with both in-memory tables and records.
    """

    if bacterial_domain_file is None and bacterial_id_file is None:
        raise InputFormatError(
            'Provide either `bacterial_domain_file` or `bacterial_id_file`.',
        )

    if bacterial_id_file is not None and bacterial_id_type is None:
        raise InputFormatError(
            '`bacterial_id_type` is required when `bacterial_id_file` is used.',
        )

    filtered_counts = filter_count_matrix_file(
        count_matrix_file,
        zscore_threshold = zscore_threshold,
        output_file = filtered_counts_output,
    )

    expressed_proteins = [
        str(protein_id)
        for protein_id, values in filtered_counts.iterrows()
        if values.notna().any()
    ]

    all_sequences = read_fasta_sequences(human_fasta_file)
    if subset_human_fasta_by_expression:
        selected_sequences = select_sequences_by_uniprot_ids(
            all_sequences,
            expressed_proteins,
        )
    else:
        selected_sequences = all_sequences

    if selected_fasta_output is not None:
        write_fasta_sequences(selected_sequences, selected_fasta_output)

    resolved_bundle = resource_bundle or load_default_dmi_resource_bundle()
    if bacterial_domain_file is not None:
        bacterial_domain_table = pd.read_csv(
            bacterial_domain_file,
            sep = '\t',
        )
        bacterial_domains = read_bacterial_domain_table(bacterial_domain_file)
    else:
        bacterial_domain_table = fetch_bacterial_domain_table_from_file(
            bacterial_id_file,
            id_type = bacterial_id_type,
            separator = bacterial_id_separator,
            id_column = bacterial_id_column,
            include_location = bool(bacterial_location_filters),
            output_file = bacterial_domain_output,
        )
        if bacterial_location_filters:
            bacterial_domain_table = filter_bacterial_domain_table_by_location(
                bacterial_domain_table,
                location_filters = bacterial_location_filters,
            )
            if bacterial_domain_output is not None:
                bacterial_domain_table.to_csv(
                    bacterial_domain_output,
                    sep = '\t',
                    index = False,
                )
        bacterial_domains = bacterial_domain_dataframe_to_mapping(
            bacterial_domain_table,
        )

    interactions = predict_domain_motif_interactions_from_data(
        human_sequences = selected_sequences,
        elm_regex = (
            read_elm_regex_table(elm_regex_file)
            if elm_regex_file is not None
            else resolved_bundle.elm_regex
        ),
        motif_domains = (
            read_motif_domain_table(motif_domain_file)
            if motif_domain_file is not None
            else resolved_bundle.motif_domains
        ),
        bacterial_domains = bacterial_domains,
    )
    interaction_frame = interactions_to_dataframe(interactions)

    if dmi_output is not None:
        write_domain_motif_interactions(interactions, dmi_output)

    return DMIWorkflowResult(
        filtered_counts = filtered_counts,
        bacterial_domain_table = bacterial_domain_table,
        selected_sequences = selected_sequences,
        interactions = interactions,
        interaction_frame = interaction_frame,
    )
