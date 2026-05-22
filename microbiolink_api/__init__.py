#!/usr/bin/env python

"""User-friendly Python API for the MicrobioLink workflow."""

from __future__ import annotations

import importlib.metadata

from microbiolink_api.dmi import DMIResourceBundle
from microbiolink_api.dmi import DomainMotifInteraction
from microbiolink_api.dmi import extract_uniprot_id
from microbiolink_api.dmi import interactions_to_dataframe
from microbiolink_api.dmi import load_default_dmi_resource_bundle
from microbiolink_api.dmi import predict_domain_motif_interactions
from microbiolink_api.dmi import predict_domain_motif_interactions_from_data
from microbiolink_api.dmi import read_bacterial_domain_table
from microbiolink_api.dmi import read_elm_regex_table
from microbiolink_api.dmi import read_fasta_sequences
from microbiolink_api.dmi import read_motif_domain_table
from microbiolink_api.dmi import select_sequences_by_uniprot_ids
from microbiolink_api.dmi import write_domain_motif_interactions
from microbiolink_api.dmi import write_fasta_sequences
from microbiolink_api.exceptions import InputFormatError
from microbiolink_api.exceptions import MicrobioLinkAPIError
from microbiolink_api.expression import filter_count_matrix_file
from microbiolink_api.expression import filter_counts_by_zscore
from microbiolink_api.expression import read_count_matrix
from microbiolink_api.microbiome import bacterial_domain_dataframe_to_mapping
from microbiolink_api.microbiome import fetch_bacterial_domain_table_from_file
from microbiolink_api.microbiome import fetch_bacterial_domain_table_from_ids
from microbiolink_api.microbiome import filter_bacterial_domain_table_by_location
from microbiolink_api.microbiome import read_microbiome_identifiers
from microbiolink_api.workflows import DMIWorkflowResult
from microbiolink_api.workflows import run_dmi_workflow


try:
    __version__ = importlib.metadata.version('microbiolink')
except importlib.metadata.PackageNotFoundError:
    __version__ = '0.0.1'


__all__ = [
    '__version__',
    'DMIWorkflowResult',
    'DMIResourceBundle',
    'DomainMotifInteraction',
    'InputFormatError',
    'MicrobioLinkAPIError',
    'bacterial_domain_dataframe_to_mapping',
    'extract_uniprot_id',
    'fetch_bacterial_domain_table_from_file',
    'fetch_bacterial_domain_table_from_ids',
    'filter_bacterial_domain_table_by_location',
    'filter_count_matrix_file',
    'filter_counts_by_zscore',
    'interactions_to_dataframe',
    'load_default_dmi_resource_bundle',
    'predict_domain_motif_interactions',
    'predict_domain_motif_interactions_from_data',
    'read_bacterial_domain_table',
    'read_count_matrix',
    'read_elm_regex_table',
    'read_fasta_sequences',
    'read_microbiome_identifiers',
    'read_motif_domain_table',
    'run_dmi_workflow',
    'select_sequences_by_uniprot_ids',
    'write_domain_motif_interactions',
    'write_fasta_sequences',
]
