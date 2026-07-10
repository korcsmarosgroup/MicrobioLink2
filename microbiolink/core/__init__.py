#!/usr/bin/env python

"""User-friendly Python API for the MicrobioLink workflow."""

from __future__ import annotations

import importlib.metadata

from microbiolink.core.dmi import BidirectionalDomainMotifInteraction
from microbiolink.core.dmi import DMIResourceBundle
from microbiolink.core.dmi import DomainMotifInteraction
from microbiolink.core.dmi import bidirectional_interactions_to_dataframe
from microbiolink.core.dmi import extract_uniprot_id
from microbiolink.core.dmi import interactions_to_dataframe
from microbiolink.core.dmi import load_default_3did_dmi_resource_bundle
from microbiolink.core.dmi import load_default_dmi_resource_bundle
from microbiolink.core.dmi import load_default_elm_dmi_resource_bundle
from microbiolink.core.dmi import merge_dmi_resource_bundles
from microbiolink.core.dmi import predict_bidirectional_domain_motif_interactions
from microbiolink.core.dmi import predict_bidirectional_domain_motif_interactions_from_data
from microbiolink.core.dmi import predict_domain_motif_interactions
from microbiolink.core.dmi import predict_domain_motif_interactions_from_data
from microbiolink.core.dmi import predict_reverse_domain_motif_interactions
from microbiolink.core.dmi import predict_reverse_domain_motif_interactions_from_data
from microbiolink.core.dmi import read_bacterial_domain_table
from microbiolink.core.dmi import read_elm_regex_table
from microbiolink.core.dmi import read_fasta_sequences
from microbiolink.core.dmi import read_motif_domain_table
from microbiolink.core.dmi import read_protein_domain_table
from microbiolink.core.dmi import resolve_dmi_resource_bundle_by_name
from microbiolink.core.dmi import select_sequences_by_uniprot_ids
from microbiolink.core.dmi import write_bidirectional_domain_motif_interactions
from microbiolink.core.dmi import write_domain_motif_interactions
from microbiolink.core.dmi import write_fasta_sequences
from microbiolink.core.ddi import DDIResourceBundle
from microbiolink.core.ddi import DomainDomainInteraction
from microbiolink.core.ddi import ddi_interactions_to_dataframe
from microbiolink.core.ddi import load_default_3did_ddi_resource_bundle
from microbiolink.core.ddi import load_default_ddi_resource_bundle
from microbiolink.core.ddi import load_default_domine_all_ddi_resource_bundle
from microbiolink.core.ddi import load_default_domine_hc_ddi_resource_bundle
from microbiolink.core.ddi import merge_ddi_resource_bundles
from microbiolink.core.ddi import predict_domain_domain_interactions
from microbiolink.core.ddi import predict_domain_domain_interactions_from_data
from microbiolink.core.ddi import write_domain_domain_interactions
from microbiolink.core.exceptions import InputFormatError
from microbiolink.core.exceptions import MicrobioLinkError
from microbiolink.core.expression import filter_count_matrix_file
from microbiolink.core.expression import filter_counts_by_zscore
from microbiolink.core.expression import read_count_matrix
from microbiolink.core.microbiome import bacterial_domain_dataframe_to_mapping
from microbiolink.core.microbiome import fetch_bacterial_domain_table_from_file
from microbiolink.core.microbiome import fetch_bacterial_domain_table_from_ids
from microbiolink.core.microbiome import filter_bacterial_domain_table_by_location
from microbiolink.core.microbiome import read_microbiome_identifiers
from microbiolink.core.workflows import DMIWorkflowResult
from microbiolink.core.workflows import run_dmi_workflow


try:
    __version__ = importlib.metadata.version('microbiolink')
except importlib.metadata.PackageNotFoundError:
    __version__ = '0.0.2'


__all__ = [
    '__version__',
    'BidirectionalDomainMotifInteraction',
    'DDIResourceBundle',
    'DMIWorkflowResult',
    'DMIResourceBundle',
    'DomainDomainInteraction',
    'DomainMotifInteraction',
    'InputFormatError',
    'MicrobioLinkError',
    'bacterial_domain_dataframe_to_mapping',
    'bidirectional_interactions_to_dataframe',
    'ddi_interactions_to_dataframe',
    'extract_uniprot_id',
    'fetch_bacterial_domain_table_from_file',
    'fetch_bacterial_domain_table_from_ids',
    'filter_bacterial_domain_table_by_location',
    'filter_count_matrix_file',
    'filter_counts_by_zscore',
    'interactions_to_dataframe',
    'load_default_3did_ddi_resource_bundle',
    'load_default_3did_dmi_resource_bundle',
    'load_default_ddi_resource_bundle',
    'load_default_dmi_resource_bundle',
    'load_default_domine_all_ddi_resource_bundle',
    'load_default_domine_hc_ddi_resource_bundle',
    'load_default_elm_dmi_resource_bundle',
    'merge_ddi_resource_bundles',
    'merge_dmi_resource_bundles',
    'predict_bidirectional_domain_motif_interactions',
    'predict_bidirectional_domain_motif_interactions_from_data',
    'predict_domain_domain_interactions',
    'predict_domain_domain_interactions_from_data',
    'predict_domain_motif_interactions',
    'predict_domain_motif_interactions_from_data',
    'predict_reverse_domain_motif_interactions',
    'predict_reverse_domain_motif_interactions_from_data',
    'read_bacterial_domain_table',
    'read_count_matrix',
    'read_elm_regex_table',
    'read_fasta_sequences',
    'read_microbiome_identifiers',
    'read_motif_domain_table',
    'read_protein_domain_table',
    'resolve_dmi_resource_bundle_by_name',
    'run_dmi_workflow',
    'select_sequences_by_uniprot_ids',
    'write_bidirectional_domain_motif_interactions',
    'write_domain_domain_interactions',
    'write_domain_motif_interactions',
    'write_fasta_sequences',
]
