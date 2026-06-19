# Results Guide

## Recommended Reading Order

For a quick read:

1. `../05_results/benchmark_overview.tsv`
2. `../05_results/approach_summary.tsv`
3. `../05_results/approach_union_summary.tsv`
4. `../05_results/phisto_unique_pairs.tsv`
5. One or more of the detailed TP tables:
   `forward_dmi_true_positive_details.tsv`,
   `reverse_dmi_true_positive_details.tsv`,
   `ddi_true_positive_details.tsv`

## Input and Annotation Tables

### `../03_resolved_annotations/bacterial_uniprot_annotations.tsv`

One row per bacterial accession from the PHISTO benchmark panel.

Columns:

- `requested_accession`: accession present in the PHISTO export.
- `current_accession`: accession returned by UniProt for that query.
- `entry_name`: UniProt entry name used in FASTA headers.
- `pfam_count`: number of Pfam domains resolved for that protein.
- `pfams`: semicolon-separated Pfam list.
- `sequence_length`: sequence length in amino acids.
- `resolution_mode`: identifier resolution path used by the script.

### `../03_resolved_annotations/human_uniprot_annotations.tsv`

Same schema as the bacterial annotation table, but for human proteins.

### `../04_benchmark_inputs/bacterial_sequences.fasta`

FASTA file containing bacterial proteins with resolved sequence.

### `../04_benchmark_inputs/human_sequences.fasta`

FASTA file containing human proteins with resolved sequence.

### `../04_benchmark_inputs/bacterial_domains.tsv`

Two-column protein-to-Pfam table used as the bacterial domain input for forward
DMI and as the bacterial domain source for DDI.

Columns:

- `protein`
- `pfam`

The `pfam` field is semicolon-separated when multiple domains are present.

### `../04_benchmark_inputs/human_domains.tsv`

Same schema as `bacterial_domains.tsv`, but for the host side.

## Benchmark Definition Tables

### `../05_results/phisto_unique_pairs.tsv`

This is the deduplicated true-positive benchmark panel derived from the raw
PHISTO CSV. One row equals one unique pathogen-human protein pair.

Columns:

- `pathogen_accession`
- `human_accession`
- `pathogen_name`
- `pathogen_protein_name`
- `human_protein_name`
- `supporting_rows`: number of raw PHISTO rows collapsing into this pair.
- `taxonomy_ids`: semicolon-separated taxonomy identifiers observed for the pair.
- `methods`: semicolon-separated PHISTO experimental methods observed for the
  pair.
- `pubmed_ids`: semicolon-separated supporting PMIDs observed for the pair.

### `../05_results/phisto_method_counts.tsv`

Frequency table of PHISTO experimental methods across the raw PHISTO export.

Columns:

- `experimental_method`
- `count`

### `../05_results/phisto_pmid_counts.tsv`

Frequency table of PMIDs across the raw PHISTO export.

Columns:

- `pubmed_id`
- `count`

## Summary Tables

### `../05_results/benchmark_overview.tsv`

Run-level overview and panel composition.

Columns:

- `metric`
- `value`

Important rows:

- `phisto_rows`
- `phisto_unique_pairs`
- `resolved_bacterial_annotations`
- `resolved_human_annotations`
- `bacterial_with_pfam`
- `human_with_pfam`
- `bacterial_with_sequence`
- `human_with_sequence`

### `../05_results/approach_summary.tsv`

Primary benchmark summary. One row per benchmarked method-resource combination.

Columns:

- `approach`: one of `forward_dmi`, `reverse_dmi`, or `ddi`.
- `resource_name`: resource bundle name.
- `benchmark_unique_pairs`: total PHISTO pair count, fixed at `9027`.
- `analyzable_true_pairs`: PHISTO pairs for which the required structural input
  existed for that method.
- `resource_rows`: number of resource interaction rules for DDI rows. This is
  empty for DMI rows because DMI uses merged motif-domain rule collections.
- `unique_predicted_pairs`: number of unique benchmark-panel protein pairs
  predicted by that method across the tested protein cross-product.
- `raw_prediction_rows`: raw number of supporting hit rows before pair
  deduplication.
- `true_positives_recovered`: unique PHISTO benchmark pairs recovered.
- `false_positive_pairs_within_panel_cross_product`: predicted benchmark-panel
  pairs that are not in the PHISTO TP set.
- `recall_total`: recovered TPs divided by all `9027` PHISTO pairs.
- `recall_on_analyzable_pairs`: recovered TPs divided by the subset with the
  necessary structural inputs for that method.
- `precision_unique_pairs`: recovered TP pairs divided by `unique_predicted_pairs`.

Interpretation note:

- `recall_total` answers "how much of PHISTO is covered overall?"
- `recall_on_analyzable_pairs` answers "how much is covered once missing Pfam or
  sequence annotations are excluded?"

### `../05_results/approach_union_summary.tsv`

Compact summary of unique TP recovery for the three most relevant aggregate
sets:

- `forward_dmi`
- `reverse_dmi`
- `ddi_3did_plus_domine_v2_all`
- `union_all_three`

Columns:

- `set_name`
- `unique_true_positive_pairs`

## Detailed True-Positive Tables

These tables contain only benchmark pairs that were actually recovered. They do
not include the large number of predicted non-TP pairs captured only in the
summary metrics.

### `../05_results/forward_dmi_true_positive_details.tsv`

Recovered PHISTO pairs supported by forward DMI logic.

Columns:

- `direction`: always `forward`
- `bacterial_accession`
- `human_accession`
- `motif`: matched host-side motif identifier
- `domain`: matched bacterial Pfam domain
- `resource`: motif source, `ELM` or `3did`
- `motif_match_count`: number of matches of that motif in the host protein
- `first_start`: first match start coordinate, zero-based
- `first_end`: first match end coordinate, zero-based exclusive

### `../05_results/reverse_dmi_true_positive_details.tsv`

Recovered PHISTO pairs supported by reverse DMI logic, where the motif is on
the bacterial protein and the matching domain is on the host protein.

Schema is the same as the forward DMI TP table.

### `../05_results/ddi_true_positive_details.tsv`

Recovered PHISTO pairs supported by DDI logic.

Columns:

- `resource_name`: DDI resource set that recovered the pair
- `bacterial_accession`
- `human_accession`
- `bacterial_domain`
- `human_domain`

Interpretation note:

- The same TP pair may appear more than once if it is recovered by multiple
  DDI resource sets.

## Resource Tables

### `../02_reference_resources/dmi/`

Contains the exact DMI rule tables used by the benchmark:

- packaged ELM motif classes and motif-domain mappings
- packaged 3did-derived structural DMI motif classes and motif-domain mappings

### `../02_reference_resources/ddi/`

Contains the exact DDI rule tables used or derived for the benchmark:

- `pfam_interactions_3did_current.tsv`: current local 3did Pfam-Pfam table
- `domine_v2_all_pfam_pairs.tsv`: all DOMINE v2 Pfam-Pfam pairs derived from
  `INTERACTION.txt`
- `domine_v2_hc_pfam_pairs.tsv`: DOMINE v2 high-confidence subset derived from
  `INTERACTION.txt`
- `domine_hc_legacy_pfam_pairs.tsv`: local legacy MicrobioLink DOMINE-style
  HC baseline previously used in the repository

## Metadata Files

### `../07_metadata/run_command.tsv`

Exact command line used for the benchmark run.

### `../07_metadata/software_environment.tsv`

Minimal run environment metadata.

### `../07_metadata/file_list.txt`

Flat list of all files in the documented bundle.

### `../07_metadata/line_counts.txt`

Line counts for bundle files. Useful for quick integrity checks.

## Interactive Rerun Assets

### `../06_scripts/package_benchmark_workflow.py`

Package-driven support module for the notebook workflow.

Key roles:

- rebuild the PHISTO unique-pair panel from the raw export,
- export accession lists,
- optionally redownload FASTA and Pfam annotations through the packaged generic
  UniProt helpers,
- rerun forward DMI, reverse DMI, and DDI with batched package calls,
- write regenerated summary tables and true-positive detail tables into
  `../10_reproduced_run/`.

Important note:

- This script is not the original provenance snapshot for the documented
  outputs in `../05_results/`. It is the package-oriented rerun layer added for
  the `MicrobioLink-2.1-beta` branch.

### `../09_notebooks/phisto_extended_package_walkthrough.ipynb`

Interactive notebook for stepping through the benchmark with the packaged
MicrobioLink 2.1 API.

The notebook covers:

1. rebuilding the PHISTO panel from the raw export,
2. exporting accession lists,
3. optionally redownloading FASTA and Pfam tables with the new generic package
   helpers,
4. previewing forward DMI, reverse DMI, and DDI schemas,
5. rerunning the benchmark through the package into
   `../10_reproduced_run/05_results/`.

### `../10_reproduced_run/`

Reserved output location for package-based reruns.

Why this exists:

- It prevents the interactive notebook from overwriting the canonical
  documented outputs in `../05_results/`.
- It keeps notebook-generated accession lists, raw UniProt downloads,
  regenerated benchmark inputs, results, and metadata in one place.
