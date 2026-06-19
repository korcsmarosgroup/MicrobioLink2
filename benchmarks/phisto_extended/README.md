# PHISTO Extended MicrobioLink Benchmark

## Purpose

This folder is a self-contained, documented bundle of the broad PHISTO
bacteria-human benchmark that was run against three upstream MicrobioLink
interaction approaches with extended resources:

1. forward DMI using packaged `ELM + 3did` motif-domain resources,
2. reverse DMI using packaged `ELM + 3did` motif-domain resources,
3. DDI using `3did`, `DOMINE v2`, and combined `3did + DOMINE v2` Pfam-Pfam
   resources.

The folder contains the raw exported benchmark data, the exact reference
resources used, the resolved UniProt-derived benchmark inputs, the result
tables, the exact benchmark script snapshots, an interactive notebook that
walks through the package-based rerun path, and documentation describing how
to interpret the outputs.

## Start Here

The most useful entry points are:

- `05_results/approach_summary.tsv`
- `05_results/approach_union_summary.tsv`
- `05_results/benchmark_overview.tsv`
- `05_results/phisto_unique_pairs.tsv`
- `00_docs/results_guide.md`
- `09_notebooks/phisto_extended_package_walkthrough.ipynb`
- `09_notebooks/README.md`
- `08_figures/figure_1_main_benchmark_overview.svg`
- `08_figures/figure_2_context_and_resource_effects.svg`
- `08_figures/README.md`

## Folder Layout

```text
phisto_extended/
├── README.md
├── 00_docs/
├── 01_raw_sources/
├── 02_reference_resources/
├── 03_resolved_annotations/
├── 04_benchmark_inputs/
├── 05_results/
├── 06_scripts/
├── 07_metadata/
├── 08_figures/
├── 09_notebooks/
└── 10_reproduced_run/
```

Directory meanings:

- `00_docs`: human-readable documentation, manifest, and resource provenance.
- `01_raw_sources`: the PHISTO and DOMINE raw source files used for the run.
- `02_reference_resources`: the exact DMI and DDI resource tables used by the
  benchmark.
- `03_resolved_annotations`: resolved UniProt accessions, Pfam annotations, and
  sequence metadata for the PHISTO protein panel.
- `04_benchmark_inputs`: the FASTA and protein-domain tables actually used by
  the benchmark.
- `05_results`: summary tables and true-positive detail tables for the three
  approaches.
- `06_scripts`: the exact benchmark script snapshot used to create the outputs.
- `07_metadata`: run command, software metadata, file list, and line counts for
  the documented reference run.
- `08_figures`: publication-style SVG figures and figure-level notes.
  The SVG files are vector graphics suitable for editing in Inkscape or
  Illustrator before manuscript submission.
- `09_notebooks`: interactive package-driven walkthrough of the benchmark.
- `10_reproduced_run`: reserved output area for reruns launched from the
  notebook, kept separate from the documented reference outputs.

## Benchmark Scope

The benchmark source is a PHISTO bacteria-only export:

- `9,333` PHISTO table rows
- `9,027` unique pathogen-human protein pairs
- `2,715` unique bacterial accessions
- `3,736` unique human accessions

UniProt resolution status:

- `2,715 / 2,715` bacterial accessions resolved
- `3,736 / 3,736` human accessions resolved
- `0` unresolved bacterial accessions
- `0` unresolved human accessions

Coverage of usable structural inputs:

- `1,232` bacterial proteins with Pfam annotation
- `1,287` bacterial proteins with sequence
- `3,555` human proteins with Pfam annotation
- `3,675` human proteins with sequence

Important panel caveat:

- `8,556 / 9,027` unique PHISTO pairs, or about `94.78%`, come from PMID
  `20711500`.
- This means the broad PHISTO panel is dominated by one high-throughput
  bacteria-human interaction screen.

## Key Results

Main result table:

- `05_results/approach_summary.tsv`

Headline numbers:

| Approach | Resource set | TP pairs recovered | Recall on full PHISTO set | Recall on analyzable pairs |
| --- | --- | ---: | ---: | ---: |
| Forward DMI | ELM + 3did | 144 | 1.60% | 1.63% |
| Reverse DMI | ELM + 3did | 756 | 8.37% | 17.97% |
| DDI | 3did current | 166 | 1.84% | 4.17% |
| DDI | DOMINE v2 HC | 22 | 0.24% | 0.55% |
| DDI | DOMINE v2 all | 217 | 2.40% | 5.45% |
| DDI | 3did + DOMINE v2 HC | 179 | 1.98% | 4.50% |
| DDI | 3did + DOMINE v2 all | 292 | 3.23% | 7.34% |

Combined coverage across the three best outputs:

- Forward DMI recovered `144` unique PHISTO pairs.
- Reverse DMI recovered `756` unique PHISTO pairs.
- Best DDI resource (`3did + DOMINE v2 all`) recovered `292` unique PHISTO
  pairs.
- The union of those three result sets recovered `1,023` unique PHISTO pairs,
  which is `11.33%` of the full PHISTO benchmark.

Overlap between the best-performing result sets:

- Forward DMI intersect Reverse DMI: `30` pairs
- Forward DMI intersect best DDI: `35` pairs
- Reverse DMI intersect best DDI: `124` pairs

Resource-level observations:

- Reverse DMI substantially outperformed forward DMI on this benchmark.
- `3did + DOMINE v2 all` was the best DDI resource combination.
- The extended structural resources improve recall but remain low-precision on
  the full cross-product protein panel.

## Interpretation

The benchmark shows that extending MicrobioLink with `3did` and `DOMINE`
improves true-positive recovery relative to the more limited legacy resource
combinations, especially for:

- reverse DMI, which is the strongest of the three tested approaches on this
  PHISTO panel,
- DDI when `3did` is combined with the full DOMINE v2 interaction table.

The benchmark also shows persistent limitations:

- DMI precision is very low on a broad, screen-heavy benchmark because many
  motif matches occur by sequence pattern alone.
- DDI recall is constrained by bacterial Pfam coverage in the PHISTO panel.
- The PHISTO benchmark itself is highly skewed toward one large high-throughput
  study, so these results should not be read as a balanced estimate across
  literature-curated mechanistic HMIs.

## Reproducibility

The exact script snapshot used for the run is:

- `06_scripts/benchmark_phisto_extended.py`
- `06_scripts/package_benchmark_workflow.py`

The exact command used is recorded in:

- `07_metadata/run_command.tsv`
- `07_metadata/figure_build_command.tsv`

The Python version and run date are recorded in:

- `07_metadata/software_environment.tsv`

## Detailed Documentation

See:

- `00_docs/results_guide.md` for table-by-table interpretation guidance
- `00_docs/file_manifest.tsv` for a curated file inventory
- `00_docs/resource_provenance.tsv` for source and version provenance
- `09_notebooks/README.md` for the interactive rerun workflow
- `08_figures/README.md` for figure-level usage notes
