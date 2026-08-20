# API Reference

The programmatic surface of MicrobioLink: the ten `microbiolink.workflow` modules, generated from their docstrings. Each module page is split into **Essential** functions (the ones you call to run that step of the pipeline) and **Helpers** (internal machinery, shown for completeness).

The reference covers `microbiolink.workflow` only. The command-line entry points (`cli/`) are covered in the [Tutorials](../tutorials/index.md), and shared utilities (`utils/`) are internal. For how the modules fit together, see [Pipeline Concepts](../concepts.md).

## Modules

In pipeline order:

1. [Z-score filter](zscore_filter.md) — `microbiolink.workflow.zscore_filter`
2. [Membrane filter](membrane_filter.md) — `microbiolink.workflow.membrane_filter`
3. [FASTA download](fasta_download.md) — `microbiolink.workflow.fasta_download`
4. [Domain download](domain_download.md) — `microbiolink.workflow.domain_download`
5. [Domain–domain interactions (DDI)](ddi.md) — `microbiolink.workflow.ddi`
6. [Domain–motif interactions (DMI)](dmi.md) — `microbiolink.workflow.dmi`
7. [IDR (disorder) filter](idr_filter.md) — `microbiolink.workflow.idr_filter`
8. [Monte Carlo significance filter](monte_carlo.md) — `microbiolink.workflow.monte_carlo`
9. [TieDie network propagation](tiedie.md) — `microbiolink.workflow.tiedie`
10. [Enrichment analysis](enrichment.md) — `microbiolink.workflow.enrichment`
