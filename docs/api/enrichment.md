# Enrichment analysis

Functional enrichment analysis: reduce a network/HMI table to human targets, run Enrichr, plot.

## Essential

The functions you call to run this module.

::: microbiolink.workflow.enrichment
    options:
      members:
        - run_enrichment_analysis
        - run_enrichment
        - plot_enrichment

## Helpers

Internal machinery, documented for completeness — you don't normally call these directly.

::: microbiolink.workflow.enrichment
    options:
      show_root_heading: false
      members:
        - extract_hmi_targets
        - extract_tiedie_nodes
