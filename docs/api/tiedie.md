# TieDie network propagation

TieDie tied-diffusion network propagation: build inputs, run TieDie, format the output network.

## Essential

The functions you call to run this module.

::: microbiolink.workflow.tiedie
    options:
      members:
        - run_tiedie_pipeline

## Helpers

Internal machinery, documented for completeness — you don't normally call these directly.

::: microbiolink.workflow.tiedie
    options:
      show_root_heading: false
      members:
        - fetch_omnipath_networks
        - combine_host_microbe_interactions
        - contextualise_networks
        - build_pathway_sif
        - build_upstream_heats
        - build_downstream_heats
        - run_tiedie
        - assemble_network
        - assemble_node_table
