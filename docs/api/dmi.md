# Domain–motif interactions (DMI)

Predicting domain-motif interactions between bacterial and human proteins.

## Essential

The functions you call to run this module.

::: microbiolink.workflow.dmi
    options:
      members:
        - predict_domain_motif_interactions
        - load_motif_regexes

## Helpers

Internal machinery, documented for completeness — you don't normally call these directly.

::: microbiolink.workflow.dmi
    options:
      show_root_heading: false
      members:
        - MotifResource
        - RoleBasedMatch
