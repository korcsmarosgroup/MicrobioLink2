# IDR (disorder) filter

Filtering domain-motif interactions by intrinsic disorder and binding likelihood.

## Essential

The functions you call to run this module.

::: microbiolink.workflow.idr_filter
    options:
      members:
        - filter_by_disorder

## Helpers

Internal machinery, documented for completeness — you don't normally call these directly.

::: microbiolink.workflow.idr_filter
    options:
      show_root_heading: false
      members:
        - iupred_profile
        - aiupred_profile
        - cached_profile
        - validate_method
        - disordered_mask
        - disordered_region
