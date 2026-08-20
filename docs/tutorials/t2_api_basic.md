# A default MicrobioLink run, for first-time users (notebook)

**What it will cover:** the default forward-only run through MicrobioLink,
driven from a Jupyter notebook using the public API. It works through the
full depth of the pipeline &mdash; membrane filter, FASTA download, domain
download, DMI, IDR filter, Monte Carlo, TieDie and enrichment (run twice:
once on the initial host proteins and once on the final TieDie network)
&mdash; and drops only the Z-score Filter and DDI modules. It still uses the
`[idr]`, `[tiedie]` and `[enrichment]` extras.

**Who it is for:** first-time users who want a sensible default run in a
notebook, taking the forward path a newcomer is most likely to want without
also setting up the reverse DMI and Z-score branches.

See the [Tutorials overview](index.md) for how this fits the other three, and
[Pipeline Concepts](../concepts.md) for the modules it uses.

!!! note "Content coming in Phase 2"
    This is a placeholder slot. The step-by-step walkthrough is added in
    **Phase 2** of the documentation work.
