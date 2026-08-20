# MicrobioLink

**A tool for predicting host-microbe interactions and their downstream effects on host cells.**

MicrobioLink is a computational pipeline that predicts host-microbe protein-protein
interactions and analyses their downstream effects on host cellular signalling
pathways. By integrating multi-omic data with network-biology approaches, it reveals
how microbial proteins engage host proteins and how those interactions ripple through
host signalling to affect homeostasis and disease.

It is built for researchers investigating complex interspecies interactions and their
implications in health and disease, such as inflammatory bowel disease (IBD).

!!! note "New to MicrobioLink?"
    This page is the canonical overview. Once you know what MicrobioLink does, head to
    [Get Started](get-started.md) to install it, or [Pipeline Concepts](concepts.md) to
    understand how the modules fit together.

## What MicrobioLink does

Starting from host transcriptomics and a set of microbial proteins, MicrobioLink
builds a predicted host-microbe interaction network and traces its influence into the
host cell:

- **Predicts host-microbe protein-protein interactions** using structural data and
  domain-motif interactions between microbial and host proteins.
- **Refines predictions with disorder and binding-site analysis** to raise confidence
  in the interactions that survive filtering.
- **Integrates multi-omic data**, combining transcriptomic and proteomic inputs into a
  single model of host-microbe interaction.
- **Analyses downstream signalling**, modelling how microbial interactions propagate
  through host cellular processes and which pathways are impacted.
- **Exports networks for visualisation** in Cytoscape for interactive exploration of
  the results.

## The pipeline modules

The refactored pipeline is organised into ten modules that can be run end to end or
composed individually. Each takes the output of the previous step, so you can enter the
pipeline wherever your data allows.

1. **Z-score Filter** — removes lowly expressed genes from a count matrix using a
   user-defined z-score cut-off.
2. **Membrane Protein Filter** — restricts proteins to the cell surface: human proteins
   via OmniPath InterCell, microbial proteins via UniProt (secreted, outer- or
   plasma-membrane).
3. **FASTA Download** — retrieves protein sequences from UniProt for IDs, gene symbols,
   or a whole UniProt proteome.
4. **Domain Download** — finds the Pfam domains carried by each protein.
5. **Domain-Domain Interactions (DDI)** — predicts interactions between microbial and
   host domains using the DOMINE and 3did resources.
6. **Domain-Motif Interactions (DMI)** — connects domains to linear motifs using ELM
   and 3did, in the forward (microbial domain to host motif) or reverse direction, or
   both.
7. **Intrinsic Disorder Region (IDR) Prediction** — filters motif interactions by the
   likelihood the motif lies in a disordered, binding-competent region, using IUPred /
   AIUPred.
8. **Monte Carlo Simulation** — tests whether a motif's position is more likely to bind
   the domain than a random position, adding a significance score.
9. **TieDie** — propagates signal across the host network from the predicted
   interactions to a set of differentially expressed genes, producing the final
   downstream-effect network.
10. **Functional Enrichment Analysis** — runs enrichment over the TieDie network or the
    host targets of microbial proteins, returning a table and an enrichment plot.

For how the modules combine, forward versus reverse MicrobioLink, and which are
optional, see [Pipeline Concepts](concepts.md).

## Inputs and outputs

MicrobioLink works from three kinds of input:

- **Host transcriptomics** — a gene count matrix with gene symbols and expression
  values.
- **A gene endpoint set** — the target or differentially expressed genes to analyse.
- **Microbial proteins** — a list of UniProt IDs or a UniProt proteome (UP) ID.

The result is a predicted, filtered host-microbe interaction network and its downstream
signalling effects, exportable to Cytoscape for interactive visualisation.

## Where to next

<div class="grid cards" markdown>

- **[Get Started](get-started.md)**

    Install MicrobioLink and run your first analysis.

- **[Pipeline Concepts](concepts.md)**

    Understand the ten modules and how they connect.

- **[Tutorials](tutorials/index.md)**

    Worked, end-to-end examples.

- **[API Reference](api/index.md)**

    The programmatic surface for running MicrobioLink without the CLI.

</div>

## Citation

If you use MicrobioLink in your research, please cite Gul et al., 2025 (STAR Protocols)
and Gul et al., 2022 (Journal of Extracellular Vesicles).
