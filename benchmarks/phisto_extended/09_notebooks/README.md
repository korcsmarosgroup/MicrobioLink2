# Interactive Walkthrough

This folder contains the package-driven walkthrough notebook for the broad
PHISTO benchmark:

- `phisto_extended_package_walkthrough.ipynb`

The notebook is designed to be read top to bottom and covers:

1. rebuilding the unique PHISTO true-positive panel from the raw export,
2. exporting accession lists for the benchmark panel,
3. optionally redownloading FASTA and Pfam annotations with the new generic
   MicrobioLink helpers,
4. previewing the forward DMI, reverse DMI, and DDI schemas from the packaged
   API,
5. rerunning the full PHISTO benchmark through the packaged API into
   `../10_reproduced_run/`.

Important operational note:

- The documented benchmark outputs in `../05_results/` are the canonical
  reference bundle shipped with this branch.
- The notebook writes any regenerated package-based outputs into
  `../10_reproduced_run/` so the original documented bundle stays untouched.
- Full reverse DMI on the PHISTO panel is large. The support module batches the
  package calls to avoid holding the entire prediction cross-product in memory
  at once.
