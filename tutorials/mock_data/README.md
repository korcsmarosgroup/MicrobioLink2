# Mock-Data Tutorials

This folder contains small runnable tutorials that demonstrate the packaged
MicrobioLink code on fully mocked inputs.

The notebooks are designed to be portable:

- they find the repository root dynamically,
- they write all generated files under `tutorials/mock_data/runs/`,
- they do not depend on user-specific absolute paths.

## Tutorials

- `01_mock_dmi.ipynb`
  Runs the domain-motif interaction step on a tiny mock human FASTA and mock
  bacterial Pfam table. This notebook only needs the base package install.

- `02_mock_dmi_aiupred_monte_carlo.ipynb`
  Extends the mock DMI example with an optional AIUPred run and a deterministic
  Monte Carlo motif filter. The AIUPred section requires `pip install ".[idr]"`.
  The Monte Carlo section still runs without AIUPred because it creates a small
  residue-level mock structure-score table.

## Recommended setup

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -e .
```

For the AIUPred notebook:

```bash
pip install -e ".[idr]"
```
