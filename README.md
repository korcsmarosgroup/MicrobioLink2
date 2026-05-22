# MicrobioLink

This repository is the clean installable package-only distribution for
MicrobioLink.

It contains only:

- the runtime Python packages `microbiolink` and `microbiolink_api`,
- the packaged ELM resource files used by `microbiolink_api`,
- the package metadata needed for `pip install`,
- a small `tutorials/mock_data` folder with runnable mock-data notebooks,
- the license file.

It does not contain:

- documentation site files,
- tests,
- generated example inputs,
- development tooling configuration.

## Install

Create and activate an environment first, then install from the repository
root:

```bash
python3 -m venv .venv
source .venv/bin/activate

pip install -U pip
pip install .
```

## Optional installs

Install AIUPred / IDR support:

```bash
pip install ".[idr]"
```

Install TieDIE workflow support:

```bash
pip install ".[workflow]"
```

Notes:

- The `idr` extra installs `biopython`, `torch`, `iupred`, and constrains
  `numpy<2` for AIUPred compatibility.
- Python 3.10 or 3.11 is recommended when using the `idr` extra and
  `microbiolink-get-human-fasta`.

## What gets installed

### CLI commands

- `microbiolink-zscore-filter`
- `microbiolink-dmi`
- `microbiolink-download-bacterial-proteins`
- `microbiolink-get-human-fasta`
- `microbiolink-aiupred`
- `microbiolink-idr-prediction`
- `microbiolink-idr-prediction-score`
- `microbiolink-tiedie-input-processing`
- `microbiolink-processing-tiedie-output`
- `microbiolink-enrichr-ranking`

### Python packages

- `microbiolink`
- `microbiolink_api`

## Minimal usage

### CLI

Predict domain-motif interactions from a human FASTA and bacterial domain table:

```bash
microbiolink-dmi \
  --fasta_file human_proteins.fasta \
  --elm_regex_file elm_classes.tsv \
  --motif_domain_file elm_interaction_domains.tsv \
  --bacterial_domain_file bacterial_domains.tsv \
  --output_file dmi_output.tsv
```

### Python API

Use the library-style API with packaged ELM resources:

```python
from microbiolink_api import predict_domain_motif_interactions

interactions = predict_domain_motif_interactions(
    fasta_file = 'human_proteins.fasta',
    bacterial_domain_file = 'bacterial_domains.tsv',
)
```

## Expected user data

Typical user-provided inputs are:

- a host count matrix,
- a host FASTA file,
- a bacterial protein list or proteome ID list,
- or a precomputed bacterial domain table.

## Package layout

```text
.
├── LICENSE
├── README.md
├── pyproject.toml
├── microbiolink/
├── microbiolink_api/
└── tutorials/
```
