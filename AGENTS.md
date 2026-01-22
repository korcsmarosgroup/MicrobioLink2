# CLAUDE.md

This file provides guidance for Claude Code when working on the microbiolink project.

## Project Overview

**microbiolink** is a bioinformatics pipeline that predicts interactions between microbe and host proteins. It analyzes:
- Protein domains
- Short linear motifs (SLiMs)
- Intrinsically disordered regions (IDRs)
- Other structural properties

The pipeline then optimizes the resulting microbe-host interaction network using host transcriptomics data.

### Prior Knowledge Sources
- **UniProt**: Protein sequences and annotations
- **OmniPath**: Signaling network and pathway data
- **ELM** (Eukaryotic Linear Motif): Motif-domain interactions

Note: These databases were previously stored in the repository but have been removed. Database access will be implemented without local storage.

### Network Optimization
Currently uses **TieDIE** for network optimization. This will be replaced with a more general solution in future development.

## Repository Structure

```
microbiolink/
├── microbiolink/           # Main Python module (all scripts here)
│   ├── DMI.py              # Domain-Motif Interaction prediction
│   ├── AIUPred.py          # AIUPred disorder/binding prediction
│   ├── aiupred_lib.py      # AIUPred library functions
│   ├── iupred2a.py         # IUPred2A interface
│   ├── idr_prediction.py   # IDR prediction
│   ├── idr_prediction_score.py
│   ├── tiedie_input_processing.py   # Prepare TieDIE inputs
│   ├── processing_tiedie_output.py  # Process TieDIE results
│   ├── z-score_filter_terminal.py   # Expression z-score filtering
│   ├── download_bacterial_proteins.py
│   ├── get_human_fasta.py
│   ├── enrichr_id_database_ranking.py
│   ├── requirements.txt
│   └── microbiolink_env.yml
├── README.md
├── LICENSE                 # BSD 2-Clause
└── CLAUDE.md
```

## Pipeline Workflow

1. **Input**: Ready gene count matrix (transcriptomics data processing is handled in separate packages)
2. **Z-score filtering**: Filter expression data
3. **DMI prediction**: Predict domain-motif interactions between bacterial and host proteins
4. **IDR filtering**: Filter interactions based on disordered regions (AIUPred)
5. **Network optimization**: Run TieDIE to identify downstream signaling effects
6. **Output**: Optimized host-microbe interaction network

## Development Context

### Current State
- Python scripts exist but lack proper packaging
- No `__init__.py`, `pyproject.toml`, or standard tooling
- Scripts use argparse CLI interfaces

### Planned Development

#### Immediate Goals
1. **Apply cookiecutter/cruft template**: Standardize project structure
2. **Package for PyPI**: Create proper Python package with `pyproject.toml`
3. **Refactor scripts**: Meet coding standards and provide intuitive API

#### Future Goals
- Replace TieDIE with a more general network optimization solution
- Implement database access (UniProt, OmniPath, ELM) without local storage
- Transcriptomics processing will remain in separate packages

## Technical Details

### Dependencies
Core: `pandas`, `numpy`, `scipy`, `mygene`, `omnipath`, `pyfasta`, `biopython`

### Python Version
Currently targets Python 3.9

### License
BSD 2-Clause (KorcsmarosLab)

## Build & Test Commands

(To be added after packaging setup)

## Code Style

Follow PEP8 and Google style guide with the following project-specific rules (see `CODING_STYLE.md` for full details):

### Naming
- Classes: `PascalCase`, resource names as single word (e.g. `ProteinatlasAnnotation`)
- Functions/variables: `snake_case`, resource names as single word (e.g. `proteinatlas_annotations`)
- Resource names never contain underscores (underscore separates primary/secondary resources)

### Blank Lines
Use vertical whitespace extensively:
- 2 blank lines between functions/classes
- 1 blank line before and after blocks
- 1 blank line after function signature, before return

### Indentation of Argument Lists
Break long argument lists with each element on its own line, trailing comma on all except `**kwargs`:
```python
a = function(
    foo = 'foo',
    bar = 'bar',
)
```

### Comprehensions
Multi-line comprehensions: value, `for`, `if` each on new line:
```python
foo = [
    bar
    for bar in baz
    if bar[0]
]
```

### Docstrings
Napoleon (Google) style, opening and closing quotes on their own lines. Use type hints instead of types in docstrings.

### Operators and Spacing
- Spaces around all operators: `a = 1`, `foo = 'bar'`
- Space after commas in sequences

### Quotes
Single quotes preferred.

### Imports
Group in order (blank line between groups):
1. Standard library
2. Typing imports
3. Third party dependencies
4. Package imports
