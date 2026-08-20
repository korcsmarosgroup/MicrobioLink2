<p align="center">
  <img src="microbiolink_logo.png" alt="MicrobioLink" width="360">
</p>

<h1 align="center">MicrobioLink</h1>

<p align="center"><strong>Predicting host-microbe interactions and their downstream effects on host cells.</strong></p>

<p align="center">
  <a href="https://microbiolink2.readthedocs.io/"><img src="https://readthedocs.org/projects/microbiolink2/badge/?version=latest" alt="Documentation Status"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-BSD_2--Clause-blue.svg" alt="License: BSD 2-Clause"></a>
  <img src="https://img.shields.io/badge/python-3.12%2B-blue.svg" alt="Python 3.12+">
</p>

MicrobioLink is a computational pipeline that predicts host-microbe protein-protein
interactions and traces their downstream effects on host cellular signalling. By
integrating multi-omic data with network-biology approaches, it reveals how microbial
proteins engage host proteins and how those interactions ripple through host signalling
in health and disease, such as inflammatory bowel disease (IBD).

## 📖 Documentation

**Full documentation is at [microbiolink2.readthedocs.io](https://microbiolink2.readthedocs.io/).**

- [Introduction](https://microbiolink2.readthedocs.io/) — what MicrobioLink does and its key features.
- [Get Started](https://microbiolink2.readthedocs.io/get-started/) — every install method and its extras.
- [Pipeline Concepts](https://microbiolink2.readthedocs.io/concepts/) — the ten modules, and forward vs reverse MicrobioLink.
- [API Reference](https://microbiolink2.readthedocs.io/api/) — the programmatic surface.

## Quick install

MicrobioLink requires **Python 3.12+**. Install the latest code straight from GitHub:

```bash
# with uv (recommended)
uv pip install "git+https://github.com/korcsmarosgroup/MicrobioLink2.git"

# or with pip
pip install "git+https://github.com/korcsmarosgroup/MicrobioLink2.git"
```

Some modules need optional extras — `enrichment`, `idr`, and `tiedie`:

```bash
uv pip install "microbiolink[idr,tiedie,enrichment] @ git+https://github.com/korcsmarosgroup/MicrobioLink2.git"
```

See [Get Started](https://microbiolink2.readthedocs.io/get-started/) for wheel/sdist installs and the full details.

## Citation

If you use MicrobioLink in your research, please cite Gul et al., 2025 (STAR Protocols)
and Gul et al., 2022 (Journal of Extracellular Vesicles).

## License

Released under the [BSD 2-Clause License](LICENSE).
