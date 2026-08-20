# Get Started

This page covers how to install **MicrobioLink** and where to go next. MicrobioLink
requires **Python 3.12 or newer**.

Installing the package registers a set of `microbiolink-*` console scripts on your `PATH`
(one per pipeline module). They are documented in
[Pipeline Concepts](concepts.md) and the [Tutorials](tutorials/index.md) — this page only
gets the package onto your machine.

We use [uv](https://docs.astral.sh/uv/) as the primary workflow tool; plain `pip`
equivalents are shown alongside where they differ.

## Install from GitHub

The quickest way to install the latest code is directly from the repository. No local
checkout or build step is required.

=== "uv"

    ```bash
    uv pip install "git+https://github.com/korcsmarosgroup/MicrobioLink2.git"
    ```

=== "pip"

    ```bash
    pip install "git+https://github.com/korcsmarosgroup/MicrobioLink2.git"
    ```

To add MicrobioLink as a dependency of an existing uv project instead:

```bash
uv add "git+https://github.com/korcsmarosgroup/MicrobioLink2.git"
```

## Install from a built wheel or sdist

If you have built the distribution locally (`uv build` produces a wheel and an sdist under
`dist/`), install the artifact directly. Adjust the version in the filename to match what
`uv build` produced.

=== "uv"

    ```bash
    # Wheel (recommended)
    uv pip install dist/microbiolink-2.1.0-py3-none-any.whl

    # Or the source distribution
    uv pip install dist/microbiolink-2.1.0.tar.gz
    ```

=== "pip"

    ```bash
    # Wheel (recommended)
    pip install dist/microbiolink-2.1.0-py3-none-any.whl

    # Or the source distribution
    pip install dist/microbiolink-2.1.0.tar.gz
    ```

## Optional extras

Some pipeline modules depend on extra packages that are **not** installed by default.
Request them with the extras syntax in square brackets. The available extras are:

| Extra | Enables | Notes |
| --- | --- | --- |
| `enrichment` | Functional enrichment (`gget`, `matplotlib`) | Standard PyPI dependencies. |
| `idr` | Intrinsically disordered region filtering (`iupred`) | Installed from a git source. |
| `tiedie` | TieDIE network diffusion (`tiedie`, `networkx`) | Installed from a git source. |

Combine extras in a comma-separated list. Because the `idr` and `tiedie` extras resolve to
git-based dependencies, they only work with local or VCS installs (not a future PyPI
release).

=== "uv"

    ```bash
    # From GitHub, with the idr and tiedie extras
    uv pip install "microbiolink[idr,tiedie] @ git+https://github.com/korcsmarosgroup/MicrobioLink2.git"

    # From a local wheel, with the enrichment extra
    uv pip install "dist/microbiolink-2.1.0-py3-none-any.whl[enrichment]"
    ```

=== "pip"

    ```bash
    # From GitHub, with the idr and tiedie extras
    pip install "microbiolink[idr,tiedie] @ git+https://github.com/korcsmarosgroup/MicrobioLink2.git"

    # From a local wheel, with the enrichment extra
    pip install "dist/microbiolink-2.1.0-py3-none-any.whl[enrichment]"
    ```

## Verify the installation

Confirm the package imports and the console scripts resolve:

```bash
python -c "import microbiolink; print('MicrobioLink import ok')"
microbiolink-ddi --help
```

## Next steps

- [Pipeline Concepts](concepts.md) — the ten modules and how forward and reverse runs fit
  together.
- [Tutorials](tutorials/index.md) — worked, end-to-end examples that put the console scripts
  to use.
