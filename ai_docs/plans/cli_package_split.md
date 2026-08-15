# CLI Package Split — Implementation Plan

Refactor the single `microbiolink/cli.py` (~955 lines, one argparse entry point per module) into a
`microbiolink/cli/` sub-package with **one file per module**, mirroring the existing `microbiolink/workflow/`
layout so a developer finds "the CLI for module X" at `cli/<module>.py`, right next to `workflow/<module>.py`.

This is a **pure structural move**: no behaviour, argument, help text, or entry-point *name* changes. Every
console script keeps its current name and semantics; only the dotted path each script points at changes.

## Context

`microbiolink/cli.py` holds all ten argparse entry points (Modules 1–10) plus their private helpers. The
functions already map **1:1 onto `workflow/`** — same module boundaries, same names. There are no imports of
`microbiolink.cli` anywhere except the ten `[project.scripts]` entries in `pyproject.toml`, and there is no
test suite referencing it, so the blast radius is `cli.py` itself plus those ten lines.

This refactor was agreed in conversation; it fulfils the "separate cli folder which contains argparse and any
functions which are cli specific boilerplate" intent from @ai_docs/plans/microbiolink_refactoring.md (Overall
Layout). The core `workflow/` modules are untouched — the CLI layer stays a thin argparse-only wrapper over
them, per the same plan.

## Decisions

1. **`cli.py` becomes the `cli/` package; one file per module, named to match `workflow/`.** File names are
   identical to their `workflow/` counterparts: `zscore_filter.py`, `membrane_filter.py`, `fasta_download.py`,
   `domain_download.py`, `ddi.py`, `dmi.py`, `idr_filter.py`, `monte_carlo.py`, `tiedie.py`, `enrichment.py`.
2. **Each module's public entry point is renamed to `main`.** So `microbiolink.cli.ddi:main`,
   `microbiolink.cli.monte_carlo:main`, etc. This makes every entry point uniform and self-locating, and
   removes the old name mismatches (`download_fasta`/`download_domains`/`monte_carlo_filter`).
3. **Cross-module helpers live in `microbiolink/cli/_common.py`.** These are the only helpers with more than
   one consumer:
   - `_add_human_identifier_arguments`, `_add_microbial_identifier_arguments` — used by `fasta_download` and
     `domain_download`.
   - `_resolve_human_identifiers`, `_resolve_microbial_identifiers` — used by `fasta_download` and
     `domain_download`.
   - `_read_domain_mapping` — used by `ddi` and `dmi`.
   The leading `_` module name marks `_common` internal to the package; siblings import from it directly
   (e.g. `from ._common import _add_human_identifier_arguments`). The helper functions **keep their leading
   underscore** — names are moved verbatim, no renames. Single-consumer helpers likewise keep their leading
   underscore inside their own module file.
4. **`cli/__init__.py` stays minimal — a module docstring only, no re-exports.** Entry points reference each
   `module:main` directly, so there is deliberately no central registry of every command (recreating one
   would undercut the point of the split).
5. **`pandas` is imported where used, not package-wide.** The current top-level `import pandas as pd` moves
   into only the module files that touch pandas: `_common.py`, `membrane_filter.py`, `idr_filter.py`,
   `monte_carlo.py`, `tiedie.py`, `enrichment.py`. Files that only call `result.to_csv(...)` on a returned
   frame (`ddi.py`, `dmi.py`) or don't touch pandas at all (`zscore_filter.py`, `fasta_download.py`,
   `domain_download.py`) do not import it. Lazy `from .workflow import ...` / `from .utils import ...`
   imports inside function bodies are preserved verbatim — note these become `from ..workflow import ...` /
   `from ..utils import ...` (one level deeper) now that the modules sit inside the `cli/` sub-package.
6. **Each module file gets an `if __name__ == "__main__": raise SystemExit(main())` guard** (new — the old
   `cli.py` had none). This makes every command runnable as `python -m microbiolink.cli.<module> ...` during
   development without reinstalling, alongside the installed console script. `_common.py` (no entry point)
   gets no guard.
7. **The change lands as one atomic commit** on the `refactoring` branch: create the `cli/` package (ten
   modules + `_common.py`), delete `cli.py`, and repoint the ten `[project.scripts]` paths together. A
   per-module commit series is not possible — `cli.py` and a `cli/` package cannot coexist, so the file must
   convert to a package in a single step.
8. **The nine per-module plans (`module_1_*` … `module_10_*`) are left unchanged.** They document the old
   `microbiolink.cli:<fn>` path as dated, checked-off build records; this plan is the source of truth for the
   new entry-point paths. History is not rewritten.

## Target layout

```
microbiolink/cli/
  __init__.py              # docstring only
  _common.py               # _add_human/microbial_identifier_arguments,
                           #   _resolve_human/microbial_identifiers, _read_domain_mapping   (imports pandas)
  zscore_filter.py         # main
  membrane_filter.py       # _add_identifier_source_arguments, _add_membrane_filter_target_arguments,
                           #   _build_membrane_filter_parser, _resolve_membrane_filter_identifiers, main
                           #   (imports pandas)
  fasta_download.py        # _build_fasta_download_parser, main
  domain_download.py       # _build_domain_download_parser, _write_domain_table, main
  ddi.py                   # _build_ddi_parser, main
  dmi.py                   # _build_dmi_parser, main
  idr_filter.py            # _add_idr_filter_source_arguments, _add_idr_filter_scoring_arguments,
                           #   _build_idr_filter_parser, main   (imports pandas)
  monte_carlo.py           # _add_monte_carlo_source_arguments, _add_monte_carlo_disorder_arguments,
                           #   _add_monte_carlo_test_arguments, _build_monte_carlo_parser, main
                           #   (imports pandas)
  tiedie.py                # _add_tiedie_source_arguments, _add_tiedie_algorithm_arguments,
                           #   _build_tiedie_parser, _read_expressed_genes, main   (imports pandas)
  enrichment.py            # _add_enrichment_source_arguments, _build_enrichment_parser, main
                           #   (imports pandas)
```

### Function-to-file map (from the current `cli.py`)

| Current function | Destination | New name |
| --- | --- | --- |
| `zscore_filter` | `cli/zscore_filter.py` | `main` |
| `_add_identifier_source_arguments` | `cli/membrane_filter.py` | unchanged |
| `_add_membrane_filter_target_arguments` | `cli/membrane_filter.py` | unchanged |
| `_build_membrane_filter_parser` | `cli/membrane_filter.py` | unchanged |
| `_resolve_membrane_filter_identifiers` | `cli/membrane_filter.py` | unchanged |
| `membrane_filter` | `cli/membrane_filter.py` | `main` |
| `_add_human_identifier_arguments` | `cli/_common.py` | unchanged |
| `_add_microbial_identifier_arguments` | `cli/_common.py` | unchanged |
| `_resolve_human_identifiers` | `cli/_common.py` | unchanged |
| `_resolve_microbial_identifiers` | `cli/_common.py` | unchanged |
| `_build_fasta_download_parser` | `cli/fasta_download.py` | unchanged |
| `download_fasta` | `cli/fasta_download.py` | `main` |
| `_build_domain_download_parser` | `cli/domain_download.py` | unchanged |
| `_write_domain_table` | `cli/domain_download.py` | unchanged |
| `download_domains` | `cli/domain_download.py` | `main` |
| `_build_ddi_parser` | `cli/ddi.py` | unchanged |
| `_read_domain_mapping` | `cli/_common.py` | unchanged |
| `ddi` | `cli/ddi.py` | `main` |
| `_build_dmi_parser` | `cli/dmi.py` | unchanged |
| `dmi` | `cli/dmi.py` | `main` |
| `_add_idr_filter_source_arguments` | `cli/idr_filter.py` | unchanged |
| `_add_idr_filter_scoring_arguments` | `cli/idr_filter.py` | unchanged |
| `_build_idr_filter_parser` | `cli/idr_filter.py` | unchanged |
| `idr_filter` | `cli/idr_filter.py` | `main` |
| `_add_monte_carlo_source_arguments` | `cli/monte_carlo.py` | unchanged |
| `_add_monte_carlo_disorder_arguments` | `cli/monte_carlo.py` | unchanged |
| `_add_monte_carlo_test_arguments` | `cli/monte_carlo.py` | unchanged |
| `_build_monte_carlo_parser` | `cli/monte_carlo.py` | unchanged |
| `monte_carlo_filter` | `cli/monte_carlo.py` | `main` |
| `_add_tiedie_source_arguments` | `cli/tiedie.py` | unchanged |
| `_add_tiedie_algorithm_arguments` | `cli/tiedie.py` | unchanged |
| `_build_tiedie_parser` | `cli/tiedie.py` | unchanged |
| `_read_expressed_genes` | `cli/tiedie.py` | unchanged |
| `tiedie` | `cli/tiedie.py` | `main` |
| `_add_enrichment_source_arguments` | `cli/enrichment.py` | unchanged |
| `_build_enrichment_parser` | `cli/enrichment.py` | unchanged |
| `enrichment` | `cli/enrichment.py` | `main` |

### `pyproject.toml` entry-point rewrite

`[project.scripts]` — names unchanged, paths re-pointed at `module:main`:

```toml
microbiolink-ddi = "microbiolink.cli.ddi:main"
microbiolink-enrichment = "microbiolink.cli.enrichment:main"
microbiolink-dmi = "microbiolink.cli.dmi:main"
microbiolink-download-domains = "microbiolink.cli.domain_download:main"
microbiolink-download-fasta = "microbiolink.cli.fasta_download:main"
microbiolink-idr-filter = "microbiolink.cli.idr_filter:main"
microbiolink-membrane-filter = "microbiolink.cli.membrane_filter:main"
microbiolink-monte-carlo = "microbiolink.cli.monte_carlo:main"
microbiolink-tiedie = "microbiolink.cli.tiedie:main"
microbiolink-zscore-filter = "microbiolink.cli.zscore_filter:main"
```

`[tool.hatch.build.targets.wheel]` already ships `packages = ["microbiolink"]`, so the new `cli/`
sub-package is picked up automatically — no packaging-config change beyond the scripts.

## Migration checklist

### 1. Create the package
- [ ] Create `microbiolink/cli/` with `__init__.py` (docstring only).
- [ ] Create `cli/_common.py` with the five shared helpers (names unchanged), `import pandas as pd`,
      and the lazy `from ..utils import uniprot_client` inside the resolver bodies.

### 2. Move each module's entry point
- [ ] For each of the ten modules, create `cli/<module>.py` with its private helpers (underscore kept),
      the entry point renamed to `main`, `import argparse` (+ `pathlib.Path` for `domain_download`,
      + `pandas` where used), `from ._common import ...` where a shared helper is needed
      (`fasta_download`, `domain_download`, `ddi`, `dmi`), and a trailing
      `if __name__ == "__main__": raise SystemExit(main())` guard.
- [ ] Fix the depth of every lazy in-body import: `from .workflow` → `from ..workflow`,
      `from .utils` → `from ..utils`.

### 3. Delete the old file + rewire packaging
- [ ] Delete `microbiolink/cli.py`.
- [ ] Rewrite the ten `[project.scripts]` paths to `microbiolink.cli.<module>:main`.

### 4. Verification (see Verification section)
- [ ] `uv run ruff format` / `uv run ruff check` over `microbiolink/cli/` and `pyproject.toml`; `uv run ty check`.
- [ ] Import every module — `python -c "import microbiolink.cli.<module>"` for all ten — to catch
      top-level import errors (argparse, pandas, `_common` wiring).
- [ ] `uv pip install -e .` (or `uv sync`) so the console scripts regenerate, then `--help` on all ten.
- [ ] **One real end-to-end run** — `microbiolink-zscore-filter -i
      case_study_input/input/human_transcriptomics/colon_BEST4_enterocyte_CD.csv -zscore <cut> -o <tmp>` —
      to prove the deepened lazy `..workflow`/`..utils` imports resolve at runtime (a `--help` alone does not
      execute the function body).
- [ ] `grep -rn "microbiolink\.cli\b\|from \.cli\|import cli" microbiolink` returns nothing.

### 5. Commit
- [ ] Single atomic commit on `refactoring` (Decision 7): new `cli/` package, deleted `cli.py`, repointed
      `[project.scripts]`.

## Verification

The refactor changes no transformation logic, so the case-study regression fixtures are not re-diffed for
correctness. Confidence comes from: (a) `ruff`/`ty` clean; (b) importing all ten `cli.<module>` modules
(top-level imports incl. `_common`); (c) `--help` on all ten (parser construction through the new entry
points); (d) **one real `microbiolink-zscore-filter` run** against a `case_study_input/` count matrix — the
only step that actually executes a function body and therefore the deepened lazy `..workflow`/`..utils`
imports, which the mechanical `.`→`..` change makes representative of all ten; and (e) the grep confirming the
old module path is fully retired.

## Future Work

- **Shared CLI-IO helpers** — `_read_domain_mapping` / `_read_expressed_genes` and the first-column symbol
  reads scattered across modules could consolidate into a small `cli/_io.py` if more commands accrue; kept in
  `_common.py` / their module for now to avoid over-abstracting a stable set of ten commands.
