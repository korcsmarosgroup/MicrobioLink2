# Plan: MicrobioLink Documentation — Static Site (Phase 1)

Status: proposed
Date: 2026-08-19
Parent: `documentation_build.md` (umbrella). Sibling: `documentation_tutorials.md`.

## Scope

The static documentation site — everything that does **not** require running the MicrobioLink
pipeline. Ships green on Read the Docs on its own, before any tutorial content exists. Owns:
Introduction, Get Started, Pipeline Concepts, API Reference, README trim, the Tutorials landing page,
the four tutorial **nav slots** (placeholders), and the RTD + nbconvert stopgap wiring the tutorials
plan will later use.

Out of scope (see `documentation_tutorials.md`): the curated dataset, any pipeline execution, and the
T1–T4 tutorial content.

## Dependency-light by design

Per the umbrella's verified facts, this build needs only the **base `microbiolink` package + zensical
+ mkdocstrings + nbconvert** — no torch, no `[idr]`/`[tiedie]` git extras. mkdocstrings uses griffe's
static AST analysis, and the heavy deps are imported lazily inside functions, so nothing heavy is
pulled to render the API reference.

## Repository layout (created by this plan)

```
docs/
  index.md                     Introduction (canonical)
  get-started.md
  concepts.md
  tutorials/
    index.md                   Tutorials landing (who each tutorial is for) + links to 4 slots
    t1_api_full.md             placeholder slot (filled by tutorials plan)
    t2_api_basic.md            placeholder slot
    t3_cli_full.md             placeholder slot
    t4_cli_basic.md            placeholder slot
  api/                         mkdocstrings pages, one per module
  requirements.txt             zensical, mkdocstrings, nbconvert
zensical.toml (or mkdocs.yml)  theme, nav, mkdocstrings config
.readthedocs.yaml              RTD build + nbconvert pre_build stopgap
```

## `.readthedocs.yaml` (sketch)

```yaml
version: 2
build:
  os: ubuntu-24.04
  tools:
    python: "3.12"
  jobs:
    pre_build:
      # THROWAWAY STOPGAP — remove when Zensical renders .ipynb natively.
      # No-op until the tutorials plan adds the .ipynb files; wiring lives here from Phase 1.
      - jupyter nbconvert --to markdown docs/tutorials/t1_api_full.ipynb docs/tutorials/t2_api_basic.ipynb || true
python:
  install:
    - method: pip
      path: .                                   # base microbiolink only (light deps)
    - requirements: docs/requirements.txt       # zensical, mkdocstrings, nbconvert
```

## Steps

### Step 1 — Scaffold Zensical + RTD
- Add zensical + mkdocstrings + nbconvert to `docs/requirements.txt` (installed via `uv`).
- Create the `docs/` skeleton and Zensical config with the five-section nav.
- Add `.readthedocs.yaml` with the nbconvert `pre_build` stopgap (guarded so it is a no-op until the
  notebooks land).
- Connect the RTD project to the `refactoring` branch as the `latest` version.

### Step 2 — Introduction (`index.md`)
- Canonical overview + key features, seeded from the current README and extended with the refactored
  functionality described in `../microbiolink_refactoring.md`.
- Four links out: Get Started, Concepts, Tutorials, API.

### Step 3 — Get Started (`get-started.md`)
- Install methods derived from `../package_build.md`: wheel/sdist, `pip install git+…`, optional extras
  `[idr]` and `[tiedie]`, and the ten `microbiolink-*` console scripts.

### Step 4 — Pipeline Concepts (`concepts.md`)
- The 10 modules, forward vs reverse MicrobioLink, and why the optional modules exist. This is the
  single home for the conceptual prose the four tutorials link back to.

### Step 5 — API Reference
- Add `__all__` to each `workflow/*.py` naming the essential (user-facing) functions.
- Configure mkdocstrings pages per module, split **Essential** (`__all__`) vs **Helpers**.
- **Docstring sense-check** over the essential functions so rendered wording reads well.

### Step 6 — Tutorials landing + slots
- Write `docs/tutorials/index.md` describing what each tutorial is and who it is for.
- Add the four placeholder pages / nav slots so the site structure is complete and links resolve.

### Step 7 — Trim the README
- Reduce README to pitch + badges + quick install + a prominent link to the hosted docs. (Also
  clears the stale-install-section flag noted in `../package_build.md`.)

### Step 8 — Verify
- Local `zensical` build succeeds; nav renders all five sections; all internal links resolve.
- API pages render from docstrings; Essential/Helper split matches `__all__`.
- RTD build of `refactoring` (`latest`) goes green using only the light dependency set.

## Acceptance criteria

- Zensical site builds on RTD from `refactoring` (`latest`) with only base `microbiolink` + zensical
  + mkdocstrings + nbconvert — no torch, no git extras.
- Five sections present and navigable; Introduction is the canonical overview and links to all four.
- Every `workflow` module declares `__all__`; API docs show Essential vs Helper accordingly.
- README trimmed to a pointer.
- Tutorial slots exist and resolve (content filled by `documentation_tutorials.md`).
