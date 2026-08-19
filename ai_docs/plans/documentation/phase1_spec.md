> **Spec — Documentation Phase 1: static site.** Synthesized from `ai_docs/plans/documentation/documentation_site.md` and `documentation_build.md`, with the Phase-1 grilling decisions and verified spike config (Zensical 0.0.56, mkdocstrings-python 2.0.7) folded in. Config sketches in the older `documentation_site.md` (its `.readthedocs.yaml`, `docs/requirements.txt`) are **stale** — the decisions below supersede them. RTD versions were bumped to newest-pinned (`ubuntu-26.04` + Python `3.14`) and **re-verified green on stable Python 3.14.6** (`zensical build` clean; `workflow/` API rendered via static griffe with pandas/torch absent).

## Problem Statement

MicrobioLink has no hosted documentation. A prospective user lands on a dense README that mixes an overview, a feature list, install instructions, and usage notes, with no rendered API reference and no guided path from "what is this" to "how do I install it" to "what do the ten pipeline modules do." The refactored package (ten `workflow/` modules, ten `microbiolink-*` console scripts, optional `[idr]`/`[tiedie]` extras) is undocumented outside its docstrings. Anyone evaluating or adopting MicrobioLink has to read source to understand it.

## Solution

A static documentation site, generated with **Zensical** and hosted on **Read the Docs**, that ships **green on its own before any tutorial content exists**. It presents five navigable sections — Introduction, Get Started, Pipeline Concepts, API Reference, and a Tutorials landing page with four placeholder slots — and an API reference auto-generated from the `workflow/` module docstrings, split into user-facing "Essential" functions and internal "Helpers." The README is trimmed to a pitch-and-pointer that sends readers to the hosted docs. The build is deliberately dependency-light: it needs only the base `microbiolink` package plus `zensical`, `mkdocstrings-python`, and `nbconvert` — no torch, no git extras — because the heavy dependencies are imported lazily inside functions and mkdocstrings reads source statically via griffe.

This is **Phase 1 of two**. It owns the static site and the nbconvert stopgap wiring; Phase 2 (`documentation_tutorials.md`) fills the four tutorial slots with a curated dataset and seeded pipeline runs. Phase 1 does not block on Phase 2 — that is the point of the split.

## User Stories

1. As a prospective MicrobioLink user, I want a hosted documentation site, so that I can evaluate the tool without reading its source code.
2. As a researcher new to MicrobioLink, I want an Introduction page that explains what the tool does and its key features, so that I can decide in a minute whether it fits my problem.
3. As a reader on the Introduction page, I want prominent links to Get Started, Pipeline Concepts, Tutorials, and API Reference, so that I can jump to whichever entry point matches my need.
4. As someone ready to try MicrobioLink, I want a Get Started page listing every install method (wheel/sdist, `pip install git+…`, and the `[idr]`/`[tiedie]` extras), so that I can install the exact configuration I need.
5. As a user who has just installed MicrobioLink, I want the Get Started page to point me to Pipeline Concepts and the Tutorials as next steps, so that I know where to go after installing (the ten `microbiolink-*` console scripts are documented under Concepts / the CLI tutorials, not on Get Started).
6. As a user learning the pipeline, I want a single Pipeline Concepts page describing all ten modules, forward vs reverse MicrobioLink, and why the optional modules exist, so that I have one authoritative conceptual reference instead of scattered explanations.
7. As a developer integrating MicrobioLink programmatically, I want an API Reference generated from docstrings, so that I can see each function's signature, parameters, and returns without opening source.
8. As a developer, I want the API reference split into "Essential" functions (the ones I call to run a module) and "Helpers" (internal machinery), so that I am not distracted by functions I never touch.
9. As a developer, I want the API reference limited to the `workflow/` modules, so that CLI and utility internals do not clutter the programmatic surface.
10. As a reader of the API reference, I want docstring wording that reads clearly (Google-style Args/Returns rendered as structured Parameters/Returns), so that the generated pages are actually usable.
11. As a visitor, I want a Tutorials landing page describing what each of the four tutorials is and who it is for, so that I can pick the right one — even before the tutorial bodies exist.
12. As a visitor, I want the four tutorial nav slots to exist and resolve (as placeholders), so that the site structure is complete and no links are broken.
13. As a maintainer, I want the documentation to build on Read the Docs from the `refactoring` branch as the `latest` version, so that the hosted site tracks current development.
14. As a maintainer, I want the RTD build to succeed using only the light dependency set, so that builds stay within RTD Community limits (15 min, 7 GB RAM) and never pull torch or git extras.
15. As a maintainer, I want the nbconvert stopgap wired in now as a no-op, so that Phase 2 can drop notebooks into the slots without touching build configuration.
16. As a maintainer, I want the docs dependencies declared as a uv `docs` dependency group, so that the environment is reproducible and consistent with the repo's uv-only guideline.
17. As a first-time visitor arriving from GitHub, I want the README trimmed to a pitch, badges, a quick install, and a prominent docs link, so that I am guided to the full documentation rather than reading a stale wall of text.
18. As a package author, I want each `workflow/*.py` module to declare `__all__`, so that the Essential/Helper split is driven by an explicit, reviewed list rather than by naming conventions.
19. As a reviewer of the API surface, I want to approve the proposed Essential set per module before it is wired in, so that the public surface reflects intent, not a guess.
20. As a maintainer preparing for release, I want the RTD `latest` version to be repointable from `refactoring` to `case-study` (main) on merge, so that the hosted docs follow the release branch later without rework.
21. As a contributor, I want the site's information architecture and nav to be complete and stable, so that Phase 2 tutorial content drops into existing slots without restructuring.

## Implementation Decisions

**Generator & hosting**
- **Zensical** is the static-site generator (MkDocs 1.x is EOL; Zensical is the sanctioned drop-in successor). Config file is **`zensical.toml`** (TOML) — not `mkdocs.yml`. `nav` is an array of single-key tables.
- Hosted on **Read the Docs**, single **`latest`** version built from the **`refactoring`** branch now; repointed to `case-study` on merge. Per-tag versioning deferred (additive later).
- `site_url = "https://microbiolink2.readthedocs.io/"` — literal string; `zensical.toml` has no variable interpolation. `site_name = "MicrobioLink"`.
- Baked-in config: `python: "3.14"` (newest stable, RTD-supported), logo/favicon from `microbiolink_logo.png`.

**Dependencies**
- Docs deps declared as a uv **`docs` dependency group** — `zensical`, `mkdocstrings-python`, `nbconvert` — added via `uv add --group docs …` (never pip; repo guideline + active hook). Supersedes the stale `docs/requirements.txt` approach.
- The build is dependency-light **by design**: base `microbiolink` + the three docs tools only. No torch, no `[idr]`/`[tiedie]` git extras. Heavy deps (`iupred`, `tiedie`) import lazily inside functions; mkdocstrings/griffe is static AST analysis — verified in the spike that `workflow/` API pages render with pandas *not* installed.

**Information architecture (five sections)**
- **Introduction** (`docs/index.md`) — canonical overview + key features (seeded from the current README, extended with refactored functionality per `../microbiolink_refactoring.md`); four links out to Get Started, Concepts, Tutorials, API.
- **Get Started** (`docs/get-started.md`) — install methods derived from `../package_build.md`: wheel/sdist, `pip install git+…`, and the `[idr]`/`[tiedie]` extras. Ends with a **Next steps** section linking to Pipeline Concepts and the Tutorials. It does **not** enumerate the ten `microbiolink-*` console scripts — those belong in Concepts / the CLI tutorials.
- **Pipeline Concepts** (`docs/concepts.md`) — the ten modules, forward vs reverse MicrobioLink, why the optional modules exist. Written once; tutorials link back here. Top-level nav section.
- **Tutorials landing** (`docs/tutorials/index.md`) — describes each of T1–T4 and who it is for, plus the four placeholder slots (`t1_api_full`, `t2_api_basic`, `t3_cli_full`, `t4_cli_basic`) so the structure is complete and links resolve. Slot **content** is Phase 2.
- **API Reference** (`docs/api/`, one page per module) — mkdocstrings pages over the ten `workflow/` modules only. `utils/` and `cli/` are excluded (CLI is covered by the CLI tutorials).

**API reference generation**
- mkdocstrings-python handler with `paths = ["."]` (repo root; flat layout), `docstring_style = "google"`, `show_root_heading = true`, `show_source = false`. Docstrings are **Google style**.
- Essential vs Helper split driven by **`__all__` per `workflow/*.py`**: an Essential page selects members explicitly; a full/Helper page omits the member block. Zero runtime effect. Verified in the spike: member selection drives the split, and underscore-privates are hidden by default.
- `__all__` process: the implementing agent **proposes** the essential set per module (derived from what each `cli/*.py` entry calls plus pipeline-entry functions), the **user approves**, then it is wired in. Not yet done.
- A **docstring sense-check** pass over the essential functions so rendered wording reads well.

Verified config from the spike (encodes decisions precisely; inline where the spec builds it):

```toml
# zensical.toml — nav is an array of single-key tables
[project]
site_url  = "https://microbiolink2.readthedocs.io/"
site_name = "MicrobioLink"
nav = [
  { "Introduction" = "index.md" },
  # Get Started, Concepts, Tutorials, API …
]

[project.plugins.mkdocstrings.handlers.python]
paths = ["."]
options.docstring_style = "google"
options.show_root_heading = true
options.show_source = false
```

```markdown
<!-- API page: Essential = selected members; omit the block for full/Helper -->
::: microbiolink.workflow.membrane_filter
    options:
      members:
        - filter_membrane_proteins
```

**Read the Docs build**
- RTD builds Zensical via a full **`build.commands`** override (uv-driven), **not** the mkdocs-native `python.install` style. The nbconvert step is a throwaway stopgap command line (Zensical cannot render `.ipynb` today), guarded so it is a no-op until the Phase 2 notebooks exist. Delete the stopgap when Zensical ships native `.ipynb` support.

```yaml
version: 2
build:
  os: ubuntu-26.04
  tools: { python: "3.14" }
  commands:
    - pip install uv
    - uv sync --group docs
    - uv run jupyter nbconvert --to markdown docs/tutorials/t1_api_full.ipynb docs/tutorials/t2_api_basic.ipynb || true  # throwaway stopgap; no-op until notebooks exist (Phase 2)
    - uv run zensical build
    - mkdir -p $READTHEDOCS_OUTPUT/html && cp -r site/* $READTHEDOCS_OUTPUT/html/
```

**README**
- Trim to a **minimal pointer**: logo + one-paragraph pitch + badges + quick install (uv/pip + `[idr]`/`[tiedie]` extras) + docs link + citation + license. Convert the current HTML to Markdown. Features / Usage / Troubleshooting move into the docs Introduction. This also clears the stale-install-section flag noted in `../package_build.md`.

**Repo convention**
- Internal docs live under `ai_docs/` (`plans/`, `decisions/`, `agents/`); `docs/` is reserved exclusively for the published Zensical site. No Zensical exclude config is needed because `docs/` contains only publishable pages.

**Human-only action (cannot be done by the implementing agent)**
- Create/connect the **`microbiolink2` Read the Docs project** to the GitHub repo and point `latest` at the `refactoring` branch. (The `ready-for-agent` work stops at a green *local* build + committed config; the cloud build is verified by the maintainer.)

## Testing Decisions

The build **is** the test. For a documentation site the natural, highest seam is the `zensical build` command, not invented Python unit tests. A good test here asserts externally observable build behavior, not internal structure:

1. **Green local build with light deps only.** `uv run zensical build` (with just the `docs` group + base `microbiolink` installed, no torch/git extras) exits 0. This is the primary regression guard: if any heavy dependency leaks into the render path, the build breaks.
2. **All five sections render.** The generated `site/` contains pages for Introduction, Get Started, Pipeline Concepts, the Tutorials landing, and the API Reference.
3. **Links resolve.** No broken internal links; the four tutorial slots exist and resolve as placeholders.
4. **API pages render from docstrings** and reflect the `__all__` Essential/Helper split — Google-style Args/Returns appear as structured Parameters/Returns (proven feasible in the spike with pandas absent).
5. **RTD cloud build green** on `refactoring` (`latest`) — verified by the maintainer after connecting the RTD project (the human-only action above).

No automated CI/pytest link-checker is in scope for Phase 1 (the maintainer opted for the build-is-the-seam approach); the checks above are run against the build output. Prior art: none in-repo — this is the first docs build; the spike (Zensical 0.0.56, mkdocstrings-python 2.0.7) is the reference that these assumptions were verified against.

## Out of Scope

- **All Phase 2 work** (`documentation_tutorials.md`): the curated example dataset, any end-to-end pipeline execution, and the T1–T4 tutorial **content**. Phase 1 creates only the landing page and empty slots.
- Native `.ipynb` rendering in Zensical (tracked upstream; the nbconvert stopgap stands in until then).
- mkdocstrings cross-references / backlinks (not yet in Zensical's preliminary handler).
- RTD per-tag versioning (deferred; additive later).
- A Contribution Guidelines page (future work).
- API docs for `utils/` and `cli/` modules (CLI is covered by the CLI tutorials in Phase 2).
- Any CI/pre-commit automation of the docs build.

## Further Notes

- **Phase 2 is not yet grilled.** Only Phase 1 has been designed and verified. Spec'ing the tutorials requires a separate grilling pass first.
- **Config precedence:** where the older `documentation_site.md` and this spec disagree (the `.readthedocs.yaml` sketch, the `docs/requirements.txt` approach, `mkdocs.yml` naming), **this spec wins** — it carries the post-grilling decisions and the verified spike config.
- **uv only:** every dependency operation uses `uv add` / `uv sync`, never pip (repo guideline + active hook). The one `pip install uv` line in the RTD commands is RTD bootstrapping itself, which is expected.
- The `__all__` propose→approve step is a genuine checkpoint inside implementation — the agent must pause for user approval of the Essential set per module before wiring it in.
- **RTD versions:** pinned to newest stable — `os: ubuntu-26.04`, `python: "3.14"` (both RTD-supported; latest stable Python is 3.14.7). Pinned rather than the `ubuntu-lts-latest` / `latest` aliases, for reproducible builds. The docs build is version-insensitive (zensical + static griffe + nbconvert, no heavy deps), and the toolchain was **re-verified green on stable Python 3.14.6** — `zensical build` clean, `workflow/` API rendered with pandas/torch absent. Note the original spike ran on 3.12; the implementing agent should confirm the *full* site still builds green on 3.14 before closing the scaffold ticket. Package floor is `requires-python = ">=3.12"`, so 3.14 is in range.
