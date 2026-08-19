# Plan: Building the MicrobioLink Documentation (umbrella)

Status: proposed
Date: 2026-08-19

Umbrella plan for the MicrobioLink documentation. The plain-language brief is `new_documentation.md`;
this umbrella holds the cross-cutting decisions, verified facts, information architecture, and the
two-phase split. Detailed steps live in the two child plans:

- **`documentation_site.md`** — the static site (scaffold, config, RTD, Introduction, Get Started,
  Concepts, API Reference, README trim, tutorial nav slots). Ships green **without any tutorial
  existing yet**.
- **`documentation_tutorials.md`** — the curated example dataset, a seeded end-to-end pipeline run,
  and all four tutorials (T1/T2 frozen notebooks, T3/T4 CLI Markdown). Drops into the slots the
  site plan creates.

**Why two children:** the real dependency boundary is *"requires running the pipeline end-to-end on
the curated dataset."* All four tutorials share that heavy prerequisite (T3/T4 need it to capture CLI
output just as T1/T2 do), while the static site needs none of it. Isolating the execution-heavy work
keeps it from blocking the site from going live.

## Decisions (locked)

1. **Generator: stay on Zensical.** Maintained forward path — MkDocs 1.x is EOL and MkDocs 2.0
   removes plugins / breaks Material for MkDocs. Zensical is the sanctioned "drop-in replacement".
2. **Notebooks via a throwaway nbconvert stopgap.** Zensical cannot render `.ipynb` today
   (`mkdocs-jupyter` only backlogged — [#52](https://github.com/zensical/zensical/issues/52),
   [#96](https://github.com/zensical/zensical/issues/96)). A temporary nbconvert pre-build step
   renders committed notebooks to Markdown; **delete when Zensical ships native support.**
3. **mkdocstrings "preliminary" accepted** (v0.0.11; no cross-refs/backlinks yet —
   [#237](https://github.com/zensical/zensical/issues/237)).
4. **Notebook execution: pre-execute locally & commit outputs, seeded.** RTD runs **nbconvert only**.
   Authored against the *stabilized* refactored API; re-run on API change.
5. **Shared "Pipeline Concepts" page + 4 lean tutorials** that link back instead of re-explaining.
6. **CLI tutorials (T3/T4) are Markdown + committed snippets/screenshots.** Only T1/T2 are frozen
   notebooks.
7. **Small curated tutorial dataset** tuned to yield non-empty DDIs/DMIs/enrichment — not the 91 MB
   `case_study_input/`.
8. **API split via `__all__` per module** → mkdocstrings "Essential" vs "Helpers". Zero runtime effect.
9. **Docs Introduction is canonical; README is a pointer.**
10. **Concepts is a top-level nav section**; Introduction links to all four (Get Started, Concepts,
    Tutorials, API).
11. **Single "latest" version on RTD now**; per-tag versioning deferred (additive later).
12. **RTD builds "latest" from `refactoring`** now, repointed to `case-study` (main) on merge.

## Verified facts (de-risk the build)

- **Docs build is light.** Heavy deps (`iupred`, `tiedie`) are imported *lazily inside functions*,
  never at module scope; mkdocstrings/griffe is static (AST). RTD needs only the base `microbiolink`
  package + zensical + mkdocstrings + nbconvert — **no torch, no git extras**.
- **Monte Carlo is seedable** (`seed: int = 0` → `np.random.default_rng`), backing seeded-frozen
  notebooks.
- **Interaction resources (3did, ELM, DOMINE)** are committed `.tsv` in `microbiolink/data/` →
  deterministic.
- **RTD Community limits:** 15 min build, 7 GB RAM, 5 GB disk (soft) — why build-time execution of
  the tutorials is infeasible (hence pre-execute-and-commit).

## Site information architecture

```
Introduction            (canonical overview; links → Get Started, Concepts, Tutorials, API)
Get Started             (install methods; derived from ../package_build.md)
Pipeline Concepts       (the 10 modules, forward/reverse, why optionals — written once)
Tutorials
  ├─ T1  API, comprehensive  (.ipynb, frozen)  all 10 modules, forward + reverse
  ├─ T2  API, basic          (.ipynb, frozen)  forward-only subset
  ├─ T3  CLI, comprehensive  (.md + snippets)  all 10 modules, forward + reverse
  └─ T4  CLI, basic          (.md + snippets)  forward-only subset
API Reference            (per module: Essential vs Helper, from mkdocstrings + __all__)
```

Owner map: the **site plan** owns Introduction, Get Started, Concepts, API Reference, the Tutorials
landing page, and the four nav *slots*. The **tutorials plan** owns the content that fills T1–T4 and
the example dataset.

## Phasing & dependencies

1. **Phase 1 — `documentation_site.md`.** Independent of pipeline execution; can be built and shipped
   green first. Creates the nav slots and the nbconvert stopgap wiring the tutorials will use.
2. **Phase 2 — `documentation_tutorials.md`.** Depends on: (a) the stabilized refactored API, (b) the
   site plan's nav slots + stopgap wiring, (c) the curated dataset it builds first. Runs the seeded
   pipeline once and produces all four tutorials from those runs.

The phases can start in parallel where convenient (e.g. dataset curation early), but Phase 1 does not
block on Phase 2 — that is the point of the split.

## Future work

- Delete the nbconvert stopgap once Zensical renders `.ipynb` natively.
- Turn on RTD per-tag versioning when tagged releases begin; repoint `latest` to `case-study` on merge.
- Add a Contribution Guidelines page.
