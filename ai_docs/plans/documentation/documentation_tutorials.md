# Plan: MicrobioLink Documentation — Tutorials & Example Data (Phase 2)

Status: proposed
Date: 2026-08-19
Parent: `documentation_build.md` (umbrella). Sibling: `documentation_site.md`.

## Scope

The curated example dataset and all four tutorials. This is the execution-heavy phase: it runs the
MicrobioLink pipeline end-to-end (seeded) once and produces every tutorial artifact from those runs.
Fills the nav slots created by `documentation_site.md`.

- **T1** — API, comprehensive (`.ipynb`, frozen): all 10 modules, forward + reverse.
- **T2** — API, basic (`.ipynb`, frozen): forward-only subset.
- **T3** — CLI, comprehensive (`.md` + committed snippets): all 10 modules, forward + reverse.
- **T4** — CLI, basic (`.md` + committed snippets): forward-only subset.

## Dependencies (must be true before authoring outputs)

1. **Stabilized refactored API** — the `workflow/*` public functions and the `microbiolink-*` console
   scripts are settled, so committed outputs will not immediately drift.
2. **Site plan slots + stopgap wiring** (`documentation_site.md`): the four nav slots and the
   nbconvert `pre_build` step exist.
3. **The curated dataset** (built in Step 1 below).
4. **Full local environment** for the seeded run: base `microbiolink` + `[idr]` + `[tiedie]` extras +
   AIUPRED model data + IUPRED `iupred_data/`.

## Steps

### Step 1 — Curate & verify the example dataset
- Build `docs/tutorials/data/` — a small, clean human + bacterial protein set and a short DEG list
  (a few MB; no `.DS_Store`, no smart-quote filenames, no pre-downloaded bulk).
- **Verify end-to-end** the set yields non-empty DDIs/DMIs/enrichment — no empty teaching tables.
  This is the critical-path artifact; do it first.

### Step 2 — One seeded end-to-end pipeline run
- Run the full pipeline locally on the curated data with **seeded** Monte Carlo (`seed=…`), for both
  forward and reverse where the tutorial calls for it.
- Capture intermediate outputs at each step (the table/plot and "what the data now looks like") — the
  same runs feed both the notebooks and the CLI snippets, so the pipeline is executed once.

### Step 3 — T1 / T2 notebooks (frozen)
- Author `t1_api_full.ipynb` and `t2_api_basic.ipynb` against the stabilized API.
- Each opens with a link to **Get Started** and links to **Concepts** rather than re-explaining.
- T1 covers all 10 modules forward + reverse, flagging optional modules as optional with a note on
  why you might add them. T2 is the forward-only subset: membrane filter (human + bacterial) →
  download fasta → download domains → predict DMIs → IUPRED IDR filter → Monte Carlo → TieDie →
  enrichment, with **enrichment run twice** (initial human proteins, and the final TieDie network).
- Each step displays the key output/data shape. **Commit the `.ipynb` with outputs.**

### Step 4 — T3 / T4 CLI tutorials (Markdown + snippets)
- Author `t3_cli_full.md` and `t4_cli_basic.md`: fenced `microbiolink-*` command blocks with the
  output snippets/screenshots captured in Step 2, saved under `docs/tutorials/assets/`.
- Same open-with-Get-Started / link-to-Concepts pattern; same module coverage as T1/T2 respectively
  (T3 = comprehensive both-directions; T4 = the forward-only basic subset with enrichment twice).

### Step 5 — Verify
- nbconvert renders T1/T2 to Markdown and their committed outputs (tables/plots) display correctly in
  the Zensical build.
- T3/T4 show the committed CLI output; all links to Get Started/Concepts resolve.
- Seeded reruns reproduce the committed notebook outputs (determinism check).

## Acceptance criteria

- Curated dataset produces non-empty DDIs/DMIs/enrichment.
- T1/T2 render committed, seeded outputs (deterministic on rerun); T3/T4 show committed CLI output.
- All four tutorials open with a Get Started link and defer concepts to the Concepts page.
- Enrichment is demonstrated on both the initial human proteins and the final TieDie network in the
  basic tutorials (T2/T4).

## Maintenance note

Committed outputs are refreshed manually when a module API changes (the "light discipline" from the
umbrella's execution decision). Re-run Step 2 and re-commit affected notebooks/snippets.
