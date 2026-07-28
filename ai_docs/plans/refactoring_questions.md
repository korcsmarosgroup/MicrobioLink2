# Refactoring Plan — Grilling Q&A

Questions and answers resolved while stress-testing `ai_docs/plans/microbiolink_refactoring.md`,
before implementation begins.

## Summary of decisions (Q1–Q11)

1. **Package layout**: single flat top-level package `microbiolink/` (no `src/` layout),
   `pyproject.toml` with `hatchling`, mirroring beta's packaging setup.
2. **CLI separation**: all argparse centralized in the cli layer; core modules take plain
   Python arguments only, never import `argparse`.
3. **Module 2 origin**: human-side membrane filter ports from `workflow/get_human_fasta.py`
   `get_proteins()`; bacterial-side ports from `microbiolink_api/microbiome.py`
   `filter_bacterial_domain_table_by_location()` (corrected after initially missing it).
4. **Non-zero expression filter**: kept (not dropped), relocated to Module 3, as a safety net
   for when Module 2 is skipped.
5. **ID translation helper**: shared internal utility, single implementation.
6. **UniProt fetch layer**: shared low-level "fetch by IDs + fields" utility underlies both
   Module 3 (fasta) and Module 4 (domains).
7. **Module 3 & 4 file structure**: one file per capability, single function dispatching on
   `id_type`/`species` — no separate human/bacterial files (reversed from initial recommendation
   after confirming the underlying UniProt calls are identical).
8. **Session scope**: cross-cutting concerns only; Modules 5–10 internals deferred to
   per-module plans.
9. **Source-of-truth precedence**: flat `microbiolink/` wins on overlapping functionality;
   `microbiolink_api` supplies logic only where the flat layout has no equivalent. Corollary:
   DDI/DMI resource data (Module 5) follows `microbiolink_api`'s packaged-TSV +
   `importlib.resources` pattern, since the flat layout has no equivalent.
10. **Testing/CI tooling**: deferred entirely; only packaging essentials adopted now.
11. **Migration path**: old `workflow/*.py` scripts deleted as each is replaced (no shims);
    `case_study_input/` → `case_study_output/` used as the regression fixture per module.
    End state is a single `microbiolink/` package — `microbiolink_api` does not survive.

## Reference facts gathered from the repo

- `case-study` branch (current ground truth) has **no package structure at all**: no
  `pyproject.toml`/`setup.py`, just flat scripts in `workflow/`.
- `MicrobioLink-2.1-beta` branch already has a package layout:
  - Flat top-level package `microbiolink/` (no `src/` layout).
  - `pyproject.toml` using `hatchling` as build backend, `name = "microbiolink"`.
  - Single `microbiolink/cli.py` providing console-script entry points
    (`[project.scripts]` in `pyproject.toml`, e.g. `microbiolink-dmi = "microbiolink.cli:dmi"`).
  - Optional dependency extras: `idr`, `workflow`, `docs`, `tests`, `dev`.
  - `CODING_STYLE.md` at repo root (based on PEP8 + Google style guide, with documented
    exceptions) — this is a **different document** from the
    `TobyL98/toby_verse` python-coding-style.md referenced in the refactoring plan. Needs
    reconciling (see open question below).
  - `microbiolink/cli.py` argparse handling is **inconsistent**: some entry points build
    their `argparse.ArgumentParser` directly inside `cli.py` (e.g. `dmi()`), while others
    delegate to a `parse_args()` function still defined inside the resource module itself
    (e.g. `download_bacterial_proteins.py` keeps its own `parse_args`/`main`). This is in
    tension with the refactoring plan's stated intent that the "cli folder contains argparse
    and any cli-specific boilerplate" and core modules stay free of it.
  - Modules 2 (membrane protein filter) and 4 (Pfam domain download) do **not exist** in
    either `case-study` or as literally described in the plan — closest beta equivalents are
    `download_protein_domains.py` / `download_human_domains.py`, which need to be checked
    against the plan's spec rather than assumed equivalent.

---

## Q1 — Package structure / layout

**Question:** The plan calls for "core" (public + private functions, usable as an API) plus
a separate "cli" folder. Since `case-study` currently has no packaging at all, this refactor
also establishes the package for the first time. What layout should it use?

**Recommendation given:** `src/microbiolink/` layout with `pyproject.toml`, `src/microbiolink/core/`
and `src/microbiolink/cli/`.

**Answer:** Do not use the recommendation — instead mirror the package structure already
established on `MicrobioLink-2.1-beta`: flat top-level `microbiolink/` package (no `src/`),
`pyproject.toml` with `hatchling`, console-script entry points wired through `microbiolink/cli.py`.

**Follow-up still open:** the beta branch's `cli.py` does not cleanly separate all argparse
boilerplate from core logic (some modules keep their own `parse_args`). Need to decide whether
the refactor fixes this inconsistency (all argparse moves out of core, into cli) or preserves
the beta precedent as-is.

---

## Q2 — CLI / argparse separation

**Question:** Beta's `cli.py` inconsistently mixes argparse placement (some in `cli.py`,
some still in the core module as `parse_args()`). The refactoring plan's stated goal is that
core modules should be usable as a package API without the CLI at all. Fix this during the
refactor, or preserve beta's precedent as-is?

**Recommendation given:** Fix it — strip all argparse out of core modules, centralize entirely
in the cli layer, since every module is already being touched for this refactor anyway.

**Answer:** Confirmed — fix it. Argparse is centralized in the cli layer; core modules take
plain Python arguments only, no `parse_args`/`argparse` imports in core.

---

## Additional fact found: Module 2 (membrane protein filter) has no existing implementation

Searched both `case-study` and `MicrobioLink-2.1-beta` for a membrane/secreted protein filter
(OmniPath InterCell for human, UniProt outer-membrane/plasma-membrane for bacteria) — **no such
module exists in either branch.** Closest beta files are `download_protein_domains.py` /
`download_human_domains.py`, which do Pfam domain downloads (module 4 territory), not membrane
filtering. This means Module 2 is net-new development, not a refactor/port, despite the plan
being framed as a reformatting exercise.

Also noted: beta branch has two near-duplicate files, `z-score_filter_terminal.py` and
`z_score_filter_terminal.py` — looks like a stale leftover from a rename, flagging in case it's
relevant when porting Module 1.

---

## Q3 — Module 2 (membrane protein filter): source of the logic

**Question:** Is Module 2 genuinely new development, or does the logic already exist
somewhere not yet found?

**Answer:** The human side already exists — in `workflow/get_human_fasta.py`, function
`get_proteins()` (lines 9–89 of `case-study`). It calls `omnipath.requests.Intercell.get(...)`
with a `location_filter_list` (options seen in the argparse help: `plasma_membrane_transmembrane`,
`plasma_membrane_peripheral`, `secreted`) — this is the Module 2 human/OmniPath-InterCell logic.

**Follow-up fact found while verifying:** `get_proteins()` is heavily entangled — in one
function it does all of:
1. A simple non-zero expression filter (`expression != 0.0`) — resembles but is **not** the
   Module 1 z-score filter (different logic, different purpose, no z-score or NaN involved).
2. The Module 2 membrane/location filter via OmniPath Intercell (human only).
3. Gene-symbol → UniProt ID translation (`translate_symbol_to_uniprot`, via `mygene`).
4. Module 3's fasta fetching (`fetch_protein_sequences`, `main()`).

Also checked `workflow/download_bacterial_proteins.py` (the closest bacterial-side candidate):
it downloads Pfam + gene names by UniProt/proteome ID but does **no** membrane/secreted
filtering. So the bacterial half of Module 2 (UniProt outer-membrane/plasma-membrane/secreted
filtering) has no existing implementation anywhere — that part is genuinely net-new and needs
to be designed (e.g. via UniProt subcellular-location/keyword search), unlike the human half.

---

## Q4 — Where does the non-zero expression filter belong?

**Question:** `get_proteins()` currently drops zero-expression genes as a side effect before
membrane filtering / fasta fetching. Module 1 already does formal z-score-based filtering
(different logic — NaN below cutoff, not a non-zero check). Should the non-zero filter be
dropped as redundant, kept as a Module 2/3 input pre-processing step, or something else?

**Recommendation given:** Drop it — rely solely on Module 1's output so there's only one place
expression filtering happens.

**Answer:** Keep it, but move it — it's still required for the case where Module 2 (membrane
filter) is skipped entirely (i.e. it's not purely redundant with Module 1, it's a safety net
for a pipeline path that bypasses Module 2). It should live in **Module 3** (fasta download),
not Module 2, so gene lists always get non-zero-expression filtering before fasta fetching
regardless of whether Module 2 ran.

---

## Q5 — Shared gene-symbol → UniProt translation helper

**Question:** Both Module 2 and Module 3 need `mygene`-based gene-symbol → UniProt translation
(currently `translate_symbol_to_uniprot` inside `get_human_fasta.py`). Where should it live
once split apart — shared internal utility, duplicated per-module, or folded into Module 3 with
Module 2 calling into it?

**Recommendation given:** Shared internal utility module, single implementation, both Module 2
and Module 3 import it.

**Answer:** Confirmed — shared internal utility module.

---

## Q6 — Shared low-level UniProt fetch utility for Modules 3 & 4

**Question:** Beta already exposes `download_protein_list_with_fields()` / `download_proteome_with_fields()` parameterized by `fields`, reused by both human and bacterial domain downloads. Should the refactor follow this precedent — one shared "fetch UniProt record(s) for these IDs with these fields" utility used by both Module 3 (fasta) and Module 4 (domains) — or should each module make fully independent API calls?

**Recommendation given:** Shared low-level utility, following beta's existing precedent.

**Answer:** Confirmed — shared internal utility module.

---

## Q7 — Module 3 & 4: one file per capability, or split by species?

**Question:** Should Module 3 (fasta) and Module 4 (domains) each be a single file with a
function that dispatches on a `species` argument (matching the plan's literal wording), or
follow beta's existing pattern of a shared low-level utility plus separate `*_human_*.py` /
`*_bacterial_*.py` wrapper files?

**Initial recommendation given:** Follow beta's existing split (shared utility + per-species
files) since it's the established precedent.

**User pushback / fact-check requested:** Do human and bacterial fasta/domain downloads
actually operate identically in UniProt? If so, why are separate per-species wrapper files
needed at all?

**Investigation:** Read all of `get_protein_fasta.py`, `get_human_fasta.py`,
`get_bacterial_fasta.py`, `download_protein_domains.py`, `download_human_domains.py` in full.
Confirmed:
- The actual UniProt fetch call is **identical** for both species in both fasta and domain
  cases — `fetch_fasta_sequences()` and `download_protein_list_with_fields()` are literally the
  same functions, imported and reused as-is, no per-species branching inside them.
- What differs between the "human" and "bacterial" files is **not** the UniProt call, it's:
  (1) input resolution — human wrapper optionally reads a gene-expression file and translates
  gene symbols → UniProt IDs via `mygene`; bacterial wrapper reads a plain pre-resolved ID list;
  (2) bacterial wrapper additionally supports proteome-ID bulk fetch (`fetch_proteome_fasta`),
  a different UniProt query shape not wired up for human in this codebase but not inherently
  species-specific; (3) each file carries its own argparse/CLI flags — moot now that Q2 already
  centralizes all argparse into the cli layer.

**Recommendation reversed to:** one file per capability (one for Module 3, one for Module 4),
each exposing a single function taking `id_type` (uniprot list / proteome ID / gene symbols)
and `species` (human / microbial), branching internally only where real logic differs
(gene-symbol translation when human, proteome bulk-fetch when id_type is proteome), built on
top of the shared low-level fetch utility from Q6.

**Answer:** Confirmed — one file per capability (Module 3, Module 4), single function per file
dispatching on `id_type` and `species`, no separate human/bacterial files.

---

## Q8 — Scope checkpoint

**Question:** Keep drilling at this depth through Modules 5–10 now, or wrap up remaining
cross-cutting concerns only and leave module-specific (5–10) detail to the per-module plans
the top-level document already promises?

**Recommendation given:** Wrap up cross-cutting concerns (coding-style reconciliation, DDI/DMI
data storage location, backward-compatibility expectations, testing/CI scope), leave Modules
5–10 internals to their own per-module plans.

**Answer:** Confirmed — (b). Wrapping up cross-cutting concerns only.

---

## Correction: `microbiolink_api` exists and overturns part of the Q1 and Q3 answers

While investigating where DDI/DMI resource data should live (the next cross-cutting question),
found that beta has a **second package**, `microbiolink_api/`, not examined for Q1–Q7. It sits
alongside the flat `microbiolink/` scripts and looks like an in-progress migration toward a
core/cli split:

- `microbiolink_api/expression.py` — `filter_counts_by_zscore`, `filter_count_matrix_file`
  (Module 1 logic, dataclass/function style).
- `microbiolink_api/microbiome.py` — **`filter_bacterial_domain_table_by_location()`**, which
  filters a bacterial domain table by matching `location_filters` substrings against UniProt's
  `Subcellular location [CC]` column. **This directly contradicts the Q3 answer** — the
  bacterial half of Module 2 is *not* net-new, it already exists here. (`microbiome.py` itself
  still imports its low-level fetch functions from the old flat `microbiolink.download_bacterial_proteins`,
  so it's a partial migration, not a clean break.)
- `microbiolink_api/ddi.py`, `dmi.py` — Modules 5 & 6 as frozen dataclasses (`DomainDomainInteraction`,
  `DomainMotifInteraction`, etc.) plus functions, with **resource bundles loaded via
  `importlib.resources` from packaged TSV files in `microbiolink_api/resources/`**
  (`domine_v2_all_pfam_pairs.tsv`, `3did_dmi_classes.tsv`, `elm_classes.tsv`, etc.) — this is a
  working answer to Module 5's "resources should be stored within a data folder."
- `microbiolink_api/workflows.py` — higher-level orchestration (`run_dmi_workflow`) chaining
  filtering → domain fetch → DMI prediction.
- `microbiolink/DDI.py` (flat package, read earlier for Q6) already imports from
  `microbiolink_api.ddi` — a sign of the migration in progress that wasn't followed up on at
  the time.

## Q9 — Precedence rule between `microbiolink` (flat) and `microbiolink_api`

**Question:** Given `microbiolink_api` already covers some of the same ground as the flat
`microbiolink/` package (plus genuinely new functionality like the bacterial location filter),
which package should the refactor draw from when both have an answer, and when only one does?

**Answer:** The flat `microbiolink/` layout is always used first where functionality is the
same between the two (consistent with the plan's existing case-study-over-beta precedence
rule, extended one level down to these two beta-side packages). `microbiolink_api` is drawn on
specifically where it has functionality with no equivalent in the flat layout — the concrete
example given: `filter_bacterial_domain_table_by_location()` for Module 2's bacterial side.

**Corollary (inferred, flagged for confirmation rather than asked outright):** By the same
rule, DDI/DMI resource storage (Module 5) has no equivalent at all in the flat `microbiolink/`
package — only `microbiolink_api/resources/` implements it, via packaged TSVs loaded through
`importlib.resources`. Following Q9's precedence rule, this packaged-resource pattern should be
adopted for Module 5/6's data folder rather than inventing a new user-supplied-path convention.
This will be carried into the Module 5/6 implementation plan; flag now in case that inference is
wrong.

---

## Resolved fact: coding style guides are consistent, no conflict

Fetched `https://raw.githubusercontent.com/TobyL98/toby_verse/main/guidelines/python-coding-style.md`
(referenced by the plan) and compared against beta's `CODING_STYLE.md`. Same rules: `PascalCase`
classes / `snake_case` functions with resource names as single words, two blank lines between
classes/methods, argument-list line-breaking with trailing commas, Google/Napoleon-style
docstrings, lazy imports for heavy dependencies inside functions rather than at module level
(confirmed in beta code — `import omnipath as op` and `from mygene import MyGeneInfo` both
appear inside function bodies, not at module top). No reconciliation needed — treat as one
style guide, not two competing ones.

---

## Q10 — Testing & CI tooling scope

**Question:** Beta ships `pytest`/`ruff`/`coverage`/`pre-commit`/`bumpversion`/`.cruft.json`/
GitHub Actions CI, none of which exist on `case-study`. Does "package created the same way as
beta" extend to this tooling now, or is it deferred?

**Recommendation given:** Adopt only packaging essentials (pyproject.toml, build backend, entry
points) now; defer CI/pre-commit/versioning tooling to separate later work.

**Answer:** Confirmed — (b). Packaging essentials only for this refactor; CI/tooling deferred.

---

## Q11 — Migration path for old `workflow/` scripts and case-study validation

**Question:** Two related items: (1) once a module's logic moves into the new package, should
the old `workflow/*.py` script be deleted (clean break) or kept temporarily delegating to the
new package? (2) should `case_study_input/` → `case_study_output/` be used as the regression
fixture while porting each module (run old vs. new pipeline on the same input, diff output)?

**Recommendation given:** Delete each old script as its replacement lands (no backward-compat
shims); yes, use `case_study_input/`/`case_study_output/` as the regression fixture per module.

**Answer:** Confirmed — both. (1) Delete old `workflow/*.py` scripts as each is replaced, no
delegating shims kept around. (2) Use `case_study_input/` → `case_study_output/` as the
regression check while porting each module.

---

## Corollary (inferred from Q1 + Q9 + plan's "Overall Layout" section, not yet asked outright)

The plan's "Overall Layout" section describes exactly **one** package with a core (public +
private functions) and a cli layer inside it — not two top-level packages. Combined with Q1
(flat `microbiolink/` layout, mirroring beta) and Q9 (flat layout wins on overlap,
`microbiolink_api` supplies logic only where the flat layout has no equivalent), the end state
should be a **single** package named `microbiolink/`, organized as core + cli, with source code
selectively pulled in from both beta packages per the Q9 precedence rule — `microbiolink_api`
does not survive as a separate installable package in the refactored result, it's a source of
logic to port in, same as the flat `microbiolink/` scripts are. Flagging this inference for
confirmation since it wasn't asked explicitly.

**Answer:** Confirmed — there should be no `microbiolink_api` package after the refactor. It is
a source to port logic from, not a surviving component.
