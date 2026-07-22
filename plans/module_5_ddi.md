# Module 5 — Domain-Domain Interactions (DDI) — Implementation Plan

Detailed plan for the fifth module described in @plans/microbiolink_refactoring.md, following the
decisions in @plans/refactoring_questions.md and the precedent set by
@plans/module_1_zscore_filter.md through @plans/module_4_download_domains.md. Scope is exactly
Module 5: given Pfam domain dictionaries for a bacterial and a human protein set (Module 4's
output shape), predict domain-domain interactions using the DOMINE and 3did resources.

## Context

Modules 1–4 are implemented in `microbiolink/workflow/` (`zscore_filter.py`, `membrane_filter.py`,
`fasta_download.py`, `domain_download.py`) and `microbiolink/utils/` (`uniprot_client.py`,
`id_resolution.py`), wired into `microbiolink/cli.py`. Module 5 has **no implementation on
`case-study`** — it is spec-only there (see `plans/microbiolink_refactoring.md`'s Module 5
section). The only working implementation anywhere is `microbiolink_api/ddi.py` on
`origin/MicrobioLink-2.1-beta`, which needs to be ported in and simplified to match this repo's
established plain-function/dict/DataFrame style (Modules 1–4 use no dataclasses).

Key facts from beta's `microbiolink_api/ddi.py`:

- Core algorithm is a canonicalized Pfam-pair set-membership lookup: build a set of known
  interacting `(pfam_a, pfam_b)` pairs (sorted, undirected) from packaged TSV resources, then for
  every `(bacterial_pfam, human_pfam)` pair present in the two inputs, check membership and, if
  present, cross-join the bacterial and human proteins carrying those domains.
- Resources are three packaged 2-column, tab-separated, headerless `PfamID<TAB>PfamID` files:
  `pfam_interactions_3did_current.tsv` (3did), `domine_v2_hc_pfam_pairs.tsv` (DOMINE high
  confidence), `domine_v2_all_pfam_pairs.tsv` (DOMINE all). Beta's default merges 3did + DOMINE-all
  and supports a `resource_set` selector plus a `resource` provenance column on each output row.
- Resources are loaded via `importlib.resources.files(package).joinpath(filename)`, per file, with
  no shared loader helper — this repo has no packaged-data precedent yet (no `microbiolink/data/`
  or equivalent), so one is introduced here.
- Beta's input shape is **domain-keyed**: `dict[pfam_id -> [protein_id, ...]]`, produced by
  `microbiolink_api/dmi.py`'s `read_protein_domain_table`. **@plans/module_4_download_domains.md
  has been revised to match this direction** — `download_domains()` now returns
  `dict[pfam_id -> [uniprot_id, ...]]` per species (previously protein-keyed), specifically so
  Module 5 can consume it as-is with no invert step anywhere in this module. Note: the
  already-shipped `microbiolink/workflow/domain_download.py` on this branch still returns the old
  protein-keyed shape and needs to be updated per that plan before Module 5 can be built against
  it.

Two scope decisions were made with the user before writing this plan:

1. Default resource set is **3did + DOMINE high confidence** (not DOMINE-all, beta's default) —
   only 2 of beta's 3 TSVs need porting.
2. **No resource-set selector and no custom-resource-file override** — always use the merged
   default, no `resource_set` parameter. Output keeps a **5th `resource` column** (`'3did'`,
   `'DOMINE_hc'`, or `'3did|DOMINE_hc'`) beyond the plan's 4 required columns, since it's cheap to
   retain and useful for downstream filtering.

No `case_study_output` fixture exists for DDI (case-study has no DDI code at all), so per Q11 this
counts as new functionality requiring manual sign-off rather than a diff-based regression check.
`case_study_input/input/bacterial_protein/BT_BEV_domains.tsv` gives a real set of bacterial `Entry`
IDs to re-fetch through `microbiolink-download-domains` (its own legacy Entry/Pfam format doesn't
match Module 4's revised Pfam/Entries CLI output directly, so it's used as an ID source, not read
in as-is); the same CLI generates a matching human domain file.

## Decisions

1. **No dataclasses.** Beta's `DomainDomainInteraction`/`DDIResourceBundle` are dropped. The
   "resource bundle" is a private `dict[tuple[str, str], list[str]]` mapping a canonicalized Pfam
   pair to the source name(s) (`'3did'`/`'DOMINE_hc'`) that contain it — matching Modules 1–4's
   plain dict/DataFrame style.
2. **Core module accepts Module 4's dict shape directly, no reshaping needed.** No separate
   "from_data" vs. file-based split like beta — a single public function,
   `predict_domain_domain_interactions(bacterial_domains, human_domains)`, taking
   `dict[pfam_id -> [uniprot_id, ...]]` for each species (Module 4's revised per-species output,
   used as-is) and returning a `pd.DataFrame`. No file I/O in `workflow/ddi.py`, matching Module
   2/4's precedent of pure-return workflow functions.
3. **New `microbiolink/data/` package for resources**, following the plan's literal "resources
   ... should be stored within a data folder" wording. Contains `__init__.py` (docstring only, the
   `importlib.resources` anchor) plus the two packaged TSVs. Loaded once via a
   `functools.lru_cache`-wrapped private loader, since both files are static and otherwise would be
   re-read from disk on every call.
4. **Only 2 resource files ship**, not beta's 3 — `domine_v2_all_pfam_pairs.tsv` is dropped since
   the merged default no longer includes DOMINE-all (per the resource-set decision above).
5. **CLI reads Module 4's `Pfam`/`Entries` TSV output directly, no shared parser needed.** Since
   Module 4's revised CLI output is already Pfam-keyed (one row per Pfam ID, semicolon-joined
   Entries column), turning it back into `dict[pfam_id -> [uniprot_id, ...]]` is a trivial
   one-to-one read — no fan-out/invert logic, so it doesn't need to reuse
   `domain_download.domain_table_to_mapping` (which parses UniProt's raw per-protein fetch response,
   a different, unrelated table shape). A small private `_read_domain_mapping` helper in `cli.py`
   covers it.
6. **Output columns**: `bacterial_uniprot_id`, `bacterial_pfam_domain`, `human_uniprot_id`,
   `human_pfam_domain`, `resource` — snake_case naming consistent with Module 2's
   `['uniprot_id', 'location_annotation']` output convention, rather than beta's
   `bacterial_protein`/`human_protein`/`bacterial_domain`/`human_domain` naming.
7. **CLI writes CSV** (`result.to_csv(output_path, index=False)`), comma-separated, matching
   Module 2's output convention, rather than beta's `;`-separated default.

## Target implementation

### `microbiolink/data/` (new package, resources only)

```
microbiolink/data/__init__.py                        # docstring only
microbiolink/data/pfam_interactions_3did_current.tsv  # copied from beta via `git show`
microbiolink/data/domine_v2_hc_pfam_pairs.tsv         # copied from beta via `git show`
```

Both copied unmodified from
`origin/MicrobioLink-2.1-beta:microbiolink_api/resources/pfam_interactions_3did_current.tsv` and
`.../domine_v2_hc_pfam_pairs.tsv` (2-column, tab-separated, headerless `PfamID<TAB>PfamID`).

### `microbiolink/workflow/ddi.py` (new, core, no argparse)

```python
def _read_pfam_pair_table(filename) -> set[tuple[str, str]]:
    """Read a headerless two-column Pfam-Pfam interaction TSV into a canonical pair set."""


@functools.lru_cache(maxsize=None)
def _load_ddi_resource_pairs() -> dict[tuple[str, str], list[str]]:
    """Load the packaged 3did and DOMINE high-confidence Pfam pair resources."""


def predict_domain_domain_interactions(
    bacterial_domains: dict[str, list[str]],
    human_domains: dict[str, list[str]],
) -> pd.DataFrame:
    """Predict domain-domain interactions between bacterial and human proteins."""
```

- `_read_pfam_pair_table`: splits each line on `'\t'`, skips lines with fewer than 2 fields,
  canonicalizes via `tuple(sorted((fields[0], fields[1])))`, returns the set.
- `_load_ddi_resource_pairs`: loads both packaged TSVs via
  `importlib.resources.files(data).joinpath(filename)` (imports the sibling `data` package), builds
  `{pair: ['3did']}` / `{pair: ['DOMINE_hc']}` per file, then merges — pairs in both files get
  `['3did', 'DOMINE_hc']`.
- `predict_domain_domain_interactions`: iterates every `(bacterial_pfam, human_pfam)` combination
  from the two input dicts' keys directly (no inversion — both are already Pfam-keyed); canonicalizes
  each pair and looks it up in `_load_ddi_resource_pairs()`; if present, cross-joins the bacterial
  and human proteins for that pair, deduplicating via a `seen` set keyed on
  `(bacterial_protein, human_protein, bacterial_pfam, human_pfam)`, and appends a row with
  `resource = '|'.join(sources)`; returns
  `pd.DataFrame(rows, columns=['bacterial_uniprot_id', 'bacterial_pfam_domain', 'human_uniprot_id', 'human_pfam_domain', 'resource'])`
  (explicit columns so an empty result still has the right schema).

### `microbiolink/cli.py` (extend, argparse only)

- Add `_build_ddi_parser()`: `-b/--bacterial_domain_file` (required),
  `-hu/--human_domain_file` (required), `-o/--output_file` (required).
- Add `_read_domain_mapping(filename) -> dict[str, list[str]]`: `pd.read_csv(filename, sep='\t')`
  then, per row, `{row['Pfam']: [entry for entry in row['Entries'].split(';') if entry]}` — a
  direct read of Module 4's `Pfam`/`Entries` TSV format, already Pfam-keyed.
- Add entry point:
  ```python
  def ddi() -> int:
      """Predict domain-domain interactions between bacterial and human proteins."""
      from .workflow import ddi as ddi_workflow

      args = _build_ddi_parser().parse_args()
      bacterial_domains = _read_domain_mapping(args.bacterial_domain_file)
      human_domains = _read_domain_mapping(args.human_domain_file)

      result = ddi_workflow.predict_domain_domain_interactions(bacterial_domains, human_domains)
      result.to_csv(args.output_file, index=False)
      return 0
  ```
- `pyproject.toml`: add `microbiolink-ddi = "microbiolink.cli:ddi"` under `[project.scripts]`, and
  add `[tool.hatch.build] include = ["microbiolink/data/*.tsv"]` so the two TSVs ship in the wheel
  (mirrors beta's identical `include` line for its own resources directory).

## Migration checklist

### 1. Resource files
- [ ] Create `microbiolink/data/__init__.py` (docstring only).
- [ ] Copy `pfam_interactions_3did_current.tsv` and `domine_v2_hc_pfam_pairs.tsv` from
      `origin/MicrobioLink-2.1-beta:microbiolink_api/resources/` into `microbiolink/data/`.

### 2. Core module
- [ ] Confirm `microbiolink/workflow/domain_download.py` has been updated per
      @plans/module_4_download_domains.md to return the Pfam-keyed shape (prerequisite for this
      module).
- [ ] Create `microbiolink/workflow/ddi.py` with `_read_pfam_pair_table`, `_load_ddi_resource_pairs`,
      `predict_domain_domain_interactions`.

### 3. CLI wiring
- [ ] Add `_build_ddi_parser`, `_read_domain_mapping`, `ddi()` entry point to `cli.py`.
- [ ] Add `microbiolink-ddi = "microbiolink.cli:ddi"` and
      `[tool.hatch.build] include = ["microbiolink/data/*.tsv"]` to `pyproject.toml`.

### 4. Manual sign-off (per Q11 — no case-study ground truth exists for this module)
- [ ] Confirm a pair present in both packaged files (e.g. `PF00001` paired with itself, present in
      both TSVs per sampled rows) round-trips through `_load_ddi_resource_pairs()` with
      `resource == '3did|DOMINE_hc'`.
- [ ] Generate a bacterial domain file by running `microbiolink-download-domains` against the
      `Entry` IDs in `case_study_input/input/bacterial_protein/BT_BEV_domains.tsv` (that fixture is
      the legacy Entry/Pfam protein-keyed format, not directly usable as `--bacterial_domain_file`
      now that Module 4's output format is Pfam-keyed).
- [ ] Generate a human domain file by running `microbiolink-download-domains` against
      `case_study_input/input/human_protein/Enterocyte_Manual/enterocyte_colon_CD_expressed_genes.csv`.
- [ ] Run `microbiolink-ddi` end-to-end on these two files and show the resulting output table to
      the user for confirmation before considering Module 5 done.

## Verification

- `uv run ruff format .`
- `uv run ruff check .`
- `uv run ty check`
- Manual CLI run as in the sign-off checklist above.
