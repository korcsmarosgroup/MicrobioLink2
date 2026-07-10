# Split the download/fetch script cluster into core primitives vs thin CLI wrappers

All code changes in this plan follow [CODING_STYLE.md](../CODING_STYLE.md).

**Status:** planned, not yet executed.

## Context

Deferred from `plans/collapse_dmi_reverse_into_core.md` as "out of scope."
Candidates #4 and #5 from the 2026-07-10 architecture review of
`microbiolink/`.

Six CLI scripts mix real, reusable UniProt-fetch/gene-identifier logic with
CLI-only `parse_args`/`main` framing, the same shallow-module shape
`DMI.py`/`reverse_DMI.py` had before the previous plan:

- `download_protein_domains.py`, `download_bacterial_proteins.py`,
  `get_protein_fasta.py`, `get_bacterial_fasta.py`,
  `download_human_domains.py`, `get_human_fasta.py`.

`microbiolink/core/microbiome.py` already imports from
`microbiolink.download_bacterial_proteins` — a CLI-layer file — via a
roundabout path (`download_bacterial_proteins.py` is itself just a
re-export shim over `download_protein_domains.py`, its real canonical
location). This violates the CLI-wraps-core principle (see memory
`feedback_cli_wraps_api.md`): core is reaching into a CLI file instead of
the reverse.

Investigation found the UniProt-fetch primitives (`read_ids`,
`build_uniprot_accession_query`, `build_uniprot_stream_url`,
`download_proteome_with_fields`, `download_protein_list_with_fields`,
`fetch_fasta_sequences`, `fetch_proteome_fasta`) are organism-agnostic —
used by both bacterial-side and human-side scripts — so they don't belong
in the bacterial/microbiome-specific `core/microbiome.py`.

**Correction to the original architecture review's framing:**
`core/expression.py` is a different, unrelated concept (KDE-based z-score
threshold filtering of a whole count matrix, replacing
`z_score_filter_terminal.py`'s pipeline stage) from
`download_human_domains.py`'s `_is_expressed` (a single-value
non-NaN/non-zero check). The real duplication is narrower than originally
scoped.

## Confirmed decisions

- **New module `microbiolink/core/uniprot.py`:** houses the
  organism-agnostic UniProt REST client primitives — `read_ids`,
  `build_uniprot_accession_query`, `build_uniprot_stream_url`,
  `download_proteome_with_fields`, `download_protein_list_with_fields`,
  `download_proteome`, `download_protein_list`, `fetch_fasta_sequences`,
  `fetch_proteome_fasta`, and the constants `UNIPROT_STREAM_BASE_URL`,
  `DEFAULT_UNIPROT_FIELDS`, `UNIPROT_BATCH_SIZE`, `UNIPROT_FASTA_URL`,
  `DEFAULT_BATCH_SIZE`, `UNIPROT_FASTA_PROTEOME_URL`. This resolves the
  `read_ids` duplication automatically: `download_protein_domains.py`'s and
  `get_protein_fasta.py`'s byte-identical copies both collapse into the one
  function living here.
- **New module `microbiolink/core/human_domains.py`:** houses the human
  gene-identifier-resolution workflow — `_is_expressed`,
  `read_expressed_genes`, `translate_symbol_to_uniprot`,
  `_extract_swissprot_ids`, the reconciled `get_proteins`, and
  `fetch_protein_sequences`. Named after `download_human_domains.py` since
  most of this logic originates there and it's one real workflow:
  `read_expressed_genes` → `translate_symbol_to_uniprot` →
  `_extract_swissprot_ids` → `download_protein_list_with_fields`.
- **`get_proteins()` reconciliation (three duplicates found and fixed):**
  1. Its 4 branches each independently re-inlined the
     non-NaN/non-zero expression check instead of calling `_is_expressed`.
  2. Its genesymbol branch (no location filter) re-implemented
     `_extract_swissprot_ids` inline, less cleanly (checking
     `for ids in protein: if ids == 'Swiss-Prot'` instead of
     `if 'Swiss-Prot' not in entry: continue`).
  3. A dead-code line: `proteins.extend(translation_dict.values())` is
     immediately overwritten by `proteins = []` two lines later — never
     observable in the output.

  Fix, by branch:
  - **uniprot, no location filter:** replace the file-parsing loop with a
    direct call to `read_expressed_genes(gene_expression_file, sep)` — for
    a single-expression-column file this is behavior-identical to the old
    loop, and picks up `read_expressed_genes`'s `encoding='utf-8-sig'`
    handling that `get_proteins()` lacked (a latent BOM-handling gap, not
    just deduplication).
  - **genesymbol, no location filter:** same treatment —
    `symbols = read_expressed_genes(...)`, then
    `translation_dict = translate_symbol_to_uniprot(symbols)`, then
    `proteins = _extract_swissprot_ids(translation_dict)` — fixing the
    dead code and the duplicate Swiss-Prot-extraction logic at the same
    time.
  - **Both location-filtered branches (genesymbol and uniprot):** keep
    their own row-by-row loop — they need per-row access to write
    `location_filtered_genes.csv` (the original raw line, for genes
    passing both the expression check and the omnipath location-table
    lookup), which `read_expressed_genes()` can't provide since it only
    returns gene IDs. Only the expression-check predicate is shared:
    replace the inline NaN/zero check with a call to `_is_expressed`.
- **`--id_type` casing standardized on lowercase `'uniprot'`:**
  `download_bacterial_proteins.py`'s `choices=['Uniprot', 'UP']` becomes
  `choices=['uniprot', 'UP']`, matching `download_protein_domains.py`/
  `get_bacterial_fasta.py` (already lowercase) and the lowercase
  `'genesymbol'`/`'uniprot'` convention used elsewhere. Breaking change
  only for existing callers of `download_bacterial_proteins` passing
  `--id_type Uniprot` with a capital U.
- **`download_bacterial_proteins.py`'s own `download_proteome`/
  `download_protein_list` wrappers deleted:** they added no behavior over
  `core.uniprot`'s `download_proteome_with_fields`/
  `download_protein_list_with_fields`; `main()` calls the core functions
  directly.
- **Clean break on the existing `test_import_compatibility.py` promise:**
  this compat layer (re-exporting `download_protein_domains.py`'s names
  from `download_bacterial_proteins.py`, and `get_protein_fasta.py`'s
  `fetch_fasta_sequences`-derived `fetch_protein_sequences` from
  `get_human_fasta.py`) was manufactured by a previous internal refactor
  (commit `4d98d6a`), not a promise made to any real external consumer.
  Given the package is pre-alpha, this is a clean break too, consistent
  with the `microbiolink_api` decision: delete the re-exports;
  `tests/test_import_compatibility.py` is deleted (nothing it protects
  remains true once these names live in core).
- **All six files become thin CLI shims**, matching `DMI.py`'s shape
  exactly: only `parse_args()` (extracted for `get_human_fasta.py`, which
  currently parses inline in `main()`) and `main()`, calling
  `core.uniprot`/`core.human_domains` functions directly.

## File mapping

<!-- markdownlint-disable MD013 -->
| From | To |
| --- | --- |
| `download_protein_domains.py`: `read_ids`, `build_uniprot_accession_query`, `build_uniprot_stream_url`, `download_proteome_with_fields`, `download_protein_list_with_fields`, `download_proteome`, `download_protein_list`, `UNIPROT_STREAM_BASE_URL`, `DEFAULT_UNIPROT_FIELDS`, `UNIPROT_BATCH_SIZE` | `microbiolink/core/uniprot.py` |
| `get_protein_fasta.py`: `fetch_fasta_sequences`, `DEFAULT_BATCH_SIZE`, `UNIPROT_FASTA_URL`, `read_ids` (duplicate, dropped) | `microbiolink/core/uniprot.py` |
| `get_bacterial_fasta.py`: `fetch_proteome_fasta`, `UNIPROT_FASTA_PROTEOME_URL` | `microbiolink/core/uniprot.py` |
| `download_human_domains.py`: `_is_expressed`, `read_expressed_genes`, `_extract_swissprot_ids` | `microbiolink/core/human_domains.py` |
| `get_human_fasta.py`: `get_proteins` (reconciled), `translate_symbol_to_uniprot`, `fetch_protein_sequences` | `microbiolink/core/human_domains.py` |
<!-- markdownlint-enable MD013 -->

All six `microbiolink/*.py` files keep only `parse_args()`/`main()`.

## Test strategy

- **New `tests/test_core_uniprot.py`:** primitive-level unit tests migrated
  from `test_download_protein_domains.py`, `test_get_protein_fasta.py`,
  `test_get_bacterial_fasta.py` (the tests currently importing directly
  from the CLI modules), retargeted at `microbiolink.core.uniprot`.
- **New `tests/test_core_human_domains.py`:** primitive-level unit tests
  migrated from `test_download_human_domains.py`, plus **new** tests for
  `get_proteins`, `translate_symbol_to_uniprot`, `fetch_protein_sequences`
  — none of which have any test coverage today, since `get_human_fasta.py`
  has no test file at all.
- **Slimmed CLI test files** (keep only `main()`/CLI end-to-end tests):
  `test_download_protein_domains.py`, `test_get_protein_fasta.py`,
  `test_get_bacterial_fasta.py`, `test_download_human_domains.py`.
- **New CLI test files** (mandatory per AGENTS.md, currently missing
  entirely): `tests/test_download_bacterial_proteins.py`,
  `tests/test_get_human_fasta.py`.
- **Deleted:** `tests/test_import_compatibility.py` (clean break, see
  above).

## Deferred considerations (not solved in this plan)

- A future primitive that yields both the gene ID *and* the raw row (not
  just the ID) could let the two location-filtered branches of
  `get_proteins()` also unify with `read_expressed_genes()`-style reuse,
  removing their own row-by-row loops entirely. Not attempted now — the
  location-table lookup and file-writing side effect make this a bigger
  design change than the rest of this plan.
- `microbiolink/enrichr_id_database_ranking.py` and
  `microbiolink/processing_tiedie_output.py` also use `MyGeneInfo`
  directly (separate from `translate_symbol_to_uniprot`) — out of scope,
  not part of the six-file cluster this plan addresses.

## Related deferred work (separate future plans, not this one)

A broader investigation (prompted mid-grilling, after this plan's design
was otherwise settled) found five more clusters in `microbiolink/` with
the same shallow-CLI-mixed-with-core-logic shape, each unrelated enough to
the other to warrant its own dedicated plan and grilling session rather
than folding into this one:

- **`z_score_filter_terminal.py` → `core/expression.py`.** The highest
  value of the five: `core/expression.py` already appears to be the
  intended rewrite of this exact algorithm (its docstring says so), but
  isn't wired up yet. Not a safe blind swap — real edge-case divergences
  exist (the legacy CLI crashes on non-numeric cells, crashes on
  zero-variance columns, silently NaNs out a whole column when nothing
  exceeds the KDE mode, and divides by zero in a case `core/expression.py`
  guards against), plus the known `--zscore` `type=int`/`required=True`
  bug.
- **`idr_prediction.py` / `idr_prediction_score.py` / `AIUPred.py` →
  `core/`.** `idr_prediction.py` and `idr_prediction_score.py` are still
  byte-identical duplicate files. Found a live scientific bug in
  `motif_selection()`: `anchor_total` is initialized but never
  incremented, so `combined_score` silently reduces to
  `iupred_average / 2`. All three share the unconditional, untested,
  destructive `delete_files_in_folder` pattern already flagged in the
  original architecture review.
- **`motif_monte_carlo_filter.py` → `core/`.** Mostly a relocation: this
  file is already core-shaped (in-memory `dict[str, str]` inputs, already
  imports from `core.dmi`) and just needs moving. Notable: it isn't wired
  to any CLI entrypoint today — an orphaned standalone script.
- **`tiedie_input_processing.py` / `processing_tiedie_output.py` →
  `core/`.** Real pandas network-transformation logic, currently
  interleaves file I/O directly inside the computation functions.
- **`enrichr_id_database_ranking.py`** — needs a design decision before a
  plan can even be written: gene-list/ID-translation helpers are clean
  core candidates, but the enrichment-filtering logic is tangled inline
  inside `main()` together with matplotlib plotting, which doesn't
  cleanly fit the dataclass/DataFrame shape `core/ddi.py`/`core/dmi.py`
  use.

## Deviations from the plan

Minimal — execution matched the plan closely. Notable details worth
recording:

1. `download_protein_domains.py` and `download_bacterial_proteins.py`
   call `download_proteome_with_fields`/`download_protein_list_with_fields`
   directly rather than going through `core.uniprot`'s
   `download_proteome`/`download_protein_list` thin wrappers — the same
   simplification the plan already specified for
   `download_bacterial_proteins.py`'s own (now-deleted) duplicate
   wrappers, applied consistently to both files' internal call sites. The
   `download_proteome`/`download_protein_list` functions still exist in
   `core.uniprot` as a public-API convenience, just aren't the CLI shims'
   internal call path.
2. `get_human_fasta.py`'s `main()` signature changed from `main()` (no
   args, implicit `sys.argv` parsing, no return value) to
   `main(argv: list[str]) -> int`, matching its five siblings in this
   cluster — this was implied by "extract `parse_args()` from inline
   `main()` parsing" but not spelled out as a signature change.
   `cli.py:get_human_fasta()` updated accordingly (`return
   _as_exit_code(module.main(sys.argv[1:]))`, matching
   `download_bacterial_proteins()`/`get_protein_fasta()`'s existing
   pattern) instead of the old `module.main(); return 0`.
3. `tests/test_package_smoke.py`'s `test_cli_download_bacterial_proteins_end_to_end`
   needed fixing, not just the `--id_type Uniprot` → `--id_type uniprot`
   casing update anticipated by the plan: it patched
   `helper_module.requests.get` where `helper_module` was
   `microbiolink.download_protein_domains` — that module no longer
   imports `requests` at all, so the patch target moved to
   `microbiolink.core.uniprot`.
4. Ran `uvx ruff check` (no local `ruff` install) across all touched
   files. Found the same tension already resolved once in Commit 2 of
   `plans/collapse_dmi_reverse_into_core.md`: ruff's isort would combine
   multi-name imports into one parenthesized statement and reorder
   `import x` statements by line length rather than alphabetically. This
   contradicts the one-name-per-line, alphabetical convention used
   throughout `microbiolink/core/*.py` (confirmed the inconsistency is
   pre-existing and not introduced here: `uvx ruff check --diff` on the
   *untouched* `microbiolink/core/dmi.py` proposes the same
   reformatting). Left as-is, consistent with the earlier decision — not
   a new issue, and `AGENTS.md` only mandates markdownlint compliance, not
   ruff compliance. `ANN201`/`ANN001` (missing type annotations on test
   functions) findings are also pre-existing convention across the whole
   test suite (confirmed on untouched `tests/test_reverse_DMI.py` too) —
   not fixed, for the same reason.
5. Added tests beyond the plan's minimum for functions that had zero
   direct coverage even before this migration: `read_ids` (parsing,
   custom column, out-of-range column, UTF-8 BOM handling) and
   `download_proteome`/`download_protein_list` (the thin wrappers) in
   `tests/test_core_uniprot.py`.

## Final outcome

- Two new core modules: `microbiolink/core/uniprot.py` (9 functions/6
  constants) and `microbiolink/core/human_domains.py` (6 functions,
  `get_proteins` reconciled per the plan's three-duplicate fix).
- All six CLI files rewritten as thin shims (`parse_args()`/`main()`
  only); `core/microbiome.py` repointed from
  `microbiolink.download_bacterial_proteins` to `microbiolink.core.uniprot`
  — the CLI-reaches-into-core violation this plan set out to fix is
  resolved.
- `--id_type` casing standardized on lowercase `'uniprot'` in
  `download_bacterial_proteins.py`.
- `download_bacterial_proteins.py`'s redundant `download_proteome`/
  `download_protein_list` wrappers deleted.
- `tests/test_import_compatibility.py` deleted (clean break, as decided).
- `microbiolink/core/__init__.py` exports updated: 9 new functions from
  `core.uniprot`, 4 new functions from `core.human_domains` (constants
  not promoted to the top-level package, matching existing convention —
  CLI shims import constants directly from the submodule).
- Tests: new `tests/test_core_uniprot.py` (17 tests) and
  `tests/test_core_human_domains.py` (16 tests); 4 existing CLI test
  files slimmed to CLI-only end-to-end tests; 2 new CLI test files
  (`test_download_bacterial_proteins.py`, `test_get_human_fasta.py`) for
  the two scripts that had zero coverage before this plan;
  `test_package_smoke.py`'s bacterial-proteins smoke test fixed for the
  new import path and casing.
- Full test suite: **75/75 passed** (65 baseline − 26 removed with the
  old primitive-level/compat tests + 36 added across new and new CLI test
  files).
- Manual smoke checks: all six CLI console scripts (`microbiolink-download-bacterial-proteins`,
  `microbiolink-download-protein-domains`, `microbiolink-get-protein-fasta`,
  `microbiolink-get-bacterial-fasta`, `microbiolink-download-human-domains`,
  `microbiolink-get-human-fasta`) respond correctly to `--help` via the
  installed entry points.
- Not yet committed to git as of writing this section.
