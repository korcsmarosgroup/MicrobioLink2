# Module 3 — Downloading Fasta — Implementation Plan

Detailed plan for the third module described in @ai_docs/plans/microbiolink_refactoring.md, following the
decisions in @ai_docs/plans/refactoring_questions.md and the precedent set by
@ai_docs/plans/module_2_membrane_filter.md. Scope is exactly Module 3: resolving a set of human or
microbial protein identifiers (uniprot IDs, a proteome ID, gene symbols, or gene symbols from a
count matrix) to UniProt accessions, downloading their sequences, and writing them to FASTA files.

## Context

Modules 1 (z-score filter) and 2 (membrane protein filter) are already implemented in
`microbiolink/zscore_filter.py`, `microbiolink/membrane_filter.py`, `microbiolink/uniprot_client.py`,
`microbiolink/gene_matrix.py` and wired into `microbiolink/cli.py`. Module 3's input contract is
worded almost identically to Module 2's in the top-level plan, and its output (Module 2's output can
feed directly into Module 3) and the top-level plan's own note for Module 4 ("reuse functions where
possible" — same inputs as Module 3) both call for building this module so its identifier-resolution
logic is cleanly reusable by Module 4 later, not fused into a single fasta-only function.

Two research passes over `case-study` and `MicrobioLink-2.1-beta` confirm:
- **No FASTA capability exists in `case-study` at all** for bacterial proteins; only
  `workflow/get_human_fasta.py`'s `fetch_protein_sequences()` covers human, entangled with the
  Module 2 membrane filter and `mygene`-based gene-symbol translation.
- **Beta already split this cleanly**: `microbiolink/get_protein_fasta.py` (`fetch_fasta_sequences`,
  batch size 100), `microbiolink/get_bacterial_fasta.py` (proteome-fasta fetch), and
  `microbiolink/get_human_fasta.py` (gene-symbol translation via `mygene`, delegating the actual
  fetch to `get_protein_fasta.fetch_fasta_sequences`). None of this lives in `microbiolink_api/` —
  `microbiolink_api/workflows.py` only *consumes* an already-downloaded FASTA file, it never
  produces one, so this module is genuinely new work for the current package, informed by (not
  copied verbatim from) beta's split.
- The current `microbiolink/uniprot_client.py` is TSV-only end-to-end (`build_uniprot_stream_url`
  hardcodes `format=tsv` and always injects a `fields=` clause; `_parse_uniprot_response` assumes
  tab-separated data). It needs new FASTA-format sibling functions, not a repurposing of the
  existing ones.
- No gene-symbol → UniProt translation utility exists yet in the refactored package (Q5's shared
  utility was deferred by Module 2 since Module 2 doesn't need it — OmniPath Intercell already
  carries both `uniprot` and `genesymbol` columns). Module 3 is where this utility actually gets
  built.

## Package layout

Adding Module 3's files as more flat top-level modules alongside Modules 1 & 2 makes the package
hard to read from the outside: four utility-ish modules (`uniprot_client.py`, `gene_matrix.py`,
`id_translation.py`, `id_resolution.py`) all feeding into one `fasta_download.py`, indistinguishable
at a glance from the actual per-module entry points. This plan instead splits `microbiolink/` into
two subpackages, applied retroactively to Modules 1 & 2 as well so the package isn't left
half-migrated:

```
microbiolink/
  __init__.py
  cli.py                      (unchanged location — argparse layer, per Q2)
  workflow/
    __init__.py
    zscore_filter.py           Module 1 — moved from top level, logic unchanged
    membrane_filter.py          Module 2 — moved from top level, logic unchanged
    fasta_download.py           Module 3 — new
  utils/
    __init__.py
    uniprot_client.py           moved from top level, extended with fetch_fasta_sequences
    gene_matrix.py              moved from top level, extended (see decisions below)
```

- `workflow/` — one file per pipeline module, mirroring the top-level plan's own "Module N"
  numbering. Future modules 4+ each add their own file here.
- `utils/` — shared plumbing any workflow file can import without depending on another workflow
  file's internals: `uniprot_client.py`, `gene_matrix.py`, `id_resolution.py`.

## Decisions from grilling

1. **Output shape — confirmed with user**: a single public function taking identifiers for *both*
   species at once (either may be omitted) plus one `output_folder`, writing
   `human_proteins.fasta` and/or `microbial_proteins.fasta` into it. This matches the top-level
   plan's literal wording ("one file for microbe and one for human in a folder") more directly than
   Module 1/2's single-`--output_file`-per-call pattern, and was explicitly chosen over the
   per-species-per-call alternative when asked.
2. **FASTA batch size — confirmed with user**: a new `FASTA_BATCH_SIZE = 100` constant, distinct
   from the existing `UNIPROT_BATCH_SIZE = 1000` used for TSV/tabular fetches. This matches both
   legacy sources (case-study's `fetch_protein_sequences` used `batch_size = 100`; beta's
   `get_protein_fasta.py` used `DEFAULT_BATCH_SIZE = 100`) and is safer given FASTA records are a
   much larger payload per ID than a TSV row.
3. **Human + proteome ID: out of scope**, inherited from Module 2 decision 2 (identical reasoning —
   a full human proteome is not a realistic input; Module 2 already restricted human `id_type` to
   `'uniprot'`/`'genesymbol'` only). Module 3 keeps the same restriction: human accepts `'uniprot'`
   or `'genesymbol'`; microbial accepts `'uniprot'` or `'proteome'`.
4. **Proteome resolution is two-step, per the plan's literal wording**: "For uniprot proteome ID, it
   should get all the uniprot IDs first and then get the fastas." Rather than porting beta's direct
   `fetch_proteome_fasta()` (a `format=fasta&query=proteome:ID` call), this reuses the **existing**
   `uniprot_client.fetch_proteome_table(proteome_id, fields=['accession'])` to resolve the ID list
   (the `accession` field returns under an `Entry` column — confirmed by how
   `membrane_filter.py` already renames `Entry` → `uniprot_id`), then fetches FASTA for that flat ID
   list through the same accession-based path as everything else. This avoids introducing a second,
   never-reused proteome-fasta query shape and sidesteps a paren-encoding inconsistency spotted
   between beta's tabular (`%28%28proteome%3A...%29%29`) and fasta (`%28proteome%3A...%29`) proteome
   queries.
5. **Non-zero expression filter lands in Module 3**, per Q4. A new
   `extract_expressed_gene_symbols_from_count_matrix()` in `gene_matrix.py` (keeping only symbols
   with at least one non-NaN, non-zero value across sample columns) is used by Module 3's
   `--from-count-matrix` CLI path instead of Module 2's plain `extract_gene_symbols_from_count_matrix`
   — Module 2 explicitly keeps the unfiltered extractor (decision 4 keeps this filter out of Module
   2).
6. **Gene-symbol translation is a private helper, not its own file.** Legacy
   `translate_symbol_to_uniprot()` (via `mygene.MyGeneInfo`) has exactly one caller in the new
   design (`resolve_uniprot_ids`, decision 7), so it becomes `_translate_gene_symbols_to_uniprot()`
   inside `utils/id_resolution.py` rather than a separate `id_translation.py` module. It still fixes
   a latent bug in the legacy version: legacy's own inline unpacking of the translation dict (in
   `get_proteins()`'s unfiltered branch) expects a `'Swiss-Prot'` key on values that
   `translate_symbol_to_uniprot` doesn't actually produce that shape for; the *correct* extraction
   logic already exists separately in beta's `download_human_domains.py`
   (`_extract_swissprot_ids`). This is new behavior, not a straight port — flagged for manual
   confirmation per the top-level plan's "if new functionality, get my confirmation" rule.
7. **New shared module for ID resolution, dispatched on `id_type` alone**: `utils/id_resolution.py`.
   `resolve_uniprot_ids(identifiers, id_type) -> list[str]` branches on `id_type` —
   `'uniprot'` passes through, `'genesymbol'` calls `_translate_gene_symbols_to_uniprot()` (always
   human), `'proteome'` calls `uniprot_client.fetch_proteome_table()` (always microbial). No
   `species` parameter: per decision 3, `genesymbol` is only ever human and `proteome` only ever
   microbial, so `id_type` already determines the behavior (unlike `membrane_filter.py`, where
   human vs. bacterial are genuinely different algorithms). Module 4 is expected to reuse this
   directly, swapping in `uniprot_client.fetch_protein_table(uniprot_ids, fields=[pfam fields])`
   for the fasta fetch.
8. **Package restructured into `workflow/` + `utils/` folders** (see "Package layout" above),
   applied retroactively to Modules 1 & 2 as well as Module 3:
   - `workflow/zscore_filter.py`, `workflow/membrane_filter.py` — Modules 1 & 2, moved unchanged.
   - `workflow/fasta_download.py` — Module 3's own dispatch + FASTA writing (new).
   - `utils/uniprot_client.py` — moved; gains FASTA-fetch functions alongside the existing
     tabular ones.
   - `utils/gene_matrix.py` — moved; gains the expressed-only extractor (decision 5) and
     `read_count_matrix` (decision 9).
   - `utils/id_resolution.py` — new, per decisions 6/7.
9. **`read_count_matrix` moves from `zscore_filter.py` into `utils/gene_matrix.py`.** It's
   currently Module 1's function, but Module 2's and Module 3's CLI code both call it directly to
   read a count matrix before extracting gene symbols from it — exactly the cross-module
   reach-through the `workflow`/`utils` split (decision 8) is meant to eliminate. After the move,
   `workflow/zscore_filter.py` itself imports it from `utils/gene_matrix.py` like every other
   caller.
10. **New dependency**: `mygene>=3.2.2` added to `pyproject.toml` (matches beta's exact pin — no
    version pin currently exists anywhere in this repo, so beta's is the only prior art).
11. **CLI helper reuse**: `cli.py`'s `_add_membrane_filter_source_arguments` is renamed to
    `_add_identifier_source_arguments` (behavior unchanged) since Module 3's human-identifier CLI
    args need the exact same four flags (`--input_file`, `--from-count-matrix`, `--sep`,
    `--id_column`) — this is a rename of an already-generic helper to reflect it now has two call
    sites, not a new abstraction.

## Target implementation

### `microbiolink/utils/uniprot_client.py` (moved from top level, then extended)

```python
FASTA_BATCH_SIZE = 100


def build_uniprot_fasta_url(query: str) -> str:
    """Build a UniProt stream endpoint URL for FASTA format (no fields)."""


def _fetch_fasta_batch(identifiers: list[str]) -> list[str]:
    """Fetch one FASTA batch as raw text, retrying by splitting on HTTP 400."""


def fetch_fasta_sequences(identifiers: list[str]) -> str:
    """Fetch FASTA sequence text for a batch of UniProt accessions."""
```

- `build_uniprot_fasta_url`: `f'{UNIPROT_STREAM_BASE_URL}format=fasta&query={query}'` — reuses
  `UNIPROT_STREAM_BASE_URL`, but unlike `build_uniprot_stream_url` takes no `fields` (FASTA format
  doesn't support field selection).
- `_fetch_fasta_batch`: same accession-query building (`build_uniprot_accession_query`, already
  reusable as-is per the earlier research) and the same recursive-split-on-400 retry pattern as
  `_fetch_protein_batch`, but returns raw response text in a list instead of parsed data frames.
- `fetch_fasta_sequences`: batches `identifiers` in groups of `FASTA_BATCH_SIZE` (not
  `UNIPROT_BATCH_SIZE`), concatenates all batch texts, returns one string.
- No proteome-fasta function is added here — proteome resolution goes through the existing
  `fetch_proteome_table` (decision 4), so no new query shape is needed.

### `microbiolink/utils/gene_matrix.py` (moved from top level, then extended)

```python
def read_count_matrix(
    filename: PathLike,
    index_col: int | str = 0,
) -> pd.DataFrame:
    """Read a gene/protein count matrix from disk."""


def extract_gene_symbols_from_count_matrix(count_matrix: pd.DataFrame) -> list[str]:
    """Extract the gene symbol row index from a gene count matrix."""


def extract_expressed_gene_symbols_from_count_matrix(count_matrix: pd.DataFrame) -> list[str]:
    """Extract gene symbols with at least one non-zero, non-NaN expression value."""
```

- `read_count_matrix` is ported unchanged from `zscore_filter.py` (decision 9) — every caller
  (`workflow/zscore_filter.py`, and the CLI resolvers for Modules 2 & 3) now imports it from here
  instead of reaching into Module 1's workflow file.
- `extract_expressed_gene_symbols_from_count_matrix` (new, decision 5):
  `count_matrix.index[count_matrix.notna().any(axis=1) & (count_matrix != 0).any(axis=1)].tolist()`
  — a safety-net filter (Q4), not a replacement for Module 1's formal z-score filter.

### `microbiolink/utils/id_resolution.py` (new, core, no argparse)

```python
def _translate_gene_symbols_to_uniprot(gene_symbols: list[str]) -> dict[str, list[str]]:
    """Translate human gene symbols to UniProt Swiss-Prot accessions via MyGene.info."""


def resolve_uniprot_ids(identifiers: list[str], id_type: str) -> list[str]:
    """Resolve protein identifiers to a flat list of UniProt accessions."""
```

- `_translate_gene_symbols_to_uniprot` (private — decision 6): lazily imports
  `from mygene import MyGeneInfo` inside the function body (matches the established lazy-import
  convention for heavy/optional dependencies, e.g. `omnipath` in `membrane_filter.py`). Calls
  `MyGeneInfo().querymany(gene_symbols, scopes='symbol', fields='uniprot', species='human',
  returnall=True)`, then for each entry in `results['out']` extracts
  `entry['uniprot']['Swiss-Prot']` (normalizing to a list whether MyGene returns a single string or
  a list), skipping entries with no `uniprot`/`Swiss-Prot` mapping. Returns
  `{gene_symbol: [uniprot_id, ...]}`.
- `resolve_uniprot_ids` (public, decision 7) branches on `id_type` alone: `'uniprot'` returns
  `identifiers` unchanged; `'genesymbol'` calls `_translate_gene_symbols_to_uniprot`, flattens and
  dedupes the resulting values; `'proteome'` loops each proteome ID through
  `uniprot_client.fetch_proteome_table(proteome_id, fields=['accession'])['Entry'].tolist()` and
  concatenates. Any other `id_type` raises `ValueError`.

### `microbiolink/workflow/fasta_download.py` (new, core, no argparse)

```python
def download_fasta(
    output_folder: PathLike,
    human_identifiers: list[str] | None = None,
    human_id_type: str | None = None,
    microbial_identifiers: list[str] | None = None,
    microbial_id_type: str | None = None,
) -> dict[str, Path]:
    """Download human and/or microbial protein sequences as FASTA files."""
```

- Raises `ValueError` if both `human_identifiers` and `microbial_identifiers` are `None`.
- Creates `output_folder` if missing (`Path(output_folder).mkdir(parents=True, exist_ok=True)`).
- For each species supplied: calls `id_resolution.resolve_uniprot_ids(...)` then
  `uniprot_client.fetch_fasta_sequences(...)`, writes the text to `human_proteins.fasta` /
  `microbial_proteins.fasta` inside `output_folder`.
- Returns `{'human': Path(...), 'microbial': Path(...)}` for whichever species were requested.

### `microbiolink/cli.py` (extend existing file, argparse only)

- Update every existing lazy import to the new subpackage paths:

  ```python
  # before                                       # after
  from . import zscore_filter                    from .workflow import zscore_filter
  from . import membrane_filter as ...            from .workflow import membrane_filter as ...
  from . import uniprot_client                    from .utils import uniprot_client
  from . import gene_matrix                       from .utils import gene_matrix
  ```

- Rename `_add_membrane_filter_source_arguments` → `_add_identifier_source_arguments` (update the
  one existing call site in `_build_membrane_filter_parser`).
- Add:
  - `_add_human_fasta_arguments(parser)`: `-hi/--human_input_file`, `-hid/--human_id_type`
    (choices `uniprot`/`genesymbol`), `--human-from-count-matrix` (flag), `-hsep/--human_sep`
    (default `,`), `-hcol/--human_id_column` (int, default 1) — all optional, since the human side
    as a whole may be omitted.
  - `_add_microbial_fasta_arguments(parser)`: `-mi/--microbial_input_file`,
    `-mid/--microbial_id_type` (choices `uniprot`/`proteome`), `-msep/--microbial_sep` (default
    `,`), `-mcol/--microbial_id_column` (int, default 1) — no count-matrix option for microbial (no
    realistic bacterial count-matrix use case in this codebase).
  - `_build_fasta_download_parser()`: composes both plus `-o/--output_folder` (required).
  - `_resolve_human_fasta_identifiers(args) -> list[str] | None`: `None` if
    `args.human_input_file` is unset; otherwise `--human-from-count-matrix` routes through
    `gene_matrix.read_count_matrix` + `gene_matrix.extract_expressed_gene_symbols_from_count_matrix`,
    else `uniprot_client.read_ids(...)`.
  - `_resolve_microbial_fasta_identifiers(args) -> list[str] | None`: `None` if
    `args.microbial_input_file` is unset; otherwise `uniprot_client.read_ids(...)`.
  - `download_fasta() -> int`: builds parser, resolves both identifier lists, calls
    `fasta_download.download_fasta(args.output_folder, human_identifiers=..., human_id_type=...,
    microbial_identifiers=..., microbial_id_type=...)`, returns 0.
- `pyproject.toml`: add `microbiolink-download-fasta = "microbiolink.cli:download_fasta"` under
  `[project.scripts]`; add `"mygene>=3.2.2"` to `dependencies`.

## Migration checklist

### 1. Retrofit modules 1 & 2 into `workflow/` + `utils/`

- [ ] Create `microbiolink/workflow/__init__.py` and `microbiolink/utils/__init__.py`.
- [ ] Move `zscore_filter.py` → `workflow/zscore_filter.py`; remove `read_count_matrix` from it.
- [ ] Move `membrane_filter.py` → `workflow/membrane_filter.py`; update its
      `from . import uniprot_client` to `from ..utils import uniprot_client`.
- [ ] Move `uniprot_client.py` → `utils/uniprot_client.py`.
- [ ] Move `gene_matrix.py` → `utils/gene_matrix.py`; add `read_count_matrix` to it (ported
      unchanged from `zscore_filter.py`, per decision 9).
- [ ] Update `cli.py`'s lazy imports to the new subpackage paths (see CLI section above).
- [ ] Re-run Modules 1 & 2's existing manual sign-off checks to confirm the move is
      behavior-preserving (pure relocation, no logic changes).

### 2. Core modules (Module 3)
- [ ] Add `fetch_fasta_sequences`, `_fetch_fasta_batch`, `build_uniprot_fasta_url`,
      `FASTA_BATCH_SIZE` to `microbiolink/utils/uniprot_client.py`.
- [ ] Add `extract_expressed_gene_symbols_from_count_matrix` to `microbiolink/utils/gene_matrix.py`.
- [ ] Create `microbiolink/utils/id_resolution.py` (`_translate_gene_symbols_to_uniprot`,
      `resolve_uniprot_ids`).
- [ ] Create `microbiolink/workflow/fasta_download.py`.

### 3. CLI wiring
- [ ] Rename `_add_membrane_filter_source_arguments` → `_add_identifier_source_arguments` in
      `cli.py`.
- [ ] Add `_add_human_fasta_arguments`, `_add_microbial_fasta_arguments`,
      `_build_fasta_download_parser`, `_resolve_human_fasta_identifiers`,
      `_resolve_microbial_fasta_identifiers`, `download_fasta` to `cli.py`.
- [ ] Add `mygene>=3.2.2` to `pyproject.toml` dependencies.
- [ ] Add `microbiolink-download-fasta = "microbiolink.cli:download_fasta"` entry point.

### 4. Regression check (per Q11)

No `case_study_output/` fixture exists for FASTA at all (confirmed — only pre-existing input fasta
files, no pipeline-generated output). Both halves need a live baseline run, same situation Modules 1
and 2 hit:

- **Human**: pull `translate_symbol_to_uniprot()` and `fetch_protein_sequences()` from
  `case-study:workflow/get_human_fasta.py`. Run against
  `case_study_input/input/human_protein/Enterocyte_Manual/enterocyte_colon_CD_expressed_genes.csv`
  (plain gene-symbol list). Run the new `fasta_download.download_fasta(..., human_identifiers=...,
  human_id_type='genesymbol')` on the same symbols.
  - Note: this file has **no header row** — `uniprot_client.read_ids` always skips the first line
    as a header, so a direct CLI run against this fixture drops the first symbol. Read the file with
    plain Python for the comparison run instead of going through `read_ids`, so this pre-existing
    `read_ids` behavior (already inherited by Module 2) isn't conflated with a Module-3-specific bug.
  - [ ] Diff: resolved UniProt ID sets match between legacy and new.
  - [ ] Spot-check a handful of FASTA records for content match; note that MyGene.info's underlying
        database may have changed since any original legacy run, so an exact byte-for-byte diff may
        not be achievable — this is expected, not a bug.
- **Microbial**: no case-study fasta baseline exists (case-study has zero bacterial fasta
  capability). Pull `fetch_fasta_sequences()` from
  `origin/MicrobioLink-2.1-beta:microbiolink/get_protein_fasta.py` as the baseline instead. Run
  against `case_study_input/input/bacterial_protein/OMV_proteins.csv` (2065 UniProt IDs, also no
  header row — same caveat as above). Run the new `fasta_download.download_fasta(...,
  microbial_identifiers=..., microbial_id_type='uniprot')` on the same IDs.
  - [ ] Diff: identical FASTA record count and identical accession set between legacy and new
        output (this one *is* a like-for-like comparison, unlike the human/MyGene case).
- **Proteome path**: no fixture proteome ID is present in `case_study_input/`. Pick any small
  real bacterial proteome ID (e.g. from UniProt) for a manual live check of decision 4's two-step
  resolution (`fetch_proteome_table` → `fetch_fasta_sequences`), comparing the resolved ID count
  against the proteome's known protein count on UniProt's website.

### 5. Delete old code (per Q11)

- [ ] Delete `workflow/get_human_fasta.py` entirely — Module 2 already stopped using its
      `get_proteins()`; Module 3 supersedes the remainder (`translate_symbol_to_uniprot`,
      `fetch_protein_sequences`, `main`). This was explicitly flagged as deferred in the Module 2
      plan pending Module 3.
- [ ] `workflow/download_bacterial_proteins.py` stays — it does Pfam/gene-name TSV downloads
      (Module 4 territory), which Module 3 does not touch.
- [ ] Confirm no new code imports from `workflow/get_human_fasta.py` after this module lands.

### 6. Manual sign-off

- [ ] Run the new CLI once with only `--human_*` args, once with only `--microbial_*` args, and once
      with both together, inspecting `human_proteins.fasta` / `microbial_proteins.fasta` by eye each
      time.
- [ ] Confirm the `mygene` translation fix (decision 6) — spot-check a few gene symbols known to
      have straightforward 1:1 UniProt mappings and confirm the returned accessions are correct,
      since this is new/corrected behavior rather than a straight port.
