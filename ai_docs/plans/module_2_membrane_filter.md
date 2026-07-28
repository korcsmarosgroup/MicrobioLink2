# Module 2 — Membrane Protein Filter — Implementation Plan

Detailed plan for the second module described in @ai_docs/plans/microbiolink_refactoring.md, following
the decisions in @ai_docs/plans/refactoring_questions.md. Scope is exactly Module 2: filtering a set of
human or bacterial protein identifiers down to membrane/secreted proteins, and annotating the
survivors with a location.

## Source of truth

Two independent halves, from two different branches, neither of which matches the top-level
plan's output wording literally (see "Output shape" below):

- **Human side** — `case-study:workflow/get_human_fasta.py`, function `get_proteins()`
  (confirmed by Q3). This function is heavily entangled: in one pass it does a non-zero
  expression filter, the OmniPath Intercell membrane filter, gene-symbol→UniProt translation for
  an *unfiltered* code path, and fasta fetching. Only the membrane-filter path is Module 2's
  concern:
  ```python
  pmtm = op.requests.Intercell.get(
      parent=location_filter_list,
      scope=['generic', 'specific'],
      source=['resource_specific', 'composite'],
      entity_type='protein',
  )
  # ... for each (gene, expression) row in the input file:
  if id_type == 'genesymbol':
      if gene in list(pmtm['genesymbol']):
          uniprot = pmtm.loc[pmtm['genesymbol'] == gene, 'uniprot'].values
          if len(uniprot) > 0:
              proteins.append(uniprot[0])
  elif id_type == 'uniprot':
      if gene in list(pmtm['uniprot']):
          proteins.append(gene)
  proteins = list(set(proteins))
  return proteins
  ```
  Today this returns a bare, deduplicated `list[str]` of UniProt IDs — **no annotation is kept**,
  even for the `genesymbol` path it discards the gene symbol and keeps only the resolved UniProt
  ID. There is no proteome-ID branch at all.

- **Bacterial side** — `MicrobioLink-2.1-beta:microbiolink_api/microbiome.py`, function
  `filter_bacterial_domain_table_by_location()` (found after Q3 initially missed it; corrected by
  the "Correction" note in the questions doc). This one already operates on a DataFrame:
  ```python
  normalized_filters = [location.lower() for location in location_filters]
  location_series = domain_table[LOCATION_COLUMN].fillna('').astype(str).str.lower()
  mask = location_series.apply(
      lambda value: any(location in value for location in normalized_filters),
  )
  return domain_table.loc[mask].reset_index(drop=True)
  ```
  where `LOCATION_COLUMN = 'Subcellular location [CC]'`, populated by fetching UniProt with
  `fields=['accession', 'xref_pfam', 'gene_names', 'cc_subcellular_location']`. This is a
  substring match against UniProt's free-text subcellular-location annotation, not a fixed
  vocabulary. It depends on `fetch_bacterial_domain_table_from_ids()` in the same file, which in
  turn calls `download_protein_list_with_fields()` / `download_proteome_with_fields()` from
  `microbiolink/download_protein_domains.py` (beta).

- **Case-study's own bacterial fetch code**, `case-study:workflow/download_bacterial_proteins.py`
  (unmodified on this branch — confirmed no diff vs. `case-study`), only exposes
  `download_protein_list()` / `download_proteome()` with **hardcoded**
  `fields=accession,xref_pfam,gene_names` — no way to also request the location column. Per Q9's
  precedence rule (flat `microbiolink` wins only where functionality is *the same*), this is a
  case where it isn't the same: case-study's version cannot express "also fetch location," so the
  parameterized `*_with_fields()` versions from beta are the ones ported in, same exception class
  as the bacterial location filter itself already being sourced from beta.

## Decisions from grilling

1. **Output shape**: a `pandas.DataFrame`, one row per surviving protein, with columns
   `uniprot_id` and `location_annotation` — always keyed by UniProt ID regardless of what
   `id_type` was supplied. This matches what the legacy human code already does implicitly (it
   resolves `genesymbol` matches back to `uniprot` before returning) and what the bacterial table
   already does (`Entry` is always the UniProt accession, whether fetched via `id_type='Uniprot'`
   or `id_type='UP'`). Neither existing implementation actually produces this shape today —
   human returns a bare ID list with no annotation at all; bacterial returns the full domain
   table (`Entry`, `Pfam`, `Gene Names`, raw location text) rather than a two-column standardized
   result. Both are narrowed/reshaped to the standard two-column form as part of this port.
2. **Human + proteome ID: out of scope.** Confirmed — a full human proteome is never a realistic
   input for this filter (unlike bacterial, where "give me the whole proteome" is a normal case).
   `filter_membrane_proteins(..., species='human')` only accepts `id_type` of `'uniprot'` or
   `'genesymbol'`. `id_type='proteome'` is only valid for `species='microbial'`.
3. **No gene-symbol→UniProt translation inside Module 2.** The OmniPath Intercell table already
   carries both `uniprot` and `genesymbol` columns per record, so a `genesymbol` identifier is
   matched directly against `pmtm['genesymbol']` and the UniProt ID is read off the same matching
   row — exactly as legacy `get_proteins()` does. Q5's shared `translate_symbol_to_uniprot()`
   utility is not needed here; it's only exercised by legacy `get_proteins()`'s *unfiltered*
   branch (no `location_filter_list` at all), which does no membrane filtering and is really
   Module 3's ID-resolution concern, not Module 2's. Left for the Module 3 plan.
4. **Non-zero expression filter stays out.** Already assigned to Module 3 by Q4. Module 2 takes
   already-resolved identifiers (or a count matrix it extracts symbols from) — it has no opinion
   on expression values.
5. **Gene-symbol-from-count-matrix extraction is a separate, shared helper**, not fused into the
   filter function. The top-level plan repeats the identical input wording for Module 2 and
   Module 3 ("gene symbols from a gene count matrix ... we will then need a helper function"), so
   per Q5's "shared utility, not per-module duplication" precedent this lives in its own module
   and Module 3 imports the same function later rather than reimplementing it.
6. **Shared low-level UniProt fetch utility, sourced from beta, not case-study** (see "Source of
   truth" above for why this is a legitimate Q9 exception). This utility is written to be reused
   by Modules 3 & 4 as well, per Q6 — the field list is the only thing that varies per caller.
7. **Annotation content legitimately differs by species** — human gets deduplicated OmniPath
   `parent` location categories (e.g. `secreted`), bacterial gets raw free-text UniProt
   subcellular-location annotation (e.g. `Cell outer membrane`). This is not a bug to reconcile;
   it's what each upstream source actually provides. Both are exposed under the same
   `location_annotation` column name so downstream code has one shape to consume regardless of
   species.
   - **Note for implementation time**: the exact column name OmniPath's `Intercell.get()` uses
     for the matched location category (`parent`, `category`, or something else) needs confirming
     against the installed `omnipath` client version — it's not pinned down by legacy behavior,
     since legacy code never reads that column at all (only membership-checks `uniprot`/
     `genesymbol`). Treat this the same way Module 1 flagged its KDE grid resolution: an unstated
     legacy parameter to nail down, not guess, before merging.
8. **No new validation on `location_filters`.** Legacy code validates neither the human category
   strings (an invalid `parent` value just yields zero Intercell rows) nor the bacterial substring
   terms (any string is a valid substring to search for — the top-level plan's "outer membrane or
   plasma membrane" is example usage, not an enum). Keep it that way; consistent with Module 1's
   "no new edge-case guards" decision.
9. **File layout**: extends Q7's "one file per capability, single function dispatching on
   `id_type`/`species`" pattern to Module 2, even though Q7 was scoped to Modules 3 & 4 — the
   input contract is worded identically in the top-level plan, so the same shape applies here.
   - `microbiolink/uniprot_client.py` — shared low-level fetch utility (decision 6), reused later
     by Modules 3 & 4.
   - `microbiolink/gene_matrix.py` — shared gene-symbol-extraction helper (decision 5), reused
     later by Module 3.
   - `microbiolink/membrane_filter.py` — Module 2's own dispatch + per-species logic.

## Target implementation

### `microbiolink/uniprot_client.py` (core, no argparse)

```python
UNIPROT_STREAM_BASE_URL = 'https://rest.uniprot.org/uniprotkb/stream?'
DEFAULT_UNIPROT_FIELDS = ['accession', 'xref_pfam', 'gene_names']
UNIPROT_BATCH_SIZE = 1000


def read_ids(filename: PathLike, separator: str, id_column: int) -> list[str]:
    """Read UniProt or proteome identifiers from a delimited file."""


def fetch_protein_table(
    identifiers: list[str],
    fields: list[str] | None = None,
) -> pd.DataFrame:
    """Fetch a UniProt annotation table for a batch of protein accessions."""
```

- Ports `build_uniprot_accession_query()` / `build_uniprot_stream_url()` /
  `download_protein_list_with_fields()` / `download_proteome_with_fields()` from
  beta's `microbiolink/download_protein_domains.py`, plus the batch-and-retry-by-splitting
  wrapper (`_download_protein_identifier_frames` / `_download_protein_identifier_batch`) from
  beta's `microbiolink_api/microbiome.py`, since that retry logic is a real robustness
  improvement over both legacy scripts (a single malformed ID no longer fails an entire 1000-ID
  batch).
- Unlike both legacy sources, this returns a parsed `pd.DataFrame` directly (via
  `pd.read_csv(StringIO(...), sep='\t')`) instead of raw TSV text — every existing caller
  immediately parses the text anyway, so this removes a repeated parsing step at each call site.
- `fetch_proteome_table(proteome_id: str, fields: list[str] | None = None) -> pd.DataFrame` — the
  proteome-ID equivalent, adding a `Proteome_ID` column (matches legacy behavior in both
  case-study's and beta's proteome download loop).
- No file I/O beyond the network fetch itself, no argparse, no `print`.

### `microbiolink/gene_matrix.py` (core, no argparse)

```python
def extract_gene_symbols_from_count_matrix(count_matrix: pd.DataFrame) -> list[str]:
    """Extract the gene symbol row index from a gene count matrix as a plain list."""
```

- Trivial (`count_matrix.index.tolist()`), but named and isolated per decision 5 so Module 3 can
  import the same function rather than re-deriving it.

### `microbiolink/membrane_filter.py` (core, no argparse)

```python
LOCATION_FIELD = 'cc_subcellular_location'
LOCATION_COLUMN = 'Subcellular location [CC]'


def filter_human_membrane_proteins(
    identifiers: list[str],
    id_type: str,
    location_filters: list[str],
) -> pd.DataFrame:
    """Filter human proteins to membrane/secreted ones via OmniPath Intercell."""


def filter_bacterial_membrane_proteins(
    identifiers: list[str],
    id_type: str,
    location_filters: list[str],
) -> pd.DataFrame:
    """Filter bacterial proteins to membrane/secreted ones via UniProt location text."""


def filter_membrane_proteins(
    identifiers: list[str],
    id_type: str,
    species: str,
    location_filters: list[str],
) -> pd.DataFrame:
    """Filter protein identifiers to membrane/secreted proteins, dispatching by species."""
```

- `filter_human_membrane_proteins`: `id_type` is `'uniprot'` or `'genesymbol'` (decision 2).
  Lazily imports `omnipath` inside the function body (matches the confirmed lazy-import
  convention for heavy dependencies). Calls `Intercell.get(parent=location_filters, ...)` exactly
  as legacy does, filters to rows whose `uniprot`/`genesymbol` column is in `identifiers`, groups
  by the resolved `uniprot` value and joins the matched location categories (decision 7) into
  `location_annotation`, returns the two-column DataFrame.
- `filter_bacterial_membrane_proteins`: `id_type` is `'uniprot'` or `'proteome'`. Calls
  `uniprot_client.fetch_protein_table()` / `fetch_proteome_table()` with
  `fields=['accession', LOCATION_FIELD]`, then applies the same substring-match logic as beta's
  `filter_bacterial_domain_table_by_location()` (decision 1/7), renaming `Entry` →
  `uniprot_id` and `Subcellular location [CC]` → `location_annotation` in the returned frame.
- `filter_membrane_proteins`: the public dispatcher — validates `species` is `'human'` or
  `'microbial'` and routes to the matching function above. This is the function Module 2's CLI
  and any other module call into; the two `filter_*_membrane_proteins()` functions stay
  independently importable for direct/library use.

### `microbiolink/cli.py` (argparse only, per Q2)

Add:

```python
def membrane_filter() -> int:
    from . import gene_matrix, membrane_filter as membrane_filter_module, uniprot_client, zscore_filter

    parser = argparse.ArgumentParser(
        description='Filter human or microbial proteins down to membrane/secreted proteins.',
    )
    parser.add_argument('-i', '--input_file', required=True, help='Identifier list, or a count matrix if --from-count-matrix is set.')
    parser.add_argument('-id', '--id_type', required=True, choices=['uniprot', 'genesymbol', 'proteome'])
    parser.add_argument('-sp', '--species', required=True, choices=['human', 'microbial'])
    parser.add_argument('-lfl', '--location_filters', required=True, nargs='+')
    parser.add_argument('--from-count-matrix', action='store_true', help='Treat --input_file as a gene count matrix and extract gene symbols from it.')
    parser.add_argument('-sep', '--sep', default=',', help='Field separator for a plain identifier list (ignored with --from-count-matrix).')
    parser.add_argument('-col', '--id_column', type=int, default=1, help='One-based identifier column for a plain list (ignored with --from-count-matrix).')
    parser.add_argument('-o', '--output_file', required=True)
    args = parser.parse_args()

    if args.from_count_matrix:
        count_matrix = zscore_filter.read_count_matrix(args.input_file)
        identifiers = gene_matrix.extract_gene_symbols_from_count_matrix(count_matrix)
    else:
        identifiers = uniprot_client.read_ids(args.input_file, args.sep, args.id_column)

    result = membrane_filter_module.filter_membrane_proteins(
        identifiers,
        id_type=args.id_type,
        species=args.species,
        location_filters=args.location_filters,
    )
    result.to_csv(args.output_file, index=False)
    return 0
```

- Reuses Module 1's `read_count_matrix()` for the `--from-count-matrix` path rather than
  duplicating a CSV-reading function.
- `pyproject.toml` entry point: `microbiolink-membrane-filter = "microbiolink.cli:membrane_filter"`.

## Inputs/Outputs recap (from top-level plan, narrowed by decisions above)

- **Input**: a list of UniProt IDs, a UniProt proteome ID (microbial only, per decision 2), human
  gene symbols, or gene symbols extracted from a gene count matrix; plus `species` and
  `location_filters` chosen by the user.
- **Output**: a `pandas.DataFrame` with columns `uniprot_id` and `location_annotation`, one row
  per protein that matched at least one of the requested location filters (decision 1).

## Migration checklist

### 1. Core modules

- [ ] Create `microbiolink/uniprot_client.py`: port `read_ids()`, the accession-query/stream-URL
      builders, `fetch_protein_table()` (batched, with retry-by-splitting on HTTP 400), and
      `fetch_proteome_table()` — sourced from beta per decision 6, not case-study's hardcoded-field
      versions.
- [ ] Create `microbiolink/gene_matrix.py`: `extract_gene_symbols_from_count_matrix()`.
- [ ] Create `microbiolink/membrane_filter.py`: `filter_human_membrane_proteins()`,
      `filter_bacterial_membrane_proteins()`, `filter_membrane_proteins()`.
- [ ] Confirm the OmniPath column name for the matched location category against the installed
      `omnipath` package version (decision 7's flagged unknown) before relying on it.
- [ ] Confirm no new validation is added on `location_filters` beyond what legacy does (decision
      8).

### 2. CLI wiring

- [ ] Add `membrane_filter()` to `microbiolink/cli.py` (argparse only, per Q2).
- [ ] Update `pyproject.toml` entry point:
      `microbiolink-membrane-filter = "microbiolink.cli:membrane_filter"`.

### 3. Regression check (per Q11)

No case-study fixture exists for either half of Module 2 (`case_study_output/` has nothing
membrane/location-related) — same situation Module 1 hit, so both checks require a live baseline
run before any old code is touched:

- **Human**: pull `get_proteins()` from `case-study:workflow/get_human_fasta.py` (already
  unmodified on this branch, but pull explicitly per Module 1's precedent for an unambiguous
  baseline). Run it against
  `case_study_input/input/human_transcriptomics/colon_BEST4_enterocyte_CD.csv` with a chosen
  `location_filter_list` (e.g. `secreted`). Run the new `filter_human_membrane_proteins()` on the
  same file's gene symbols (id_type=`genesymbol`) with the same filter.
  - [ ] Diff: the new function's `uniprot_id` column, as a set, equals the legacy function's
        returned list, as a set.
  - [ ] Note: legacy never returned an annotation, so there is nothing to diff the new
        `location_annotation` column against — that column is new territory, not a ported value.
        Cross-check a handful of rows by hand against the OmniPath Intercell response instead.
- **Bacterial**: pull `filter_bacterial_domain_table_by_location()` and
  `fetch_bacterial_domain_table_from_ids()` from
  `origin/MicrobioLink-2.1-beta:microbiolink_api/microbiome.py`. Run them against
  `case_study_input/input/bacterial_protein/OMV_proteins.csv` (2065 UniProt IDs, `id_type='Uniprot'`)
  with a chosen `location_filters` (e.g. `['outer membrane']`). Run the new
  `filter_bacterial_membrane_proteins()` on the same IDs with the same filter.
  - [ ] Diff: the new function's `uniprot_id` set matches the legacy filtered table's `Entry` set.
  - [ ] Diff: the new function's `location_annotation` values match the legacy table's raw
        `Subcellular location [CC]` text for the same accessions (this one *is* a like-for-like
        value, unlike the human side).

### 4. Delete old code (per Q11) — partial only, see caveat

Unlike Module 1, no old file is fully superseded by Module 2 alone:

- `workflow/get_human_fasta.py` also contains `translate_symbol_to_uniprot()`,
  `fetch_protein_sequences()`, and `main()` — all Module 3 territory (decision 3), not ported
  here. **Do not delete this file yet.** Only stop calling/importing its `get_proteins()` for
  membrane filtering from any new code; the file itself is deleted once the Module 3 plan lands
  and supersedes the rest of it.
- `workflow/download_bacterial_proteins.py` is not the source for anything in this module
  (decision 6/9 — its hardcoded-fields functions aren't equivalent to what's needed). It stays
  until Modules 3/4 adopt `microbiolink/uniprot_client.py` in its place.
- [ ] Confirm no new code path imports from either of the above two files after this module lands.

### 5. Manual sign-off

- [ ] Run the new CLI once end-to-end for the human path (`--species human --id_type genesymbol
      --from-count-matrix`) and once for the bacterial path (`--species microbial --id_type
      uniprot`), and inspect both output files by eye.
- [ ] Confirm the OmniPath annotation column resolved in step 1 actually contains the expected
      category labels (e.g. `secreted`) — this was flagged as unverified against legacy behavior
      and needs a human look before merging.
