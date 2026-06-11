# Context

MicrobioLink currently predicts host-microbe interactions by finding ELM short linear motifs in **human** proteins and matching them to Pfam domains in **bacterial** proteins. The user wants to reverse this: find ELM motifs in **bacterial** proteins and match them to Pfam domains in **human** proteins. This answers the complementary biological question — which human protein domains could be engaged by bacterial mimicry motifs?

---

# Reversed Pipeline Overview

```
z_score_filter_terminal.py         download_bacterial_proteins.py
         ↓                                      ↓
download_human_domains.py          get_bacterial_fasta.py
(human accessions + Pfam)          (bacterial FASTA sequences)
         ↓                                      ↓
              reverse_DMI.py
   (scan bacterial FASTA for ELM motifs,
    match motifs → Pfam domains → human proteins)
         ↓
  Output: bacterial_protein;motif;start;end;human_domain;human_protein
```

Both new download steps are powered by two new **generic** modules that centralise the UniProt download logic currently split across `download_bacterial_proteins.py` and `get_human_fasta.py`.

---

# Code Style (`.agents/code_style.md`)

All new and modified files must follow these rules:

- **PEP8 + Google style guide** throughout.
- **Docstrings**: Napoleon (Google) style; opening and closing `"""` on their own lines; use type hints in signatures, not in docstring bodies.
- **Paths**: Use `pathlib.Path` for **all** filesystem code — no `os.path`, no raw strings for paths. Any existing `os.path` calls in files being modified must be updated.
- **Package management**: Run tests with `uv run pytest`, not `pytest` directly. Verify any new dependencies are added to `pyproject.toml`.
- **No throwaway scripts**: all verification is done through the `tests/` test suite.

---

# New Files

## 1. `microbiolink/download_protein_domains.py` *(new, generic)*

**Purpose:** Download Pfam domain annotations for any list of UniProt accessions, regardless of organism.

**Input:** list of UniProt accessions (passed programmatically or via a file)
**Output:** TSV with `accession`, `xref_pfam`, `gene_names`

**Key logic — extracted and generalised from `download_bacterial_proteins.py`:**

- `build_uniprot_accession_query(uniprot_ids)` — build percent-encoded OR query
- `build_uniprot_stream_url(query, fields)` — assemble full stream URL
- `download_protein_list_with_fields(uniprot_ids, fields)` — fetch TSV for an accession list
- `download_proteome_with_fields(proteome_id, fields)` — fetch TSV for a full proteome (UP identifier)
- `UNIPROT_BATCH_SIZE = 1000`, `DEFAULT_UNIPROT_FIELDS`

**CLI args:**

```
--id_list     Path to file containing accessions or proteome IDs
--sep         Field separator
--id_type     uniprot | UP
--id_column   One-based column number
--output      Output TSV path
```

**Effect on existing code:**
`download_bacterial_proteins.py` is updated to import its core functions from `download_protein_domains.py` rather than defining them inline. Its CLI and `main` are otherwise unchanged. To preserve the public API, `download_bacterial_proteins.py` re-exports the moved names (`build_uniprot_accession_query`, `build_uniprot_stream_url`, `download_protein_list_with_fields`, `download_proteome_with_fields`, `UNIPROT_BATCH_SIZE`, `DEFAULT_UNIPROT_FIELDS`) so any existing code that imports them from `download_bacterial_proteins` continues to work without modification.

---

## 2. `microbiolink/get_protein_fasta.py` *(new, generic)*

**Purpose:** Download FASTA sequences for any list of UniProt accessions, regardless of organism.

**Key logic — extracted and generalised from `fetch_protein_sequences` in `get_human_fasta.py`:**

- `fetch_fasta_sequences(uniprot_ids)` — calls UniProt stream endpoint with `format=fasta`; uses proper `response.raise_for_status()` instead of manual status check; removes the debug `print(url)` present in the original
- Batching loop driven by caller (consistent with the 100-accession batch pattern already in `get_human_fasta.py`)

**CLI args:**

```
--id_list     Path to file containing UniProt accessions
--sep         Field separator
--id_column   One-based column number
--output      Output FASTA path
--batch_size  Accessions per request (default 100)
```

**Effect on existing code:**
`get_human_fasta.py` is updated so `fetch_protein_sequences` delegates to `fetch_fasta_sequences` from `get_protein_fasta.py`. The error handling difference must be handled explicitly: the generic `fetch_fasta_sequences` uses `response.raise_for_status()` (raising on failure), but `get_human_fasta.py` currently silently skips failed batches. To avoid a breaking change in the forward pipeline, `fetch_protein_sequences` in `get_human_fasta.py` wraps the call in a try/except that catches `requests.HTTPError` and prints a warning, preserving the existing skip-and-continue behaviour.

---

## 3. `microbiolink/download_human_domains.py` *(new, human-specific)*

**Purpose:** Download human protein Pfam domain annotations from expressed genes (z-score filtered output).

**Input:** z-score filtered CSV (output of `z_score_filter_terminal.py`)
**Output:** TSV with `accession`, `xref_pfam`, `gene_names` — identical format to the bacterial domain file

**Key logic:**

- Read the z-score filtered CSV; keep rows where the expression value in the nominated column is non-NaN and non-zero
- If `--id_type genesymbol`: call `translate_symbol_to_uniprot` (imported from `get_human_fasta.py`) to convert to Swiss-Prot UniProt accessions
- If `--id_type uniprot`: use accessions directly
- Call `download_protein_list_with_fields` from `download_protein_domains.py` in batches of 1000

**CLI args:**

```
--gene_expression   Path to z-score filtered CSV
--id_type           genesymbol | uniprot
--sep               Separator used in the CSV
--output            Output TSV file path
```

---

## 4. `microbiolink/get_bacterial_fasta.py` *(new, bacterial-specific)*

**Purpose:** Download FASTA sequences for **all** bacterial proteins in a proteome or accession list — including proteins without Pfam domain annotations, which are excluded from the domain TSV but may still carry ELM motifs.

**Input:** Same upstream source as `download_bacterial_proteins.py` — either a proteome ID list or a UniProt accession list (the original `--id_list` file, before any domain filtering)
**Output:** FASTA file of bacterial protein sequences

**Key logic:**

- If `--id_type UP` (proteome): query UniProt FASTA endpoint directly using the proteome query format `(proteome:{proteome_id})`, one proteome at a time
- If `--id_type uniprot` (accession list): read accessions from the ID list file and call `fetch_fasta_sequences` from `get_protein_fasta.py` in batches of 100

**CLI args:**

```
--id_list     Path to the same ID list used for download_bacterial_proteins.py
--sep         Field separator
--id_type     uniprot | UP
--id_column   One-based column number
--output      Output FASTA file path
--batch_size  Accessions per request for uniprot mode (default 100)
```

---

## 5. `microbiolink/reverse_DMI.py`

**Purpose:** Reversed domain-motif interaction prediction — bacterial motifs × human domains.

**Inputs:**

- Bacterial FASTA (from `get_bacterial_fasta.py`)
- ELM regex file (same resource file as forward DMI)
- Motif-domain interaction file (same resource file as forward DMI)
- Human protein domain TSV (from `download_human_domains.py`)

**Output:** Semicolon-separated table:

```
# Bacterial Protein;Motif;Start;End;Human Domain;Human Protein
```

**Key logic (reversed from `DMI.py`):**

1. `read_fasta_sequences(bacterial_fasta)` → bacterial sequences *(reuse from `DMI.py`)*
2. `parse_elm_regex(elm_regex_file)` *(reuse from `DMI.py`)*
3. Filter ELM motifs: remove all entries whose identifier starts with `CLV_` (cleavage site motifs that mediate proteolytic digestion rather than domain binding, e.g. `CLV_C14_Caspase3-7`, `CLV_PCSK_FUR_1`). This is implemented as `filter_cleavage_motifs(elm_regex)` which drops any key beginning with `"CLV_"`.
4. `parse_motif_domain(motif_domain_file)` *(reuse from `DMI.py`)*
5. `parse_protein_domain(human_domain_file)` → `{pfam_id: [human_uniprot_ids]}` *(reuse from `DMI.py`)*
6. `create_uniprot_motif_dict(bacterial_sequences, elm_regex)` → bacterial proteins with motif hits *(reuse from `DMI.py` unchanged — CLV motifs already removed from elm_regex before this call)*
7. Main matching loop (reversed):
   - For each motif → look up its interacting Pfam domains
   - For each domain → look up human proteins carrying that domain
   - For each bacterial protein with that motif → emit one row

**CLI args:**

```
-fasta        Bacterial FASTA file
-motif        ELM regex file
-interaction  Motif-domain interaction file
-domain       Human protein domain file
-o            Output file
```

**Functions to reuse from `DMI.py`:**
`extract_uniprot_id`, `read_fasta_sequences`, `parse_elm_regex`, `parse_motif_domain`, `parse_protein_domain`, `create_uniprot_motif_dict`

`filter_cleavage_motifs` is defined locally in `reverse_DMI.py` only — `DMI.py` is not modified.

---

## 6. Update `microbiolink/cli.py`

Add five new entry points following the existing dynamic-import pattern:

- `download_protein_domains()` → `.download_protein_domains`
- `get_protein_fasta()` → `.get_protein_fasta`
- `download_human_domains()` → `.download_human_domains`
- `get_bacterial_fasta()` → `.get_bacterial_fasta`
- `reverse_dmi()` → `.reverse_DMI`

No changes needed to existing entry points — `download_bacterial_proteins` and `get_human_fasta` keep their current CLI signatures.

---

# Tests

All new pipeline components require a corresponding test file in `tests/`. Tests are pytest-style functions (`test_` prefix), use fixtures, and assert directly. No ad hoc scripts. Run with `uv run pytest`.

## `tests/test_download_protein_domains.py`

- `test_build_uniprot_accession_query_single` — single accession produces correctly encoded clause
- `test_build_uniprot_accession_query_multiple` — multiple accessions joined with `+OR+` and wrapped in outer parens
- `test_build_uniprot_stream_url_default_fields` — URL contains all three default fields encoded
- `test_build_uniprot_stream_url_custom_fields` — custom fields override defaults
- `test_download_protein_list_with_fields_calls_correct_url` — mock `requests.get`; assert correct URL formed and `raise_for_status` called
- `test_download_proteome_with_fields_calls_correct_url` — mock `requests.get`; assert proteome query format used
- `test_main_writes_tsv_output` — write a tmp ID file, mock HTTP, assert output TSV written with correct header

## `tests/test_get_protein_fasta.py`

- `test_fetch_fasta_sequences_returns_fasta` — mock `requests.get` returning FASTA text; assert result matches
- `test_fetch_fasta_sequences_raises_on_http_error` — mock 500 response; assert `requests.HTTPError` raised
- `test_main_writes_fasta_output` — write a tmp ID file, mock HTTP, assert FASTA written

## `tests/test_download_human_domains.py`

- `test_read_expressed_genes_uniprot_filters_nan` — rows with `NaN` expression excluded; valid UniProt IDs returned
- `test_read_expressed_genes_uniprot_filters_zero` — rows with zero expression excluded
- `test_read_expressed_genes_genesymbol` — mock `translate_symbol_to_uniprot`; assert Swiss-Prot IDs extracted, TrEMBL entries dropped
- `test_main_uniprot_mode` — end-to-end with tmp input CSV and mocked `download_protein_list_with_fields`; assert output TSV correct

## `tests/test_get_bacterial_fasta.py`

- `test_proteome_mode_builds_correct_url` — mock `requests.get`; assert proteome query format in URL
- `test_uniprot_mode_batches_correctly` — provide 150 accessions; assert two batched calls to `fetch_fasta_sequences`
- `test_main_writes_fasta` — end-to-end with tmp ID file and mocked HTTP; assert FASTA written

## `tests/test_reverse_DMI.py`

- `test_filter_cleavage_motifs_removes_clv` — dict with `CLV_*` and non-`CLV_*` keys; assert only `CLV_*` removed
- `test_filter_cleavage_motifs_empty_input` — empty dict returns empty dict
- `test_reverse_dmi_main_output_format` — run `main` with tmp fixture files (minimal FASTA, ELM regex with one non-CLV motif, motif-domain mapping, human domain TSV); assert output rows are `bacterial;motif;start;end;domain;human` format
- `test_reverse_dmi_excludes_clv_motifs` — ELM regex fixture includes a `CLV_` motif; assert no output rows contain it

## `tests/test_import_compatibility.py` *(import verification)*

- `test_import_build_uniprot_accession_query_from_download_bacterial_proteins` — `from microbiolink.download_bacterial_proteins import build_uniprot_accession_query` does not raise
- `test_import_build_uniprot_stream_url_from_download_bacterial_proteins` — same for `build_uniprot_stream_url`
- `test_import_download_protein_list_with_fields_from_download_bacterial_proteins` — same
- `test_import_download_proteome_with_fields_from_download_bacterial_proteins` — same
- `test_import_constants_from_download_bacterial_proteins` — `UNIPROT_BATCH_SIZE` and `DEFAULT_UNIPROT_FIELDS` importable from `download_bacterial_proteins`
- `test_import_fetch_protein_sequences_from_get_human_fasta` — `from microbiolink.get_human_fasta import fetch_protein_sequences` does not raise

---

# Verification

1. **Import compatibility**: `uv run pytest tests/test_import_compatibility.py` — all re-exported names from refactored modules import without error.
2. **Forward pipeline regression**: confirm `get_human_fasta` skips (not crashes) on a failed batch.
3. **New generic modules**: `uv run pytest tests/test_download_protein_domains.py tests/test_get_protein_fasta.py`
4. **New specific modules**: `uv run pytest tests/test_download_human_domains.py tests/test_get_bacterial_fasta.py tests/test_reverse_DMI.py`
5. **Full suite**: `uv run pytest` — no regressions across all modules.
