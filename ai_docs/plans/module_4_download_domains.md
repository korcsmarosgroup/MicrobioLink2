# Module 4 — Downloading the Domains — Implementation Plan

Detailed plan for the fourth module described in @ai_docs/plans/microbiolink_refactoring.md, following the
decisions in @ai_docs/plans/refactoring_questions.md and the precedent set by
@ai_docs/plans/module_2_membrane_filter.md and @ai_docs/plans/module_3_fasta_download.md. Scope is exactly
Module 4: resolving a set of human or microbial protein identifiers (same input contract as
Module 3 — uniprot IDs, a proteome ID, gene symbols, or gene symbols from a count matrix) and
fetching their Pfam domain annotations from UniProt.

## Context

Modules 1–3 are implemented in `microbiolink/workflow/` (`zscore_filter.py`, `membrane_filter.py`,
`fasta_download.py`) and `microbiolink/utils/` (`uniprot_client.py`, `id_resolution.py`), wired
into `microbiolink/cli.py`. Reading the actual code (not just the Module 3 plan document, which
turned out to be aspirational in a few places — e.g. it describes a `gene_matrix.py` that was
never built) shows Module 4 needs very little new plumbing:

- `utils/uniprot_client.fetch_protein_table(identifiers, fields=None)` already defaults to
  `DEFAULT_UNIPROT_FIELDS = ['accession', 'xref_pfam', 'gene_names']` — Pfam data (`xref_pfam`)
  comes back for free with no new fields/utils work. Same for `fetch_proteome_table`.
- `utils/id_resolution.resolve_uniprot_ids(identifiers, id_type)` is already generic over
  `id_type` (`uniprot`/`genesymbol`/`proteome`) and was explicitly designed (Module 3 plan,
  decision 7) to be reused unchanged by Module 4, swapping in `fetch_protein_table` for the fasta
  fetch.
- `workflow/fasta_download.py` is the direct structural template: a private per-species helper
  plus a public dispatcher taking parallel `human_*`/`microbial_*` optional argument groups, no
  explicit `species` parameter.

Legacy precedent (case-study's `workflow/download_bacterial_proteins.py`, still present at this
repo's top-level `workflow/` folder; beta's `microbiolink/download_protein_domains.py` /
`download_human_domains.py`; `microbiolink_api/microbiome.py`) all build a dict keyed by **Pfam**
(`pfam_id -> [uniprot_ids]`) for downstream DMI consumption (`workflow/DMI.py`'s
`parse_protein_domain`). **Revised while planning Module 5**: the top-level plan originally read
"Dictionary with key as Uniprot ID and value as pfam domains" (protein-keyed), but Module 5's DDI
lookup is keyed by Pfam domain pair, so a protein-keyed input forced an extra invert step in every
consumer. The top-level plan has been updated to match the legacy direction instead —
`pfam_id -> [uniprot_ids]` — so Module 5 (and any future Pfam-keyed consumer) can use Module 4's
output as-is. **This plan document reflects that revision; the already-implemented
`microbiolink/workflow/domain_download.py` on this branch still returns the old protein-keyed shape
and needs to be updated to match before Module 5 is built.**

No `case_study_output` fixture exists for domain download. `case_study_input/input/bacterial_protein/BT_BEV_domains.tsv`
is a pre-existing fixture in the legacy `Entry\tPfam\tGene Names` TSV shape, usable as a live
regression baseline (re-fetch its `Entry` IDs and diff Pfam sets).

## Decisions

1. **Core module returns a pure Python dict, no file I/O in `workflow/domain_download.py`.** The
   top-level plan describes Module 4's output as "a dictionary", unlike Module 3's literal "stored
   as fasta files ... in a folder" wording — matches Module 2's precedent of pure-return workflow
   functions, with `cli.py` doing all persistence.
2. **Reuse `id_resolution.resolve_uniprot_ids()` unchanged** for both species — no `species`
   parameter on the public function, mirroring `fasta_download.download_fasta`'s shape exactly.
3. **Reuse `uniprot_client.fetch_protein_table()` with no `fields` override** — defaults already
   include `xref_pfam`. No new utils code needed for the fetch step itself.
4. **CLI output format: Pfam-keyed TSV**, `Pfam`/`Entries` columns, one file per species in an
   output folder (`human_domains.tsv` / `microbial_domains.tsv`), written directly from the
   Pfam-keyed dict with a semicolon-per-Entry-plus-trailing-semicolon convention (e.g.
   `Q9ABC1;Q9ABC2;`) — no invert step anywhere, CLI or core. This is a deliberate format change from
   `BT_BEV_domains.tsv`'s legacy `Entry`/`Pfam` (protein-keyed) convention; the regression check
   below compares against an inverted view of that legacy fixture instead of a literal diff.
5. **CLI arg-builder/resolver reuse via rename**: `cli.py`'s `_add_human_fasta_arguments` /
   `_add_microbial_fasta_arguments` and `_resolve_human_fasta_identifiers` /
   `_resolve_microbial_fasta_identifiers` are renamed to generic `_add_human_identifier_arguments`
   / `_add_microbial_identifier_arguments` and `_resolve_human_identifiers` /
   `_resolve_microbial_identifiers` (behavior unchanged) since Module 4 needs the exact same
   identifier-resolution flags/logic as Module 3 — the same rename precedent Module 3 itself
   applied to Module 2's `_add_membrane_filter_source_arguments`.
6. **Old `workflow/download_bacterial_proteins.py`** (top-level, pre-refactor, still present in the
   repo) is fully superseded by this module and is deleted once Module 4 lands and is verified,
   per Q11.

## Target implementation

### `microbiolink/workflow/domain_download.py` (new, core, no argparse)

```python
def domain_table_to_mapping(domain_table: pd.DataFrame) -> dict[str, list[str]]:
    """Convert a UniProt Pfam table into a pfam_id -> uniprot_ids mapping."""


def _download_species_domains(identifiers: list[str], id_type: str) -> dict[str, list[str]]:
    """Resolve identifiers to UniProt accessions and fetch their Pfam domains."""


def download_domains(
    human_identifiers: list[str] | None = None,
    human_id_type: str | None = None,
    microbial_identifiers: list[str] | None = None,
    microbial_id_type: str | None = None,
) -> dict[str, dict[str, list[str]]]:
    """Download Pfam domains for human and/or microbial proteins."""
```

- `download_domains` mirrors `fasta_download.download_fasta`'s shape: same
  `if human_identifiers is None and microbial_identifiers is None: raise ValueError(...)` guard,
  same dict-comprehension-over-`species_requests` dispatch — except it returns
  `{species: {pfam_id: [uniprot_id, ...]}}` instead of writing files and returning paths.
- `_download_species_domains`: `uniprot_ids = id_resolution.resolve_uniprot_ids(identifiers, id_type)`;
  `domain_table = uniprot_client.fetch_protein_table(uniprot_ids)`; returns
  `domain_table_to_mapping(domain_table)`.
- `domain_table_to_mapping`: iterates rows; for each row, if `row['Pfam']` is NaN/empty the protein
  is skipped (it carries no domains, so it can never appear as a value in a Pfam-keyed mapping);
  otherwise for each `pfam in row['Pfam'].split(';')` (non-empty), appends `row['Entry']` to
  `mapping.setdefault(pfam, [])`. This converts UniProt's raw protein-row fetch response into the
  Pfam-keyed shape in a single pass — no separate invert step. **Public** (no leading underscore)
  purely as a matter of module hygiene (it's a meaningful, independently-testable conversion), not
  because Module 5 needs to call it — Module 5's CLI reads Module 4's own Pfam-keyed TSV output
  directly (already the right shape, see @ai_docs/plans/module_5_ddi.md), it does not re-parse UniProt's
  raw response format.

### `microbiolink/cli.py` (extend, argparse only)

- Rename `_add_human_fasta_arguments` → `_add_human_identifier_arguments`,
  `_add_microbial_fasta_arguments` → `_add_microbial_identifier_arguments` (update
  `_build_fasta_download_parser`'s call sites) — behavior unchanged.
- Rename `_resolve_human_fasta_identifiers` → `_resolve_human_identifiers`,
  `_resolve_microbial_fasta_identifiers` → `_resolve_microbial_identifiers` (update
  `download_fasta()`'s call sites) — behavior unchanged.
- Add `_build_domain_download_parser()`: composes the (renamed) human/microbial identifier
  arguments plus `-o/--output_folder` (required).
- Add `_write_domain_table(domains: dict[str, list[str]], output_path: Path) -> None`: writes a
  `Pfam\tEntries\n` header, then one row per pfam_id:
  `f"{pfam_id}\t{';'.join(uniprot_ids) + ';' if uniprot_ids else ''}\n"`.
- Add entry point:
  ```python
  def download_domains() -> int:
      """Download Pfam domains for human and/or microbial proteins."""
      from .workflow import domain_download

      args = _build_domain_download_parser().parse_args()
      human_identifiers = _resolve_human_identifiers(args)
      microbial_identifiers = _resolve_microbial_identifiers(args)

      results = domain_download.download_domains(
          human_identifiers=human_identifiers,
          human_id_type=args.human_id_type,
          microbial_identifiers=microbial_identifiers,
          microbial_id_type=args.microbial_id_type,
      )

      output_folder = Path(args.output_folder)
      output_folder.mkdir(parents=True, exist_ok=True)
      filenames = {'human': 'human_domains.tsv', 'microbial': 'microbial_domains.tsv'}
      for species, domains in results.items():
          _write_domain_table(domains, output_folder / filenames[species])
      return 0
  ```
- `pyproject.toml`: add `microbiolink-download-domains = "microbiolink.cli:download_domains"`
  under `[project.scripts]`. No new dependency required.

## Migration checklist

### 1. Core module
- [ ] Update `microbiolink/workflow/domain_download.py` (already exists, currently returns the old
      protein-keyed shape) so `domain_table_to_mapping` (renamed public, was
      `_domain_table_to_mapping`), `_download_species_domains`, and `download_domains` all build
      and return the Pfam-keyed shape (`pfam_id -> [uniprot_id, ...]`) instead.

### 2. CLI wiring
- [ ] Rename `_add_human_fasta_arguments`/`_add_microbial_fasta_arguments` →
      `_add_human_identifier_arguments`/`_add_microbial_identifier_arguments`.
- [ ] Rename `_resolve_human_fasta_identifiers`/`_resolve_microbial_fasta_identifiers` →
      `_resolve_human_identifiers`/`_resolve_microbial_identifiers`.
- [ ] Add `_build_domain_download_parser`, `_write_domain_table`, `download_domains()` entry point.
- [ ] Add `microbiolink-download-domains = "microbiolink.cli:download_domains"` to
      `pyproject.toml`'s `[project.scripts]`.

### 3. Regression check (per Q11)
- [ ] **Bacterial**: pull `download_protein_list()`/`download_proteome()` logic from
      `workflow/download_bacterial_proteins.py` (top-level, still present). Re-fetch the `Entry`
      IDs from `case_study_input/input/bacterial_protein/BT_BEV_domains.tsv` through both the
      legacy script and the new `download_domains(microbial_identifiers=..., microbial_id_type='uniprot')`.
      Since the legacy fixture is protein-keyed (`Entry`/`Pfam`) and the new output is Pfam-keyed,
      build an inverted view of the legacy fixture (`pfam -> [entries]`) and diff that against the
      new dict, rather than comparing files directly.
- [ ] **Human**: no case-study equivalent exists (case-study has no human-domain download). Live
      comparison against beta's `microbiolink/download_human_domains.py` is best-effort only, same
      MyGene.info-drift caveat Module 3 already hit.
- [ ] Confirm `human_domains.tsv`/`microbial_domains.tsv` header + row format is the new
      `Pfam`/`Entries` convention (no longer expected to match `BT_BEV_domains.tsv`'s literal
      `Entry`/`Pfam` layout, since the output direction changed).

### 4. Delete old code (per Q11)
- [ ] Delete `workflow/download_bacterial_proteins.py` (top-level, pre-refactor) once the
      regression check above passes.
- [ ] Confirm no other code still imports from it (`workflow/DMI.py` only reads a TSV file path via
      argparse, not a Python import, so it's unaffected).

### 5. Manual sign-off
- [ ] Run the new CLI once with only `--human_*` args, once with only `--microbial_*` args, and
      once with both, inspecting `human_domains.tsv`/`microbial_domains.tsv` by eye each time.
- [ ] Spot-check a couple of `Entry` rows against UniProt's website (Pfam section) for correctness.

## Verification

- `uv run ruff format .`
- `uv run ty check`
- Manual CLI runs as in the sign-off checklist above (requires network access to
  `rest.uniprot.org`).
