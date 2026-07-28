# Module 6 — Domain-Motif Interactions (DMI) — Implementation Plan

Detailed plan for the sixth module described in @ai_docs/plans/microbiolink_refactoring.md, following the
decisions in @ai_docs/plans/refactoring_questions.md and the precedent set by
@ai_docs/plans/module_1_zscore_filter.md through @ai_docs/plans/module_5_ddi.md. Scope is exactly Module 6: given
a human or bacterial FASTA file (Module 3's output) and the partner species' Pfam domain
dictionary (Module 4's output shape), predict domain-motif interactions using the ELM and 3did
resources, in the forward direction (motif in a human sequence, matched against bacterial
domains), the reverse direction (motif in a bacterial sequence, matched against human domains), or
both.

## Context

Modules 1–5 are implemented in `microbiolink/workflow/` (`zscore_filter.py`, `membrane_filter.py`,
`fasta_download.py`, `domain_download.py`, `ddi.py`), `microbiolink/utils/` (`uniprot_client.py`,
`id_resolution.py`), and `microbiolink/data/` (the Module 5 DDI resource TSVs), all wired into
`microbiolink/cli.py`. Module 6 has **no implementation on `case-study`** beyond a forward-only
legacy script, `workflow/DMI.py`, which is not packaged and mixes argparse with logic — it's used
below as a fixture source, not ported directly. The only working *bidirectional* implementation
anywhere is `microbiolink_api/dmi.py` on `origin/MicrobioLink-2.1-beta`, which needs to be ported
in and simplified to match this repo's established plain-function/dict/DataFrame style (Modules
1–5 use no dataclasses).

Key facts from beta's `microbiolink_api/dmi.py`:

- Core algorithm: for each motif class, regex-search it across every sequence in one species'
  FASTA (`re.finditer`), then for each Pfam domain known to interact with that motif class, cross-
  join every motif hit against every partner-species protein carrying that domain. Forward mode
  searches motifs in **human** sequences and matches them to **bacterial** domains; reverse mode
  searches motifs in **bacterial** sequences and matches them to **human** domains, excluding ELM
  motif classes whose name starts with `CLV_` (cleavage-site motifs — eukaryotic-protease-specific,
  not meaningful on bacterial sequences).
- Resources are four packaged TSVs: `elm_classes.tsv` / `elm_interaction_domains.tsv` (ELM) and
  `3did_dmi_classes.tsv` / `3did_dmi_interaction_domains.tsv` (3did structural motifs). Beta's
  default merges ELM + 3did, matching this module's spec text ("Known DMIs from ... ELM and 3did").
  Both `*_classes.tsv` files use the same column layout — motif ID at index 1, regex at index 4 —
  and both `*_interaction_domains.tsv` files use motif ID at index 0, Pfam ID at index 1, so one
  pair of parsing functions covers both resources (mirrors Module 5's single
  `_read_pfam_pair_table` covering both its resource files).
- **Bug found while inspecting the packaged data**: beta's `_parse_elm_regex_lines` skips exactly
  one leading line (`next(iterator, None)`) plus any further `#`-prefixed comment lines, then treats
  every remaining line as data. `elm_classes.tsv` has a 5-line `#`-prefixed comment preamble
  *followed by* a real tab-separated header row (`"Accession"	"ELMIdentifier"	...`) that does not
  start with `#`. Beta's parser skips the 5 comment lines correctly but then parses the header row
  itself as a data row, inserting one spurious entry (`elm_regex['ELMIdentifier'] = 'Regex'`) into
  every resource bundle it loads. It's functionally harmless (a literal-string "Regex" pattern never
  matches an uppercase amino-acid sequence, and the phantom key never appears in
  `elm_interaction_domains.tsv`, so it produces zero interactions) but is fixed in this port rather
  than carried forward, per decision 6 below.
- Beta's input shapes: FASTA sequences as `dict[fasta_header -> sequence]` (via
  `read_fasta_sequences`, header-keyed, not yet resolved to UniProt accessions — the UniProt
  accession is extracted per-match via `extract_uniprot_id`, which splits `header` on `'|'` and
  takes index 1, matching UniProt's `>sp|ACCESSION|...` FASTA header convention that
  `uniprot_client.fetch_fasta_sequences` already produces) and domains as `dict[pfam_id ->
  [uniprot_id, ...]]` — identical to Module 5, i.e. Module 4's output shape used as-is, no reshaping.
- No `case_study_output` fixture exists for DMI (case-study's `workflow/DMI.py` is forward-only,
  standalone, and was never run to produce a committed output file). Per Q11 this counts as new
  functionality requiring manual sign-off rather than a diff-based regression check, same as Module
  5. Real fixtures exist to exercise it end-to-end though:
  `case_study_input/input/human_transcriptomics/protein_sequences.fasta` (real human FASTA, correct
  `>sp|ACCESSION|...` header shape) for the forward-mode motif side, and
  `case_study_input/input/bacterial_protein/BT_BEV_domains.tsv` (same fixture Module 5 used) as an
  `Entry` ID source to re-fetch bacterial domains and, for reverse mode, bacterial FASTA, through
  Module 3/4's CLIs.

Two scope decisions were confirmed with the user before writing this plan:

1. **Keep the `CLV_`-prefix exclusion for reverse mode** — ported faithfully from beta even though
   the top-level plan text doesn't mention it, since it's a real biological constraint (cleavage
   motifs are protease-specific and not meaningful when searched in a bacterial sequence).
2. **Output keeps 3 extra columns beyond the plan's literal 5** — `start`, `end` (the motif match's
   position in whichever protein carries the motif) and `resource` (`'ELM'` or `'3did'`) — cheap to
   retain and useful for downstream filtering (Module 7 needs sequence positions for IDR scoring),
   mirroring Module 5's identical call to keep an extra `resource` column.

## Decisions

1. **No dataclasses.** Beta's `DomainMotifInteraction` / `BidirectionalDomainMotifInteraction` /
   `DMIResourceBundle` are dropped. Resources are two plain dicts (`motif_regex`, `motif_domains`)
   plus a `motif_sources` dict (`motif_id -> 'ELM'/'3did'`) — matching Modules 1–5's plain
   dict/DataFrame style, and Module 5's precedent of dropping beta's equivalent
   `DDIResourceBundle`.
2. **Core module accepts Module 3's FASTA-parse shape and Module 4's dict shape directly.** A
   single public function, `predict_domain_motif_interactions(mode, human_sequences=None,
   bacterial_domains=None, bacterial_sequences=None, human_domains=None)`, taking
   `dict[fasta_header -> sequence]` for whichever FASTA side(s) `mode` needs and
   `dict[pfam_id -> [uniprot_id, ...]]` for whichever domain side(s) `mode` needs. No file I/O in
   `workflow/dmi.py` — `cli.py` reads the FASTA and TSV files and passes in the parsed dicts,
   matching Module 5's precedent.
3. **`mode` controls which inputs are required, mirroring beta's `mode` parameter** (`'forward'`,
   `'reverse'`, `'both'`) rather than two separate public functions. `'forward'` requires
   `human_sequences` + `bacterial_domains`; `'reverse'` requires `bacterial_sequences` +
   `human_domains`; `'both'` requires all four. Missing required inputs for the requested mode raise
   `ValueError` (no custom exception classes — matches `id_resolution.resolve_uniprot_ids`'s
   precedent of `ValueError` over beta's `InputFormatError`, which doesn't exist in this package).
4. **New `microbiolink/utils/fasta.py`** for `read_fasta_sequences(filename) ->
   dict[str, str]` and `extract_uniprot_id(fasta_header: str) -> str`, ported from beta's
   `dmi.py`. This is generic FASTA I/O, not DMI-specific logic, so it goes in `utils/` (which
   "has functions that repeat across modules" per the top-level plan's package layout) rather than
   `workflow/dmi.py` — Module 7 (IDR prediction) will need protein sequences too, so this is written
   to be reused, not duplicated, the same reasoning Q5/Q6 applied to the shared translation/fetch
   utilities.
5. **Extend `microbiolink/data/` with the four ELM/3did resource files**, following Module 5's
   packaged-TSV + `importlib.resources` pattern (Q9's corollary — this already-established location
   is reused rather than inventing a second resource convention for Module 6).
6. **Fix beta's header-row parsing bug during the port** (see Context) rather than carry it
   forward: the shared `_read_motif_regex_table` / `_read_motif_domain_table` helpers skip leading
   `#`-prefixed comment lines *and then* unconditionally skip exactly one header row, instead of
   beta's "skip one line, then skip further comment lines" order which misses a header row that
   comes after the comments. This fixes `elm_classes.tsv` (5 comment lines + 1 header row) while
   still working correctly for `3did_dmi_classes.tsv` and both `*_interaction_domains.tsv` files
   (0 comment lines + 1 header row each).
7. **`_find_motif_matches` keeps taking header-keyed FASTA dicts, not pre-resolved UniProt-ID-keyed
   dicts**, resolving the accession per match via `utils.fasta.extract_uniprot_id` — matches beta's
   contract exactly and means the CLI only has to call `fasta.read_fasta_sequences(path)` with no
   extra reshaping step before passing sequences in.
8. **Post-hoc merge of same-position matches across resources, not exact-row dedup.** ELM and 3did
   motif ID namespaces don't overlap (ELM IDs look like `LIG_SH3_1`; 3did IDs are prefixed
   `3DID_...`), so `motif_sources` in this port is a single flat `motif_id -> str` mapping (no merge
   conflict when loading resources — beta's own merge logic resolves each motif to exactly one
   source too, via last-write-wins with no conflict check, unlike the conflict check it does perform
   on regex patterns). But that means the *same* biological site, independently matched by an ELM
   class and a 3did class, produces two structurally different rows (different motif-side
   `annotation`, different `resource`) that a plain `seen`-set dedup would never catch. Two rows are
   treated as the same underlying interaction if they agree on `dmi_type`, both protein IDs, the
   domain-side `annotation` (already a shared Pfam ID, needs no merging), and the motif's
   `start`/`end` — an exact position match across independently-curated resources is strong evidence
   both are flagging the same real site. `_merge_duplicate_matches` groups rows on that key and
   pipe-joins the motif-side `annotation` and `resource` values (deduplicated, order-preserving)
   into one row, mirroring Module 5's `resource='3did|DOMINE_hc'` merge for Pfam pairs found in both
   DDI resources — except here only the motif-side `annotation` and `resource` are merged, not the
   whole row, since the domain-side annotation is already resource-agnostic.
9. **Output columns**: `dmi_type` (`'forward'`/`'reverse'`), `bacterial_uniprot_id`,
   `bacterial_annotation`, `human_uniprot_id`, `human_annotation`, `start`, `end`, `resource` —
   snake_case naming consistent with Module 5's `bacterial_uniprot_id`/`human_uniprot_id`
   convention, rather than beta's `host_protein`/`microbial_protein`/`domain`/`motif` naming. Column
   named `dmi_type` rather than `type` to avoid shadowing the Python builtin and for pandas-column
   clarity.
10. **CLI writes CSV** (`result.to_csv(output_path, index=False)`), matching Modules 2/4/5's output
    convention.
11. **CLI reuses `_read_domain_mapping` unchanged** (already added to `cli.py` for Module 5) for
    both `--bacterial_domain_file` and `--human_domain_file` — no new domain-parsing code needed,
    since Module 6 consumes the exact same Module-4-CLI-output TSV shape Module 5 does.
12. **Domain-relevance prefilter, consolidated into one `MotifResource` NamedTuple keyed by motif
    ID.** Scanning every packaged motif class against every sequence regardless of whether its
    compatible domains are even present in that call's `partner_domains` wastes `re.finditer` work
    whenever the curated ELM/3did universe is larger than the partner species' actual Pfam content
    (the common case). A new `_build_motif_index(motif_regex, motif_domains, motif_sources,
    partner_domains) -> dict[str, MotifResource]` does two jobs in one pass: it keeps only motif IDs
    with at least one compatible domain present in `partner_domains` (dropping the rest before any
    regex runs), and consolidates each surviving motif's regex/compatible-domains/source into one
    `MotifResource(regex, domains, source)` NamedTuple. Keyed by **motif ID**, not by regex text —
    regex strings aren't guaranteed unique across ~700+ curated motif classes, and collisions would
    silently merge two different motifs' domain lists. `MotifResource.domains` is itself pre-filtered
    to just the domains present in `partner_domains`, so the row-building loop in
    `_predict_directional_interactions` indexes `partner_domains[domain_name]` directly instead of
    `.get(domain_name, [])`. `_find_motif_matches` and `_predict_directional_interactions` take this
    single `motif_index` in place of the three parallel `motif_regex`/`motif_domains`/`motif_sources`
    dicts described in decision 1. `_build_motif_index` runs per call (once per direction), not
    inside the `lru_cache`-wrapped `_load_dmi_resources`, since it depends on the call-specific
    `partner_domains` input, not just the packaged resource data. `MotifResource` is a `NamedTuple`,
    not a dataclass — a narrow exception to decision 1's "no dataclasses", justified here purely for
    readability of a function that would otherwise take three redundant parallel dicts; it carries no
    methods and behaves like a plain tuple.

## Target implementation

### `microbiolink/data/` (extend existing package, resources only)

```
microbiolink/data/elm_classes.tsv                        # copied from beta via `git show`
microbiolink/data/elm_interaction_domains.tsv             # copied from beta via `git show`
microbiolink/data/3did_dmi_classes.tsv                    # copied from beta via `git show`
microbiolink/data/3did_dmi_interaction_domains.tsv        # copied from beta via `git show`
```

All four copied unmodified from
`origin/MicrobioLink-2.1-beta:microbiolink_api/resources/{elm_classes,elm_interaction_domains,
3did_dmi_classes,3did_dmi_interaction_domains}.tsv`.

### `microbiolink/utils/fasta.py` (new)

```python
def read_fasta_sequences(filename: PathLike) -> dict[str, str]:
    """Read a FASTA file into a header-to-sequence mapping."""


def extract_uniprot_id(fasta_header: str) -> str:
    """Extract the UniProt accession from a FASTA header."""
```

- `read_fasta_sequences`: ported near-verbatim from beta's `dmi.read_fasta_sequences` — streams the
  file, accumulates sequence lines under the current `>`-prefixed header, no external FASTA
  library dependency (matches Module 3's own no-new-dependency fetch/write approach).
- `extract_uniprot_id`: splits on `'|'`, returns index 1; raises `ValueError` (not beta's
  `InputFormatError`, which has no equivalent in this package) if the header has fewer than 2
  `'|'`-separated fields.

### `microbiolink/workflow/dmi.py` (new, core, no argparse)

```python
def _read_motif_regex_table(filename) -> dict[str, str]:
    """Read a motif-class TSV (ELM or 3did) into a motif_id -> regex mapping."""


def _read_motif_domain_table(filename) -> dict[str, list[str]]:
    """Read a motif-domain TSV (ELM or 3did) into a motif_id -> Pfam ID list mapping."""


@functools.lru_cache(maxsize=None)
def _load_dmi_resources() -> tuple[dict[str, str], dict[str, list[str]], dict[str, str]]:
    """Load and merge the packaged ELM and 3did motif-regex/motif-domain resources.

    Returns:
        (motif_regex, motif_domains, motif_sources) — motif_sources maps
        motif_id to 'ELM' or '3did'.
    """


def _filter_reverse_motifs(
    motif_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    motif_sources: dict[str, str],
) -> tuple[dict[str, str], dict[str, list[str]], dict[str, str]]:
    """Drop CLV_-prefixed (cleavage-site) motif classes for reverse-mode DMI."""


class MotifResource(NamedTuple):
    """One motif class's regex, domain-relevance-filtered compatible domains, and source."""

    regex: str
    domains: list[str]
    source: str


def _build_motif_index(
    motif_regex: dict[str, str],
    motif_domains: dict[str, list[str]],
    motif_sources: dict[str, str],
    partner_domains: dict[str, list[str]],
) -> dict[str, MotifResource]:
    """Index motif classes with at least one compatible domain present in partner_domains."""


def _find_motif_matches(
    sequences: dict[str, str],
    motif_index: dict[str, MotifResource],
) -> dict[str, list[tuple[str, int, int]]]:
    """Find every (motif_id, start, end) regex match per UniProt accession."""


def _predict_directional_interactions(
    motif_sequences: dict[str, str],
    motif_index: dict[str, MotifResource],
    partner_domains: dict[str, list[str]],
    motif_is_bacterial: bool,
) -> list[tuple]:
    """Predict one directional set of domain-motif interaction rows."""


def _merge_duplicate_matches(rows: list[tuple]) -> list[tuple]:
    """Merge rows that share protein pair, domain, and motif position across resources."""


def predict_domain_motif_interactions(
    mode: str,
    human_sequences: dict[str, str] | None = None,
    bacterial_domains: dict[str, list[str]] | None = None,
    bacterial_sequences: dict[str, str] | None = None,
    human_domains: dict[str, list[str]] | None = None,
) -> pd.DataFrame:
    """Predict forward and/or reverse domain-motif interactions."""
```

- `_read_motif_regex_table` / `_read_motif_domain_table`: shared by both ELM and 3did files (see
  decision 6) — skip leading `#`-comment lines, skip exactly one header row, then for each
  remaining line, strip `"` quotes and split on `'\t'`; regex table takes `fields[1]` (motif ID) /
  `fields[4]` (regex) when `len(fields) > 4`; domain table takes `fields[0]` (motif ID) /
  `fields[1]` (Pfam ID) when `len(fields) > 1`, appending to `motif_domains.setdefault(...)`.
- `_load_dmi_resources`: loads all four packaged TSVs via
  `importlib.resources.files(data).joinpath(filename)` (imports the sibling `data` package, same as
  `ddi.py`), builds `motif_sources` per resource (`'ELM'` for the ELM pair, `'3did'` for the 3did
  pair) over `set(motif_regex) | set(motif_domains)`, then merges ELM and 3did into single
  `motif_regex` / `motif_domains` / `motif_sources` dicts (later resource's `motif_sources` entries
  win on key collision, matching beta's own merge order — collisions aren't expected in practice
  per decision 8).
- `_filter_reverse_motifs`: returns new dicts filtered to keys not starting with `'CLV_'`.
- `_build_motif_index` (decision 12): for each `motif_id, domains` in `motif_domains`, computes
  `relevant_domains = [d for d in domains if d in partner_domains]`; if empty, skips the motif
  entirely (it can never contribute a row for this call's `partner_domains`). Otherwise builds
  `MotifResource(regex=motif_regex[motif_id], domains=relevant_domains,
  source=motif_sources.get(motif_id, 'ELM'))` and stores it under `motif_id`. Returns the resulting
  `dict[str, MotifResource]`, typically far smaller than the full `motif_domains` universe.
- `_find_motif_matches`: for each `(header, sequence)`, resolves `uniprot_id =
  fasta.extract_uniprot_id(header)`, then for each `(motif_id, resource)` in `motif_index`, appends
  every `re.finditer(resource.regex, sequence)` hit as `(motif_id, match.start(), match.end())`;
  skips the accession entirely if it had zero hits. Only scans the motifs `_build_motif_index` kept,
  not the full packaged resource set.
- `_predict_directional_interactions`: calls `_find_motif_matches(motif_sequences, motif_index)`;
  for each `(uniprot_id, hits)` in the result, for each `(motif_id, start, end)` hit, looks up
  `resource = motif_index[motif_id]` and, for each `domain_name` in `resource.domains` (already
  filtered, so `partner_domains[domain_name]` is indexed directly, no `.get()` needed), for each
  partner protein carrying that Pfam, emits one row. `motif_is_bacterial` picks which side is which:
  `False` → forward (`dmi_type='forward'`, motif protein is `human_uniprot_id`/`human_annotation`,
  partner is `bacterial_uniprot_id`/`bacterial_annotation`); `True` → reverse (motif protein is
  `bacterial_uniprot_id`/`bacterial_annotation`, partner is `human_uniprot_id`/`human_annotation`).
  `resource = resource.source` per row.
- `predict_domain_motif_interactions`: validates `mode in {'forward', 'reverse', 'both'}`
  (`ValueError` otherwise); loads `_load_dmi_resources()` once; for `mode in {'forward', 'both'}`,
  requires `human_sequences` and `bacterial_domains` (`ValueError` if either is `None`), builds
  `motif_index = _build_motif_index(motif_regex, motif_domains, motif_sources, bacterial_domains)`,
  and calls `_predict_directional_interactions(human_sequences, motif_index, bacterial_domains,
  motif_is_bacterial=False)`; for `mode in {'reverse', 'both'}`, requires `bacterial_sequences` and
  `human_domains`, filters resources via `_filter_reverse_motifs`, builds `motif_index =
  _build_motif_index(filtered_regex, filtered_domains, filtered_sources, human_domains)`, and calls
  `_predict_directional_interactions(bacterial_sequences, motif_index, human_domains,
  motif_is_bacterial=True)`. Rows from both directions are combined and passed through
  `_merge_duplicate_matches` (decision 8), then returned as
  `pd.DataFrame(rows, columns=OUTPUT_COLUMNS)` with `OUTPUT_COLUMNS = ['dmi_type',
  'bacterial_uniprot_id', 'bacterial_annotation', 'human_uniprot_id', 'human_annotation', 'start',
  'end', 'resource']` (explicit columns so an empty result still has the right schema).
- `_merge_duplicate_matches`: for each row, picks `domain_annotation` = `bacterial_annotation` if
  `dmi_type == 'forward'` else `human_annotation` (the Pfam-ID side), and `motif_annotation` = the
  other side. Groups rows by `(dmi_type, bacterial_uniprot_id, human_uniprot_id, start, end,
  domain_annotation)`, preserving first-seen group order; within each group, collects
  `motif_annotation` and `resource` values into order-preserving deduplicated lists (matching beta's
  `_unique_preserve_order` helper) and `'|'.join`s each; rebuilds one output row per group with the
  merged, pipe-joined `annotation`/`resource` values on the motif side and the unchanged shared
  Pfam ID on the domain side.

### `microbiolink/cli.py` (extend, argparse only)

- Add `_build_dmi_parser()`: `-m/--mode` (required, `choices=['forward', 'reverse', 'both']`),
  `-hf/--human_fasta_file` (optional — required for forward/both, validated in the workflow call),
  `-b/--bacterial_domain_file` (optional), `-bf/--bacterial_fasta_file` (optional),
  `-hu/--human_domain_file` (optional), `-o/--output_file` (required).
- Reuse `_read_domain_mapping` (already defined for `ddi()`) unchanged for
  `--bacterial_domain_file` / `--human_domain_file`.
- Add entry point:
  ```python
  def dmi() -> int:
      """Predict domain-motif interactions between bacterial and human proteins."""
      from .utils import fasta
      from .workflow import dmi as dmi_workflow

      args = _build_dmi_parser().parse_args()
      human_sequences = fasta.read_fasta_sequences(args.human_fasta_file) if args.human_fasta_file else None
      bacterial_domains = _read_domain_mapping(args.bacterial_domain_file) if args.bacterial_domain_file else None
      bacterial_sequences = fasta.read_fasta_sequences(args.bacterial_fasta_file) if args.bacterial_fasta_file else None
      human_domains = _read_domain_mapping(args.human_domain_file) if args.human_domain_file else None

      result = dmi_workflow.predict_domain_motif_interactions(
          args.mode,
          human_sequences=human_sequences,
          bacterial_domains=bacterial_domains,
          bacterial_sequences=bacterial_sequences,
          human_domains=human_domains,
      )
      result.to_csv(args.output_file, index=False)
      return 0
  ```
- `pyproject.toml`: add `microbiolink-dmi = "microbiolink.cli:dmi"` under `[project.scripts]`. No
  `[tool.hatch.build] include` change needed — it already covers `microbiolink/data/*.tsv`.

## Migration checklist

### 1. Resource files
- [x] Copy `elm_classes.tsv`, `elm_interaction_domains.tsv`, `3did_dmi_classes.tsv`,
      `3did_dmi_interaction_domains.tsv` from `origin/MicrobioLink-2.1-beta:microbiolink_api/resources/`
      into `microbiolink/data/`.

### 2. Utils
- [x] Create `microbiolink/utils/fasta.py` with `read_fasta_sequences`, `extract_uniprot_id`.

### 3. Core module
- [x] Create `microbiolink/workflow/dmi.py` with `_read_motif_regex_table`,
      `_read_motif_domain_table`, `_load_dmi_resources`, `_filter_reverse_motifs`,
      `MotifResource`, `_build_motif_index`, `_find_motif_matches`,
      `_predict_directional_interactions`, `_merge_duplicate_matches`,
      `predict_domain_motif_interactions`.
- [x] Confirm the header-row parsing fix (decision 6) — sanity-check that `_load_dmi_resources()`
      does **not** contain a spurious `'ELMIdentifier'` key (beta's bug, see Context).
- [x] Confirm `_build_motif_index` (decision 12) excludes a motif whose compatible domains don't
      intersect `partner_domains`, and that a surviving `MotifResource.domains` list only contains
      domains actually present as `partner_domains` keys.

### 4. CLI wiring
- [x] Add `_build_dmi_parser`, `dmi()` entry point to `cli.py`, reusing `_read_domain_mapping`.
- [x] Add `microbiolink-dmi = "microbiolink.cli:dmi"` to `pyproject.toml`.

### 5. Manual sign-off (per Q11 — no case-study ground truth exists for this module)
- [x] Confirm a known ELM motif/domain pair (e.g. pick one `motif_id` from
      `elm_interaction_domains.tsv` and confirm its regex from `elm_classes.tsv` and its Pfam ID
      round-trip through `_load_dmi_resources()` with `resource == 'ELM'`).
- [x] Confirm a `CLV_`-prefixed motif is present in `_load_dmi_resources()`'s output but absent
      after `_filter_reverse_motifs()`.
- [x] Confirm `_merge_duplicate_matches` merges two synthetic rows that share `dmi_type`, both
      protein IDs, `start`/`end`, and the domain-side annotation but differ in motif-side
      `annotation`/`resource` into a single row with pipe-joined `annotation`/`resource`, and leaves
      rows with any other differing field unmerged.
- [x] Run forward mode using
      `case_study_input/input/human_transcriptomics/protein_sequences.fasta` directly as
      `--human_fasta_file`, and a bacterial domain file generated by running
      `microbiolink-download-domains` against the `Entry` IDs in
      `case_study_input/input/bacterial_protein/BT_BEV_domains.tsv` as `--bacterial_domain_file`.
- [x] Generate a bacterial FASTA file by running `microbiolink-download-fasta` against the same
      `Entry` IDs, and a human domain file by running `microbiolink-download-domains` against
      `case_study_input/input/human_protein/Enterocyte_Manual/enterocyte_colon_CD_expressed_genes.csv`
      (same fixture Module 5 used), then run reverse mode with `--bacterial_fasta_file` /
      `--human_domain_file`.
- [x] Run `mode=both` with all four inputs supplied and confirm the output contains both
      `dmi_type` values.
- [ ] Show the resulting output table(s) to the user for confirmation before considering Module 6
      done.

## Verification

- `uv run ruff format microbiolink/utils/fasta.py microbiolink/workflow/dmi.py` — scoped to the two
  newly created files, not the whole tree (`ruff format .` would also reformat pre-existing files
  that aren't currently ruff-clean, which is out of scope for this module).
- `uv run ruff check microbiolink/utils/fasta.py microbiolink/workflow/dmi.py` — same scoping.
  `cli.py`'s new `_build_dmi_parser`/`dmi()` additions are reviewed by eye against the existing
  file's style instead, since `cli.py` itself is modified, not newly created, and isn't currently
  ruff-clean in full.
- `uv run ty check`
- Manual CLI runs as in the sign-off checklist above.
