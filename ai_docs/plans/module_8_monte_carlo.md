# Module 8 — Monte Carlo Simulation — Implementation Plan

Detailed plan for the eighth module described in @ai_docs/plans/microbiolink_refactoring.md, following the
decisions in @ai_docs/plans/refactoring_questions.md and the precedent set by
@ai_docs/plans/module_1_zscore_filter.md through @ai_docs/plans/module_7_idr.md. Scope is exactly Module 8: given a
DMI table (Module 7's IDR-filtered output, **or** Module 6's raw DMI output if Module 7 was skipped)
and the FASTA sequences its motifs came from, test each motif with a Monte Carlo **sequence shuffle**
and keep only motifs that occur more often than the protein's amino-acid composition predicts.

## Context

Modules 1–7 are implemented in `microbiolink/workflow/` (`zscore_filter.py`, `membrane_filter.py`,
`fasta_download.py`, `domain_download.py`, `ddi.py`, `dmi.py`, `idr_filter.py`),
`microbiolink/utils/` (`uniprot_client.py`, `id_resolution.py`, `fasta.py`), and
`microbiolink/data/` (Module 5/6 resource TSVs), all wired into `microbiolink/cli.py`.

**The null this module tests.** A motif's presence is compared against its own protein's amino-acid
composition: *"is this motif more frequent in this protein than chance would predict, given the
protein's residue makeup?"* For one motif on one protein:

1. Look up the motif class's **regex** (the same pattern Module 6 matched with).
2. Count the motif's **observed occurrences** — `re.finditer` matches of that regex in the whole
   motif-bearing protein sequence.
3. **Shuffle**: repeatedly permute the protein's residues (a composition-preserving shuffle — the
   multiset of amino acids is fixed, only their order is destroyed) and re-count regex matches.
4. The **p-value** is `(hits + 1) / (iterations + 1)`, where a hit is a shuffle whose match count is
   `>=` the observed count. A motif **passes** if its p-value is `<= alpha`.

"Monte Carlo" here denotes this motif-shuffle over-representation test. Unlike the earlier
disorder-relocation design, **it needs no disorder or binding profile** — nothing calls IUPred /
AIUPred, and there is no `method`, threshold, device, or `idr` dependency in Module 8. This is the
main way the data processing changes versus Module 7.

**Input — Module 7's table or Module 6's table.** Module 8 accepts either:
- Module 7's output (`idr_filter.filter_by_disorder`): the eight Module 6 columns (`dmi_type`,
  `bacterial_uniprot_id`, `bacterial_annotation`, `human_uniprot_id`, `human_annotation`, `start`,
  `end`, `resource`) plus `disordered_score`, `binding_score`, `combined_score`; **or**
- Module 6's output (`dmi.predict_domain_motif_interactions`) directly, when the IDR step is
  skipped: just the eight columns.

Module 8 reads only the motif-side identifier columns; every other column (including Module 7's three
IDR scores, when present) is passed through untouched. So one code path serves both inputs; the only
difference is how many columns are carried through to the output.

**What Module 8 reads from a row.** The motif-bearing side is `human` if `dmi_type == 'forward'`,
`bacterial` if `'reverse'`. From that side it needs the **protein accession** (`human_uniprot_id` /
`bacterial_uniprot_id`) to fetch the sequence, and the **motif class id(s)** (`human_annotation` /
`bacterial_annotation`) to resolve the regex. `start`/`end` are **not** used — the over-representation
test counts matches across the whole protein, so the motif's specific position is irrelevant to the
null.

**Where the regex comes from.** The DMI table stores the motif *id* in the annotation column, not the
regex. The regexes live in the packaged ELM / 3did resources that `dmi.py` already loads via
`_load_dmi_resources()`. Module 8 resolves ids to patterns through a small **new public accessor**
added to `dmi.py` (`load_motif_regexes()`), rather than reaching into that private loader.

**Ground-truth source.** A prior standalone Monte Carlo script exists on the `development` branch
(`microbiolink/motif_monte_carlo_filter.py`), but it implements the *disorder-relocation* null and is
therefore **not** the reference for this module. Module 8 implements the motif-shuffle null described
above from scratch; the dev-branch script is not ported. It is on `development` only and is not
present in the `refactoring` working tree, so there is no legacy file to delete here (unlike Modules
1–7's per-module deletions under Q11).

## Shared code moves to `utils` (prerequisite refactor of Module 7)

Module 7 (`idr_filter.py`) owns a few generic helpers that Module 8 needs identically. Rather than
have Module 8 reach into `idr_filter`'s privates, the **disorder-free** shared pieces move into
`utils` and both modules import them from there (the Q5/Q6 "shared internal utility, single
implementation" rule). The disorder/binding profile helpers (`_iupred_profile`, `_aiupred_profile`,
`_cached_profile`) are Module-7-only and **stay in `idr_filter.py`** — Module 8 has no use for them.

What moves, and where:

- **`microbiolink/utils/fasta.py` (extend)** — add `build_sequence_lookup(sequences)` (moved from
  `idr_filter._build_sequence_lookup`): reindex a header-keyed FASTA dict to
  `{uniprot_accession: sequence}` via `extract_uniprot_id`, `{}` if `sequences is None`. A pure
  FASTA-dict transform, so `fasta.py` is its natural home.
- **`microbiolink/utils/dmi_table.py` (new)** — helpers that operate on a Module 6/7 interaction row
  or table, generic to any module consuming that shape:
  - `motif_side(row)` — `(motif_uniprot_id, start, end)`, routed by `dmi_type` (moved from
    `idr_filter._motif_side`). Module 7 uses all three fields; Module 8 uses only the accession.
  - `motif_class_ids(row)` — the motif class id(s) on the motif-bearing side
    (`human_annotation` if `dmi_type == 'forward'`, else `bacterial_annotation`), split on `'|'`
    into a list (Module 6 merges co-located classes with `'|'`). New; Module 8-specific but belongs
    with the other row-routing helpers.
  - `validate_direction_sequences(table, human_sequences, bacterial_sequences)` — raise `ValueError`
    if the table has `'forward'` rows but no `human_sequences` (symmetric for `'reverse'` /
    `bacterial_sequences`). This is the sequence-presence half of `idr_filter._validate_inputs`,
    which only reads `dmi_type` and so applies unchanged to both Module 6 and Module 7 tables.

**`idr_filter.py` is then updated** to import `fasta.build_sequence_lookup` and `dmi_table.motif_side`
/ `dmi_table.validate_direction_sequences`, and to delete its now-duplicated `_build_sequence_lookup`
/ `_motif_side` and the sequence-presence half of `_validate_inputs`. Its `method` validation (the
other half of `_validate_inputs`) stays inline in `idr_filter` — it is disorder-specific. Everything
else in Module 7 (`OUTPUT_COLUMNS`, the profile helpers, `_score_row`, `_to_output_row`,
`_score_all_rows`, `filter_by_disorder`) is unchanged. The refactor must be behaviour-preserving —
re-run Module 7's AIUPred case-study regression as a check.

## `dmi.py` — new public regex accessor

Add one small public function so Module 8 can resolve motif ids to regexes without touching the
private loader:

```python
def load_motif_regexes() -> dict[str, str]:
    """Return a motif_id -> regex mapping for every packaged motif class."""
    return {
        motif_id: resource.regex
        for motif_id, resource in _load_dmi_resources().items()
    }
```

It reuses the existing `lru_cache`'d `_load_dmi_resources()`, so it is cheap to call and stays the
single source of truth for the packaged resources.

## Decisions

- **The null is motif over-representation (composition shuffle), not disorder relocation.** Module 8
  destroys residue *order* while preserving residue *composition*, then re-matches the motif regex.
  No disorder/binding profile, no `method`, no thresholds, no `require_binding`, no device flags.
- **Drop failing rows (do not annotate-and-keep).** Like `filter_by_disorder`, Module 8 returns only
  passing rows. The three new columns (`monte_carlo_hits`, `monte_carlo_pvalue`,
  `passes_monte_carlo`) are still added, but `passes_monte_carlo` is `True` on every returned row by
  construction.
- **Rows with several `'|'`-joined motif classes take the minimum p-value.** A merged row can list
  more than one motif class at the same site. Each class's regex is tested independently against the
  same protein; the row's reported `monte_carlo_hits`/`monte_carlo_pvalue` are those of the class
  with the smallest p-value, and the row passes iff that minimum p-value `<= alpha` (the site is
  over-represented under at least its best-supported class). Single-class rows — the common case —
  are unaffected. (Alternative considered: require *all* listed classes to pass, i.e. the maximum
  p-value; rejected as over-strict since the co-located classes are usually variants of one pattern.)

The remaining decisions follow the refactoring precedence rules and Module 7's precedent:

1. **One public function**, `filter_by_monte_carlo`, mirroring `filter_by_disorder`'s "take the
   upstream table shape directly, no dataclasses" style:
   ```python
   def filter_by_monte_carlo(
       interaction_table: pd.DataFrame,
       human_sequences: dict[str, str] | None,
       bacterial_sequences: dict[str, str] | None,
       iterations: int = 1000,
       alpha: float = 0.05,
       seed: int = 0,
   ) -> pd.DataFrame
   ```
   `interaction_table` is Module 7's **or** Module 6's output. `human_sequences`/`bacterial_sequences`
   are the same header-keyed FASTA dicts every downstream module uses. No file I/O in the workflow
   module (file reading stays in `cli.py`).
2. **Output columns are the input's columns plus the three Monte Carlo columns** —
   `[*interaction_table.columns, 'monte_carlo_hits', 'monte_carlo_pvalue', 'passes_monte_carlo']`.
   This is what makes one code path serve both an 8-column Module 6 table and an 11-column Module 7
   table: whatever came in is passed through, three columns are appended. Passing rows are built as
   `row.to_dict()` updated with the three new keys, then assembled with
   `pd.DataFrame(records, columns=output_columns)` (explicit columns so an empty result keeps the
   schema, matching the other modules).
3. **Routing, lookup, and validation come from `utils`**; regexes come from `dmi.load_motif_regexes`.
   Module 8 imports `utils.fasta`, `utils.dmi_table`, and `workflow.dmi` — never `idr_filter`.
4. **A row is dropped, not raised, when it cannot be tested**: its motif protein has no sequence in
   the corresponding lookup, or none of its motif class ids resolve to a known regex (mirrors Module
   7's per-row skip). A wholesale-missing `human_sequences`/`bacterial_sequences` for a *present*
   direction is a hard `ValueError`, via `dmi_table.validate_direction_sequences`.
5. **The shuffle test is cached per `(regex, sequence)`.** Module 6's cross-join means the same
   `(motif class, protein)` recurs across many rows; a `functools.lru_cache`'d test runs each unique
   pair's `iterations` shuffles only once. Reproducibility is preserved by deriving that pair's RNG
   deterministically from `seed` plus a stable hash of the regex and sequence
   (`np.random.default_rng([seed, crc32(regex), crc32(sequence)])`), so every `(seed, input)` gives
   identical results while each protein/motif gets an independent random stream and caching stays
   sound.
6. **No new dependency, and the `idr` extra is not required.** `numpy` (base) supplies the
   permutation; `re` is stdlib. Module 8 runs without any IUPred/AIUPred install — a deliberate
   difference from Module 7.

## Target implementation

### `microbiolink/utils/fasta.py` (extend)

```python
def build_sequence_lookup(sequences: dict[str, str] | None) -> dict[str, str]:
    """Reindex a header-keyed FASTA dict by UniProt accession.

    Args:
        sequences: FASTA header -> sequence mapping, or None.

    Returns:
        A dict mapping each UniProt accession to its sequence. Empty if sequences is None.
    """
```

### `microbiolink/utils/dmi_table.py` (new)

```python
def motif_side(row: pd.Series) -> tuple[str, int, int]:
    """Resolve which protein and position carries the motif for one interaction row.

    Returns (motif_uniprot_id, start, end): the human side if row['dmi_type'] == 'forward',
    otherwise the bacterial side.
    """


def motif_class_ids(row: pd.Series) -> list[str]:
    """Motif class id(s) on the motif-bearing side, split on '|'.

    Reads human_annotation if row['dmi_type'] == 'forward', else bacterial_annotation, and
    splits Module 6's '|'-merged class ids into a list.
    """


def validate_direction_sequences(
    table: pd.DataFrame,
    human_sequences: dict[str, str] | None,
    bacterial_sequences: dict[str, str] | None,
) -> None:
    """Check that sequences are supplied for whichever DMI directions the table contains.

    Raises:
        ValueError: If the table has 'forward' rows with human_sequences=None, or 'reverse'
            rows with bacterial_sequences=None.
    """
```

### `microbiolink/workflow/dmi.py` (edit — add public accessor)

- Add `load_motif_regexes() -> dict[str, str]` (shown above). No other changes.

### `microbiolink/workflow/idr_filter.py` (edit — adopt the utils, delete duplicates)

- Delete `_build_sequence_lookup`, `_motif_side`, and the sequence-presence checks of
  `_validate_inputs`.
- Import `from ..utils import dmi_table` (and keep `fasta`); use `fasta.build_sequence_lookup`,
  `dmi_table.motif_side`, `dmi_table.validate_direction_sequences`.
- In `filter_by_disorder`, replace the `_validate_inputs(...)` call with the inline `method` check
  (kept local) + `dmi_table.validate_direction_sequences(dmi_table, human_sequences,
  bacterial_sequences)`.
- The profile helpers (`_iupred_profile`, `_aiupred_profile`, `_cached_profile`), `OUTPUT_COLUMNS`,
  `_score_row`, `_to_output_row`, `_score_all_rows`, `filter_by_disorder` keep their current
  behaviour.

### `microbiolink/workflow/monte_carlo.py` (new, core, no argparse)

```python
"""Monte Carlo over-representation filtering of domain-motif interactions by sequence shuffling."""

import functools
import re
import zlib

import numpy as np
import pandas as pd

from ..utils import dmi_table, fasta
from . import dmi

MONTE_CARLO_COLUMNS = ["monte_carlo_hits", "monte_carlo_pvalue", "passes_monte_carlo"]


def _count_matches(regex: str, sequence: str) -> int:
    """Number of re.finditer matches of regex in sequence (same matching Module 6 uses)."""


def _shuffle_sequence(sequence: str, rng: np.random.Generator) -> str:
    """A composition-preserving permutation of sequence's residues."""


def _stable_hash(text: str) -> int:
    """Deterministic, run-stable integer hash of text (zlib.crc32), for RNG seeding."""


@functools.lru_cache(maxsize=None)
def _cached_shuffle_test(
    regex: str, sequence: str, iterations: int, seed: int
) -> tuple[int, int, float]:
    """(observed_count, hits, pvalue) for one (regex, sequence).

    observed_count is the regex's match count in sequence; a hit is a shuffle whose match
    count is >= observed_count; pvalue = (hits + 1) / (iterations + 1). The RNG is derived from
    seed and a stable hash of (regex, sequence), so results are reproducible and independent per
    pair, and caching returns each pair's result after a single run.
    """


def _score_row(
    row, sequence_lookup, regexes, iterations, alpha, seed
) -> tuple[bool, int, float] | None:
    """Run the shuffle test for one row's motif class(es).

    Resolves the motif protein via dmi_table.motif_side and its class ids via
    dmi_table.motif_class_ids; tests each resolvable class's regex against the protein sequence
    and keeps the smallest-p-value result. Returns (passes, monte_carlo_hits, monte_carlo_pvalue),
    or None (row dropped) if the protein sequence is missing from sequence_lookup or no class id
    resolves to a known regex. passes is (pvalue <= alpha).
    """


def _filter_all_rows(
    interaction_table, human_lookup, bacterial_lookup, regexes, iterations, alpha, seed
) -> list[dict]:
    """Run the test on every row; collect passing rows as row.to_dict() records with the three
    Monte Carlo keys added. Picks human_lookup for 'forward' rows, bacterial_lookup otherwise."""


def filter_by_monte_carlo(
    interaction_table: pd.DataFrame,
    human_sequences: dict[str, str] | None,
    bacterial_sequences: dict[str, str] | None,
    iterations: int = 1000,
    alpha: float = 0.05,
    seed: int = 0,
) -> pd.DataFrame:
    """Filter a Module 6/7 DMI table by a Monte Carlo motif over-representation test.

    For each row, the motif class's regex is counted in its motif-bearing protein, then re-counted
    over `iterations` composition-preserving shuffles of that protein; the motif passes if its
    p-value (= (hits + 1) / (iterations + 1)) is <= alpha. Tests whether the motif is more frequent
    than the protein's amino-acid composition predicts. No disorder/binding profile is computed.

    Args:
        interaction_table: Module 7's IDR-filtered table, or Module 6's DMI table if IDR was
            skipped. Any extra columns are passed through.
        human_sequences: Human FASTA header -> sequence mapping. Required for 'forward' rows.
        bacterial_sequences: Bacterial FASTA header -> sequence mapping. Required for 'reverse'.
        iterations: Number of shuffles per motif.
        alpha: Maximum p-value a motif may have and still pass.
        seed: Random seed for reproducible sampling.

    Returns:
        interaction_table restricted to rows passing the test, with monte_carlo_hits,
        monte_carlo_pvalue, and passes_monte_carlo columns appended.

    Raises:
        ValueError: If a present direction's sequences are missing.
    """
```

- `filter_by_monte_carlo`: `dmi_table.validate_direction_sequences(interaction_table,
  human_sequences, bacterial_sequences)`; `human_lookup = fasta.build_sequence_lookup(human_sequences)`
  (and bacterial); `regexes = dmi.load_motif_regexes()`; `records = _filter_all_rows(...)`;
  `output_columns = [*interaction_table.columns, *MONTE_CARLO_COLUMNS]`;
  `return pd.DataFrame(records, columns=output_columns)`.
- `_score_row`: `protein_id, _, _ = dmi_table.motif_side(row)`;
  `sequence = sequence_lookup.get(protein_id)`, `None` if missing; `resolved = [regexes[mid] for mid
  in dmi_table.motif_class_ids(row) if mid in regexes]`, `None` if empty; take
  `min((_cached_shuffle_test(regex, sequence, iterations, seed) for regex in resolved),
  key=lambda result: result[2])`; `passes = pvalue <= alpha`.
- `_filter_all_rows`: `for _, row in interaction_table.iterrows()`, pick `human_lookup if
  row['dmi_type'] == 'forward' else bacterial_lookup`, call `_score_row`; skip `None`; when
  `passes`, append `{**row.to_dict(), 'monte_carlo_hits': hits, 'monte_carlo_pvalue': pvalue,
  'passes_monte_carlo': True}`.

### `microbiolink/cli.py` (extend, argparse only)

- `_add_monte_carlo_source_arguments(parser)`: `--interaction_file` (required — Module 7 **or**
  Module 6 output CSV), `--human_fasta_file` (optional), `--bacterial_fasta_file` (optional).
- `_add_monte_carlo_test_arguments(parser)`: `--iterations` (`int`, `1000`), `--alpha` (`float`,
  `0.05`), `--seed` (`int`, `0`).
- `_build_monte_carlo_parser()`: both helpers plus `-o/--output_file` (required).
- `monte_carlo_filter()` entry point: `pd.read_csv(args.interaction_file)`, read the two optional
  FASTA files via `fasta.read_fasta_sequences`, call `monte_carlo.filter_by_monte_carlo(...)`,
  `result.to_csv(args.output_file, index=False)`, `return 0`.

### `pyproject.toml`

- Add `microbiolink-monte-carlo = "microbiolink.cli:monte_carlo_filter"` under `[project.scripts]`.
- No dependency changes, and **no `idr` extra** — Module 8 needs neither.

## Migration checklist

### 1. Shared-code extraction (Module 7 refactor)
- [ ] Add `build_sequence_lookup` to `microbiolink/utils/fasta.py`.
- [ ] Create `microbiolink/utils/dmi_table.py` with `motif_side`, `motif_class_ids`,
      `validate_direction_sequences`.
- [ ] Update `idr_filter.py` to import from utils and delete its duplicated `_build_sequence_lookup`
      / `_motif_side` and the sequence-presence half of `_validate_inputs`.
- [ ] Re-run the Module 7 AIUPred case-study regression from the Module 7 plan to confirm the
      refactor is behaviour-preserving (same passing-row set / scores as before).

### 2. Resource accessor + packaging
- [ ] Add `load_motif_regexes()` to `microbiolink/workflow/dmi.py`.
- [ ] Add `microbiolink-monte-carlo = "microbiolink.cli:monte_carlo_filter"` to `[project.scripts]`.

### 3. Core module
- [ ] Create `microbiolink/workflow/monte_carlo.py` with `_count_matches`, `_shuffle_sequence`,
      `_stable_hash`, `_cached_shuffle_test`, `_score_row`, `_filter_all_rows`,
      `filter_by_monte_carlo`, importing only from `utils` and `workflow.dmi`.
- [ ] Confirm the p-value is `(hits + 1) / (iterations + 1)` with a hit being a shuffle whose match
      count is `>=` observed, and that the shuffle preserves amino-acid composition.
- [ ] Confirm a missing-sequence / unresolvable-motif row is dropped, while a wholesale-missing
      sequences dict for a present direction raises `ValueError`.
- [ ] Confirm output columns = input columns + the three Monte Carlo columns, for **both** a Module 6
      (8-column) and a Module 7 (11-column) input table.

### 4. CLI wiring
- [ ] Add `_add_monte_carlo_source_arguments`, `_add_monte_carlo_test_arguments`,
      `_build_monte_carlo_parser`, `monte_carlo_filter()` to `cli.py`.

### 5. Delete old code (per Q11)
- [ ] None. No Monte Carlo file exists in the `refactoring` working tree (the `development` script
      implements a different null and is not ported).

### 6. Verification — no committed fixture (manual sign-off, as Modules 5/6 and Module 7's IUPred path)
- [ ] Confirm there is no `case_study_output` Monte Carlo fixture (checked: none exists).
- [ ] Deterministic seeded unit check on synthetic sequences: a protein where the motif is repeated
      far more than its composition predicts passes (low p-value); a protein where the motif occurs
      about as often as chance fails (high p-value); confirm same-`seed` runs are identical and a
      different `seed` perturbs `monte_carlo_hits` but not the pass/fail conclusion.
- [ ] Run end-to-end on a real Module 7 output (the AIUPred-filtered table) **and** on a raw Module 6
      table with `microbiolink-monte-carlo`: confirm the appended-3-column schema in both, output
      rows ⊆ input rows, `passes_monte_carlo` all-`True`, plausible pass counts.
- [ ] Verify `dmi_type` routing directly with a synthetic 2-row table (one `forward`, one `reverse`)
      by swapping motif-rich vs. motif-poor synthetic sequences between the human and bacterial
      dicts.
- [ ] Shown to the user for confirmation.

## Verification

- `uv run ruff format` / `uv run ruff check` scoped to the new files (`utils/dmi_table.py`,
  `workflow/monte_carlo.py`); `fasta.py`, `dmi.py`, `idr_filter.py`, `cli.py`, and `pyproject.toml`
  edits reviewed by eye against surrounding style (modified, not new).
- `uv run ty check`
- The Module 7 regression re-run (step 1) plus the seeded/synthetic checks above.
