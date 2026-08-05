# Module 8 — Monte Carlo Simulation — Implementation Plan

Detailed plan for the eighth module described in @ai_docs/plans/microbiolink_refactoring.md, following the
decisions in @ai_docs/plans/refactoring_questions.md, the null-model choice recorded in
@ai_docs/decisions/monte_carlo_shuffling_strategy.md, and the precedent set by
@ai_docs/plans/module_1_zscore_filter.md through @ai_docs/plans/module_7_idr.md. Scope is exactly
Module 8: given Module 7's IDR-filtered DMI table and the FASTA sequences its motifs came from, test
each motif with a Monte Carlo **over-representation** test against a **pooled disordered-region
composition**, and keep only motifs that occur more often than that background predicts.

The null model and every alternative it was chosen over are recorded in
@ai_docs/decisions/monte_carlo_shuffling_strategy.md (the chosen approach is "option (a)"). This plan
implements that decision; read the decision document first for the *why*, this plan for the *how*.

## Context

Modules 1–7 are implemented in `microbiolink/workflow/` (`zscore_filter.py`, `membrane_filter.py`,
`fasta_download.py`, `domain_download.py`, `ddi.py`, `dmi.py`, `idr_filter.py`),
`microbiolink/utils/` (`uniprot_client.py`, `id_resolution.py`, `fasta.py`), and
`microbiolink/data/` (Module 5/6 resource TSVs), all wired into `microbiolink/cli.py`.

**The null this module tests (option (a)).** A motif's presence is compared against the amino-acid
composition of *pooled disordered regions*: *"does this motif's regex occur more often in the
protein's disordered region than a random disordered stretch of the same length would, given the
disordered proteome's residue composition?"* For one motif on one protein:

1. Compute the motif-bearing protein's **per-residue disorder profile** and mark residues disordered
   at `disorder_cutoff`. Find the maximal disordered run (**region**) containing the motif; its length
   is `D`.
2. Look up the motif class's **regex** (the same pattern Module 6 matched with) and count its
   **observed occurrences** in the region subsequence.
3. **Sample** `iterations` synthetic length-`D` sequences by drawing residues i.i.d. from the pooled
   disordered-region composition (leave-one-out: the target protein's own residues excluded), and
   re-count regex matches in each.
4. The **p-value** is `(hits + 1) / (iterations + 1)`, where a hit is a sample whose match count is
   `>=` the observed count.

Raw p-values are then **corrected for multiple testing** across all motif instances with
Benjamini–Hochberg (BH), and a motif **passes** if its BH-adjusted **q-value** is `<= alpha` — where
`alpha` is the **target false discovery rate**, not a per-test threshold. See the decision document's
"Multiple-testing correction" section for the rationale; the mechanics are detailed under *Decisions*
below.

This is a composition-matched *sampling* null, not a permutation shuffle: sampling from the pooled
frequency vector is the efficient canonical form of the identical composition-only null (see the
decision document's rejected options 3–4).

**Disorder dependency is required.** Unlike the earlier whole-protein-shuffle design, option (a)
**needs a disorder profile** — both to locate each motif's disordered region and to build the pooled
background. Module 8 computes the profile itself (Module 7 emits per-window *scores*, not per-residue
profiles or region boundaries), reusing Module 7's existing disorder-profile helpers (imported from
`idr_filter`), and therefore takes the same `method` / `disorder_cutoff` / device arguments as Module
7 and requires the same optional `idr` dependency. This reverses the previous "Module 8 needs no disorder install" property; it is
acceptable because Module 8 is a **post-IDR step** (pipeline order M7 → M8) where that dependency is
already present, and the shared profile cache makes re-profiling near-free if Module 7 ran in-process.

**Input — Module 7's table.** Module 8's intended predecessor is Module 7's output
(`idr_filter.filter_by_disorder`): the eight Module 6 columns (`dmi_type`, `bacterial_uniprot_id`,
`bacterial_annotation`, `human_uniprot_id`, `human_annotation`, `start`, `end`, `resource`) plus
`disordered_score`, `binding_score`, `combined_score`. Module 8 reads only the motif-side columns
(accession, `start`/`end`, and the motif class id(s)); every other column — including Module 7's three
IDR scores — is passed through untouched. The `start`/`end` positions **are** used now, to locate the
motif's disordered region. A row whose motif does not lie inside a disordered region is dropped
(shouldn't occur for genuine post-IDR input, but guarded). Because non-motif columns are only passed
through, a wider or narrower input table still works; the appended-columns contract is identical.

**Where the regex comes from.** The DMI table stores the motif *id* in the annotation column, not the
regex. The regexes live in the packaged ELM / 3did resources that `dmi.py` already loads via
`_load_dmi_resources()`. Module 8 resolves ids to patterns through a small **new public accessor**
added to `dmi.py` (`load_motif_regexes()`), rather than reaching into that private loader.

**Ground-truth source.** A prior standalone Monte Carlo script exists on the `development` branch
(`microbiolink/motif_monte_carlo_filter.py`), but it implements the *disorder-relocation* null
(decision document Method 1) and is therefore **not** the reference for this module. Module 8
implements option (a) from scratch; the dev-branch script is not ported. It is on `development` only
and is not present in the `refactoring` working tree, so there is no legacy file to delete here
(unlike Modules 1–7's per-module deletions under Q11).

## Shared code: row helpers move to `utils`, disorder helpers stay in `idr_filter`

Module 7 (`idr_filter.py`) owns helpers that Module 8 now needs identically. The row-routing and
sequence-lookup helpers move into `utils` so both modules import them from one place (the Q5/Q6
"shared internal utility, single implementation" rule). **The disorder-profile helpers stay in
`idr_filter.py`**: Module 8 imports them from there rather than from a new `utils/disorder.py`. This
keeps all disorder logic in its existing home and matches the decision document, which records Module
8 as reusing `idr_filter`'s disorder helpers. (Module 8 already imports one other workflow module,
`workflow.dmi`, so importing `idr_filter` for the profile helpers follows the same precedent.)

What moves to `utils`:

- **`microbiolink/utils/fasta.py` (extend)** — add `build_sequence_lookup(sequences)` (moved from
  `idr_filter._build_sequence_lookup`): reindex a header-keyed FASTA dict to
  `{uniprot_accession: sequence}` via `extract_uniprot_id`, `{}` if `sequences is None`.
- **`microbiolink/utils/dmi_reader.py` (new)** — read-side helpers for **interpreting the DMI table
  produced by `workflow/dmi.py`** (the schema in `dmi.OUTPUT_COLUMNS`). `dmi.py` *writes* that table;
  these functions *read* a row of it back, inverting `dmi._to_output_row`'s `dmi_type` column routing.
  They live here — not in `dmi.py` — because only the consuming modules (Module 7, Module 8) call them;
  the producer never does. Generic to any module that consumes a `dmi.py` table:
  - `motif_side(row)` — `(motif_uniprot_id, start, end)`, routed by `dmi_type` (moved from
    `idr_filter._motif_side`). Both Module 7 and Module 8 use all three fields.
  - `motif_class_ids(row)` — the motif class id(s) on the motif-bearing side (`human_annotation` if
    `dmi_type == 'forward'`, else `bacterial_annotation`), split on `'|'` into a list (Module 6
    merges co-located classes with `'|'`). New; used by Module 8.
  - `validate_direction_sequences(table, human_sequences, bacterial_sequences)` — raise `ValueError`
    if the table has `'forward'` rows but no `human_sequences` (symmetric for `'reverse'`). The
    sequence-presence half of `idr_filter._validate_inputs`.

What stays in `idr_filter.py`, promoted to public so Module 8 can import them:

- `iupred_profile(sequence)`, `aiupred_profile(sequence, force_cpu, gpu_num)`,
  `cached_profile(method, sequence, force_cpu, gpu_num) -> (disorder, binding)` — the existing profile
  helpers, renamed from their `_`-prefixed forms to public names (behaviour unchanged; `cached_profile`
  keeps its module-level `lru_cache`, so a protein referenced by many rows — or by both Module 7 and
  Module 8 in one process — is profiled once).
- `validate_method(method)` — raise `ValueError` unless `method in {'iupred', 'aiupred'}`; the method
  half of `_validate_inputs`, split out as a public function (the sequence-presence half moves to
  `dmi_reader.validate_direction_sequences`).
- `disordered_mask(disorder_profile, cutoff) -> np.ndarray` (new) — boolean per-residue mask
  `disorder_profile >= cutoff`. Only Module 8 uses it, but it belongs with the disorder helpers.
- `disordered_region(mask, position) -> tuple[int, int] | None` (new) — the `(start, end)` half-open
  bounds of the maximal `True` run containing `position`, or `None` if `mask[position]` is `False`.

**`idr_filter.py` is then updated** to import `fasta.build_sequence_lookup`, `dmi_reader.motif_side`,
and `dmi_reader.validate_direction_sequences`, and to delete its now-duplicated `_build_sequence_lookup`,
`_motif_side`, and `_validate_inputs` (its two halves now come from the local `validate_method` +
`dmi_reader.validate_direction_sequences`); its `_iupred_profile`, `_aiupred_profile`, and
`_cached_profile` are promoted to the public `iupred_profile` / `aiupred_profile` / `cached_profile`,
and it gains the two new region helpers plus `validate_method`. Everything else in Module 7
(`OUTPUT_COLUMNS`, `_score_row`, `_to_output_row`, `_score_all_rows`, `filter_by_disorder`) keeps its
current behaviour — `filter_by_disorder` now calls `validate_method(method)` +
`dmi_reader.validate_direction_sequences(...)`, and `_score_row` calls the renamed `cached_profile`. The
refactor must be behaviour-preserving — re-run Module 7's AIUPred case-study regression as a check.

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

- **The null is pooled-disordered over-representation (option (a)).** Sample synthetic length-`D`
  sequences from the pooled disordered-region composition and re-match the motif regex; the regex is
  treated as a **black box** (`re.finditer`), so no pattern introspection is needed. See the decision
  document for why this beats whole-protein shuffle, single-region shuffle, and whole-pool shuffle.
- **Background is pooled per motif-bearing species, with leave-one-out.** `'forward'` (human) motifs
  are judged against the pooled composition of **human** disordered regions; `'reverse'` (bacterial)
  motifs against pooled **bacterial** disordered regions — the two proteomes differ compositionally.
  The target protein's own disordered residues are excluded from its background (leave-one-out); with
  a large pool this barely moves the frequencies, but it removes self-inclusion. *Caveat:* a small
  species pool (often the bacterial side) is a less stable background — the main reason the
  proteome-wide-background upgrade in the decision document matters most for the smaller species.
- **`D` is the disordered region length; observed counts are within that region.** The motif's region
  is the maximal disordered run containing its `start`; `D` is its length and the observed count is
  `re.finditer` matches of the class regex within that region's subsequence. Sampling is length-matched
  to `D`, so region length enters the null linearly, as intended.
- **Multiple-testing correction with Benjamini–Hochberg; `alpha` is the target FDR.** Raw MC p-values
  are corrected across all motif instances tested in the run, and a row passes iff its BH-adjusted
  **q-value** is `<= alpha` (the target false discovery rate, default `0.05`). This is the resolution
  of the decision document's open question 1; the full rationale (why FDR over FWER, BH over BY, the
  p-value-floor interaction) lives there. Implementation specifics:
  - **The unit of testing is the unique motif instance `(regex, protein, disordered region)`**, *not*
    the DMI row. Module 6's fan-out repeats one instance across many (domain × partner) rows carrying
    an identical p-value; BH runs over **de-duplicated** instances (each contributes one test to the
    count `m`), and each instance's q-value is broadcast back to every row that shares it. This also
    resolves the decision document's open question 2 (multiple motifs per protein → several
    independent instances). The instance key is `(regex, protein_id, region_start, region_end)`.
  - **The q-value is computed with `scipy.stats.false_discovery_control(pvals, method="bh")`** over the
    vector of unique-instance p-values (scipy is already a dependency; no hand-rolled BH). It does not
    depend on `alpha` — `alpha` only sets the pass/fail cutoff — so q-values are computed once and
    could be re-thresholded without recomputation.
  - **`iterations` must scale with the test count.** The MC p-value floor is `1/(iterations + 1)`,
    while BH's threshold for the most significant of `m` instances is `≈ alpha/m`; resolving that tail
    needs `iterations ≳ m/alpha`. The default `1000` suits small case-study runs; large runs should
    raise it (the per-unique-null cache keeps this affordable). Documented in the CLI help and docstring.
- **Two-phase flow (score all → BH → filter), so no per-row drop during scoring.** Because BH needs the
  full p-value vector before thresholding anything, Module 8 first scores every unique instance, then
  BH-corrects, then filters — unlike `filter_by_disorder`'s row-at-a-time score-and-drop. Rows that
  *cannot* be tested (see decision 4) are still dropped, but pass/fail is decided only after correction.
- **Return only passing rows, annotated.** Like `filter_by_disorder`, Module 8 returns passing rows
  only. The four new columns (`monte_carlo_hits`, `monte_carlo_pvalue`, `monte_carlo_qvalue`,
  `passes_monte_carlo`) are added, with `passes_monte_carlo` `True` on every returned row by
  construction.
- **Rows with several `'|'`-joined motif classes take the best (minimum-q) class.** Each class's regex
  is a separate instance tested independently and enters the BH pool on its own; the row's reported
  `monte_carlo_hits`/`monte_carlo_pvalue`/`monte_carlo_qvalue` are those of the class with the smallest
  q-value, and the row passes iff that minimum q-value `<= alpha`. Because BH is order-preserving, the
  minimum-q class is the minimum-p class, so this is a stable choice. Single-class rows — the common
  case — are unaffected. (The best-of-several selection is not itself separately corrected; accepted,
  as co-located classes are usually variants of one pattern, and every class tested is already counted
  in `m`. Alternative — require *all* listed classes to pass — rejected as over-strict.)

The remaining decisions follow the refactoring precedence rules and Module 7's precedent:

1. **One public function**, `filter_by_monte_carlo`, mirroring `filter_by_disorder`'s "take the
   upstream table shape directly, no dataclasses" style, plus the same disorder arguments Module 7
   takes (minus `binding_cutoff`, which option (a) does not use):
   ```python
   def filter_by_monte_carlo(
       interaction_table: pd.DataFrame,
       human_sequences: dict[str, str] | None,
       bacterial_sequences: dict[str, str] | None,
       method: str,
       disorder_cutoff: float,
       iterations: int = 1000,
       alpha: float = 0.05,
       seed: int = 0,
       force_cpu: bool = False,
       gpu_num: int = 0,
   ) -> pd.DataFrame
   ```
   `interaction_table` is Module 7's output. `human_sequences`/`bacterial_sequences` are the same
   header-keyed FASTA dicts every downstream module uses. No file I/O in the workflow module (file
   reading stays in `cli.py`).
2. **Output columns are the input's columns plus the four Monte Carlo columns** —
   `[*interaction_table.columns, 'monte_carlo_hits', 'monte_carlo_pvalue', 'monte_carlo_qvalue',
   'passes_monte_carlo']`. Whatever came in is passed through; four columns are appended. Passing rows
   are built as `row.to_dict()` updated with the four new keys, then assembled with
   `pd.DataFrame(records, columns=output_columns)` (explicit columns so an empty result keeps the
   schema, matching the other modules).
3. **Routing, lookup, and direction validation come from `utils`; disorder helpers come from
   `idr_filter`; regexes come from `dmi.load_motif_regexes`.** Module 8 imports `utils.fasta`,
   `utils.dmi_reader`, `workflow.dmi`, and `workflow.idr_filter` (for `cached_profile`,
   `disordered_mask`, `disordered_region`, `validate_method`).
4. **A row is dropped, not raised, when it cannot be tested**: its motif protein has no sequence in
   the corresponding lookup, its motif does not lie in a disordered region, or none of its motif class
   ids resolve to a known regex (mirrors Module 7's per-row skip). A wholesale-missing
   `human_sequences`/`bacterial_sequences` for a *present* direction is a hard `ValueError`
   (`dmi_reader.validate_direction_sequences`); an invalid `method` is a hard `ValueError`
   (`idr_filter.validate_method`).
5. **The sampling null is cached per `(regex, region-length, leave-one-out composition, iterations,
   seed)`, and each unique instance is scored once.** Module 6's cross-join means the same
   `(motif class, protein, region)` recurs across many partner rows; a `functools.lru_cache`'d sampler
   produces each unique combination's null match-count distribution only once, and the instance-level
   de-duplication (decision above) computes each instance's p-value once and enters it into the BH pool
   once. The cache key uses the *integer* leave-one-out disordered-residue count vector (exact, stable,
   and fully determining), so results are reproducible **across** runs, not just within one — a
   different input table that yields a different pool gets different keys.
   Reproducibility of the sampling itself comes from deriving the RNG deterministically from `seed`
   plus stable hashes of the regex, the composition, and `region_length`
   (`np.random.default_rng([seed, crc32(regex), crc32(composition), region_length])`).
6. **Requires the `idr` optional dependency; no new base dependency.** Option (a) computes disorder, so
   it needs the same IUPred/AIUPred install as Module 7. `numpy` (base) supplies the sampling; `re` and
   `zlib` are stdlib; BH correction uses `scipy.stats.false_discovery_control` (scipy is already a base
   dependency, `scipy>=1.13.1`). This is a deliberate change from the previous Module 8 plan, which
   needed no `idr` extra.

## Target implementation

### `microbiolink/utils/fasta.py` (extend)

```python
def build_sequence_lookup(sequences: dict[str, str] | None) -> dict[str, str]:
    """Reindex a header-keyed FASTA dict by UniProt accession. Empty dict if sequences is None."""
```

### `microbiolink/utils/dmi_reader.py` (new)

```python
"""Read-side helpers for interpreting a DMI interaction table produced by workflow/dmi.py.

workflow/dmi.py *builds* the table (columns per dmi.OUTPUT_COLUMNS, routed by dmi_type in
dmi._to_output_row); the functions here *read* a finished row of that table back — inverting the same
forward/reverse column routing — so consuming modules (idr_filter, monte_carlo) share one interpretation
of the schema. They live in utils, not in dmi.py, because only consumers call them; the producer never
reads its own output rows.
"""


def motif_side(row: pd.Series) -> tuple[str, int, int]:
    """(motif_uniprot_id, start, end): the human side if row['dmi_type'] == 'forward', else bacterial."""


def motif_class_ids(row: pd.Series) -> list[str]:
    """Motif class id(s) on the motif-bearing side, split on '|' (human_annotation if forward)."""


def validate_direction_sequences(
    table: pd.DataFrame,
    human_sequences: dict[str, str] | None,
    bacterial_sequences: dict[str, str] | None,
) -> None:
    """Raise ValueError if a present DMI direction lacks its sequences."""
```

### `microbiolink/workflow/dmi.py` (edit — add public accessor)

- Add `load_motif_regexes() -> dict[str, str]` (shown above). No other changes.

### `microbiolink/workflow/idr_filter.py` (edit — adopt the utils, publicise the disorder helpers)

- Delete `_build_sequence_lookup`, `_motif_side`, and `_validate_inputs`.
- Import `from ..utils import dmi_reader, fasta`; use `fasta.build_sequence_lookup`,
  `dmi_reader.motif_side`, `dmi_reader.validate_direction_sequences`.
- Rename `_iupred_profile` / `_aiupred_profile` / `_cached_profile` to public `iupred_profile` /
  `aiupred_profile` / `cached_profile` (behaviour unchanged; `cached_profile` keeps its `lru_cache`);
  update the internal call sites accordingly.
- Add the two new region helpers and the method validator (used by Module 8, colocated with the
  disorder logic):

  ```python
  def validate_method(method: str) -> None:
      """Raise ValueError unless method in {'iupred', 'aiupred'}."""


  def disordered_mask(disorder_profile: np.ndarray, cutoff: float) -> np.ndarray:
      """Boolean per-residue mask: disorder_profile >= cutoff."""


  def disordered_region(mask: np.ndarray, position: int) -> tuple[int, int] | None:
      """Half-open (start, end) of the maximal True run containing position, or None if mask[position] is False."""
  ```

- In `filter_by_disorder`, replace `_validate_inputs(...)` with `validate_method(method)` +
  `dmi_reader.validate_direction_sequences(dmi_table, human_sequences, bacterial_sequences)`; `_score_row`
  calls the renamed `cached_profile(...)`.
- `OUTPUT_COLUMNS`, `_score_row`, `_to_output_row`, `_score_all_rows`, `filter_by_disorder` keep their
  current behaviour.

### `microbiolink/workflow/monte_carlo.py` (new, core, no argparse)

```python
"""Monte Carlo over-representation filtering of DMIs against a pooled disordered-region composition."""

import functools
import re
import zlib

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

from ..utils import dmi_reader, fasta
from . import dmi, idr_filter

MONTE_CARLO_COLUMNS = [
    "monte_carlo_hits",
    "monte_carlo_pvalue",
    "monte_carlo_qvalue",
    "passes_monte_carlo",
]
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"  # canonical order for the composition count vector


def _count_matches(regex: str, sequence: str) -> int:
    """Number of re.finditer matches of regex in sequence (same matching Module 6 uses)."""


def _disordered_residue_counts(
    sequence: str, method: str, disorder_cutoff: float, force_cpu: bool, gpu_num: int
) -> np.ndarray:
    """Count vector over AMINO_ACIDS of the sequence's disordered residues.

    Profiles the sequence (idr_filter.cached_profile), masks residues at disorder_cutoff
    (idr_filter.disordered_mask), and tallies the masked residues into an AMINO_ACIDS-length vector.
    """


def _build_species_pools(
    interaction_table, human_lookup, bacterial_lookup, method, disorder_cutoff, force_cpu, gpu_num
) -> dict:
    """Precompute, per motif-bearing species, the pooled disordered-residue counts and per-protein
    counts (for leave-one-out).

    Iterates the unique motif proteins referenced by 'forward' rows (human) and 'reverse' rows
    (bacterial), accumulating a per-species global AMINO_ACIDS count vector and a
    {protein_id: count_vector} map. Uses idr_filter.cached_profile, so each protein is profiled once
    and the same cache serves per-row region lookup later.
    """


def _stable_hash(text: str) -> int:
    """Deterministic, run-stable integer hash of text (zlib.crc32), for RNG seeding."""


@functools.lru_cache(maxsize=None)
def _cached_null_counts(
    regex: str, region_length: int, loo_counts: tuple[int, ...], iterations: int, seed: int
) -> tuple[int, ...]:
    """Null match-count distribution for one (regex, region_length, leave-one-out composition).

    Draws `iterations` synthetic length-`region_length` sequences i.i.d. from the frequencies implied
    by loo_counts (over AMINO_ACIDS) and returns the per-sample regex match counts. The RNG is derived
    from seed and stable hashes of regex, loo_counts, and region_length, so results are reproducible
    and independent per combination; loo_counts is the integer count vector so the key fully
    determines the result across runs, and Module 6's fan-out reuses each combination's draw.
    """


def _row_instances(
    row, sequence_lookup, species_pool, regexes, method, disorder_cutoff,
    iterations, seed, force_cpu, gpu_num
) -> list[tuple[tuple, int, float]] | None:
    """Resolve one row into its testable motif instances and their raw MC p-values.

    Resolves the motif protein/position (dmi_reader.motif_side) and class ids
    (dmi_reader.motif_class_ids); locates the disordered region containing the motif start
    (idr_filter.cached_profile + disordered_mask + disordered_region); forms the leave-one-out count
    vector (species global counts minus this protein's counts). For each resolvable class regex,
    computes the observed count in the region subsequence and the raw p-value from _cached_null_counts
    (hits = null counts >= observed; pvalue = (hits + 1) / (iterations + 1)), returning one
    (instance_key, hits, pvalue) tuple per class where instance_key =
    (regex, protein_id, region_start, region_end). Returns None (row not testable → dropped) if the
    sequence is missing, the motif is not in a disordered region, or no class id resolves. No
    thresholding here — BH correction happens after every instance is scored.
    """


def _score_instances(
    interaction_table, human_lookup, bacterial_lookup, species_pools, regexes, method,
    disorder_cutoff, iterations, seed, force_cpu, gpu_num
) -> tuple[dict[tuple, tuple[int, float]], list[list[tuple] | None]]:
    """Phase 1: score every row's instances, de-duplicating to unique instances.

    Iterates the table (picking human vs bacterial lookup/pool by dmi_type), calling _row_instances.
    Returns (instance_stats, row_keys): instance_stats maps each unique instance_key -> (hits, pvalue)
    computed once (fan-out duplicates collapse); row_keys is aligned to the table, each element the
    list of instance_keys for that row (or None if the row was not testable). No BH yet.
    """


def _benjamini_hochberg(instance_stats: dict[tuple, tuple[int, float]]) -> dict[tuple, float]:
    """Phase 2: BH-adjust the unique-instance p-values into q-values.

    Runs scipy.stats.false_discovery_control(pvals, method="bh") over the instance p-values (one entry
    per unique instance, so m is the count of distinct tests) and returns instance_key -> qvalue.
    """


def _filter_rows(
    interaction_table, row_keys, instance_stats, instance_qvalues, alpha
) -> list[dict]:
    """Phase 3: broadcast q-values back to rows and keep the passing ones.

    For each testable row (row_keys entry is not None), looks up each class instance's (hits, pvalue)
    and qvalue, selects the class with the smallest qvalue (== smallest pvalue, BH being
    order-preserving), and — if that qvalue <= alpha — emits {**row.to_dict(), 'monte_carlo_hits': ...,
    'monte_carlo_pvalue': ..., 'monte_carlo_qvalue': ..., 'passes_monte_carlo': True}. Rows whose best
    q-value exceeds alpha, and rows with row_keys None, are dropped.
    """


def filter_by_monte_carlo(
    interaction_table: pd.DataFrame,
    human_sequences: dict[str, str] | None,
    bacterial_sequences: dict[str, str] | None,
    method: str,
    disorder_cutoff: float,
    iterations: int = 1000,
    alpha: float = 0.05,
    seed: int = 0,
    force_cpu: bool = False,
    gpu_num: int = 0,
) -> pd.DataFrame:
    """Filter Module 7's DMI table by a Monte Carlo over-representation test against pooled disorder.

    For each unique motif instance (regex, protein, disordered region), the class regex is counted in
    the region, then re-counted over `iterations` length-matched samples drawn from the pooled
    disordered-region composition of the motif's species (leave-one-out), giving a raw p-value
    (= (hits + 1) / (iterations + 1)). Raw p-values across all instances are Benjamini–Hochberg
    corrected, and a row passes if its best class's q-value is <= alpha (the target false discovery
    rate). Tests whether the motif is more frequent than the disordered proteome's composition predicts.

    Args:
        interaction_table: Module 7's IDR-filtered table (extra columns are passed through).
        human_sequences: Human FASTA header -> sequence mapping. Required for 'forward' rows.
        bacterial_sequences: Bacterial FASTA header -> sequence mapping. Required for 'reverse'.
        method: 'iupred' or 'aiupred' (disorder predictor).
        disorder_cutoff: Per-residue disorder score at/above which a residue is disordered.
        iterations: Number of samples per motif instance. The MC p-value floor is 1/(iterations + 1);
            for a run testing m instances at FDR alpha, resolving the BH tail needs iterations >~
            m/alpha, so raise this for large runs (the per-unique-null cache keeps it affordable).
        alpha: Target false discovery rate. A row passes iff its BH-adjusted q-value is <= alpha.
        seed: Random seed for reproducible sampling.
        force_cpu: Force CPU inference (aiupred only; ignored for iupred).
        gpu_num: GPU index (aiupred only; ignored for iupred).

    Returns:
        interaction_table restricted to passing rows, with monte_carlo_hits, monte_carlo_pvalue,
        monte_carlo_qvalue, and passes_monte_carlo columns appended.

    Raises:
        ValueError: If method is invalid, or a present direction's sequences are missing.
    """
```

- `filter_by_monte_carlo` (orchestrates the three phases): `idr_filter.validate_method(method)`;
  `dmi_reader.validate_direction_sequences(interaction_table, human_sequences, bacterial_sequences)`;
  `human_lookup = fasta.build_sequence_lookup(human_sequences)` (and bacterial);
  `regexes = dmi.load_motif_regexes()`; `species_pools = _build_species_pools(...)`;
  `instance_stats, row_keys = _score_instances(...)`;
  `instance_qvalues = _benjamini_hochberg(instance_stats)`;
  `records = _filter_rows(interaction_table, row_keys, instance_stats, instance_qvalues, alpha)`;
  `output_columns = [*interaction_table.columns, *MONTE_CARLO_COLUMNS]`;
  `return pd.DataFrame(records, columns=output_columns)`.
- `_row_instances`: `protein_id, start, _ = dmi_reader.motif_side(row)`;
  `sequence = sequence_lookup.get(protein_id)` (`None` → return `None`);
  `disorder_profile, _ = idr_filter.cached_profile(method, sequence, force_cpu, gpu_num)`;
  `region = idr_filter.disordered_region(idr_filter.disordered_mask(disorder_profile, disorder_cutoff),
  start)` (`None` → return `None`); `region_seq = sequence[region_start:region_end]`,
  `D = region_end - region_start`; `loo_counts = tuple(species_global - species_per_protein[protein_id])`;
  `resolved = [(mid, regexes[mid]) for mid in dmi_reader.motif_class_ids(row) if mid in regexes]`
  (empty → return `None`); for each `regex`, `observed = _count_matches(regex, region_seq)`,
  `null = _cached_null_counts(regex, D, loo_counts, iterations, seed)`,
  `hits = sum(count >= observed for count in null)`, `pvalue = (hits + 1) / (iterations + 1)`,
  `key = (regex, protein_id, region_start, region_end)`; return `[(key, hits, pvalue), ...]`.
- `_score_instances`: `for _, row in interaction_table.iterrows()`, pick human vs bacterial lookup/pool
  by `row['dmi_type'] == 'forward'`, call `_row_instances`; append its result (list or `None`) to
  `row_keys`; for a non-`None` result, record each `key -> (hits, pvalue)` in `instance_stats` (setdefault,
  so a fan-out duplicate is stored once). Return `(instance_stats, row_keys)`.
- `_benjamini_hochberg`: `keys = list(instance_stats)`;
  `pvals = [instance_stats[k][1] for k in keys]`;
  `qvals = false_discovery_control(pvals, method="bh")`; `return dict(zip(keys, qvals))`.
- `_filter_rows`: `for row, keys in zip(rows, row_keys)`, skip `None`; among `keys`, pick the one with
  the smallest `instance_qvalues[key]`; `hits, pvalue = instance_stats[key]`,
  `qvalue = instance_qvalues[key]`; if `qvalue <= alpha`, append `{**row.to_dict(),
  'monte_carlo_hits': hits, 'monte_carlo_pvalue': pvalue, 'monte_carlo_qvalue': qvalue,
  'passes_monte_carlo': True}`.

### `microbiolink/cli.py` (extend, argparse only)

- `_add_monte_carlo_source_arguments(parser)`: `--interaction_file` (required — Module 7 output CSV),
  `--human_fasta_file` (optional), `--bacterial_fasta_file` (optional).
- `_add_monte_carlo_disorder_arguments(parser)`: `--method` (`iupred`/`aiupred`), `--disorder_cutoff`
  (`float`), `--force_cpu` (flag), `--gpu_num` (`int`, `0`) — the same disorder arguments Module 7's
  CLI exposes, minus `--binding_cutoff`.
- `_add_monte_carlo_test_arguments(parser)`: `--iterations` (`int`, `1000`; help notes the p-value
  floor `1/(iterations + 1)` and the `iterations >~ m/alpha` rule of thumb for large runs), `--alpha`
  (`float`, `0.05`; help: "target false discovery rate; a motif passes if its BH-corrected q-value <=
  alpha"), `--seed` (`int`, `0`).
- `_build_monte_carlo_parser()`: the three helpers plus `-o/--output_file` (required).
- `monte_carlo_filter()` entry point: `pd.read_csv(args.interaction_file)`, read the two optional FASTA
  files via `fasta.read_fasta_sequences`, call `monte_carlo.filter_by_monte_carlo(...)`,
  `result.to_csv(args.output_file, index=False)`, `return 0`.

### `pyproject.toml`

- Add `microbiolink-monte-carlo = "microbiolink.cli:monte_carlo_filter"` under `[project.scripts]`.
- No new third-party dependency, but Module 8 now **requires the `idr` optional dependency** (the
  IUPred/AIUPred install) to compute disorder, exactly like Module 7. Note this in the module's
  install requirements alongside Module 7.

## Migration checklist

### 1. Shared-code extraction (Module 7 refactor)
- [x] Add `build_sequence_lookup` to `microbiolink/utils/fasta.py`.
- [x] Create `microbiolink/utils/dmi_reader.py` with `motif_side`, `motif_class_ids`,
      `validate_direction_sequences`.
- [x] Update `idr_filter.py`: import `fasta.build_sequence_lookup`, `dmi_reader.motif_side`,
      `dmi_reader.validate_direction_sequences`; delete `_build_sequence_lookup`, `_motif_side`,
      `_validate_inputs`; rename `_iupred_profile` / `_aiupred_profile` / `_cached_profile` to public
      `iupred_profile` / `aiupred_profile` / `cached_profile`; add `validate_method`, `disordered_mask`,
      and `disordered_region`.
- [ ] Re-run the Module 7 AIUPred case-study regression from the Module 7 plan to confirm the refactor
      is behaviour-preserving (same passing-row set / scores as before). *(No pre-refactor baseline
      fixture exists to diff against, so a strict same-output regression could not be performed. Two
      pieces of evidence instead: a synthetic behaviour-preserving check on `filter_by_disorder`
      (pass/fail set, appended schema, both `ValueError` paths) passed; and the refactored Module 7 ran
      cleanly on the real AIUPred case study (615 proteins → 7,904 disordered+binding interactions),
      producing a sensible result — see the Module 8 end-to-end note in §6.)*

### 2. Resource accessor + packaging
- [x] Add `load_motif_regexes()` to `microbiolink/workflow/dmi.py`.
- [x] Add `microbiolink-monte-carlo = "microbiolink.cli:monte_carlo_filter"` to `[project.scripts]` and
      note the `idr` extra requirement. *(Script added; the shared `idr` optional dependency already
      exists in `pyproject.toml` and now also covers Module 8.)*

### 3. Core module
- [x] Create `microbiolink/workflow/monte_carlo.py` with `_count_matches`,
      `_disordered_residue_counts`, `_build_species_pools`, `_stable_hash`, `_cached_null_counts`,
      `_row_instances`, `_score_instances`, `_benjamini_hochberg`, `_filter_rows`,
      `filter_by_monte_carlo`, importing from `utils.fasta`, `utils.dmi_reader`, `workflow.dmi`,
      `workflow.idr_filter` (disorder helpers), and `scipy.stats.false_discovery_control`.
      *(To satisfy the repo's 50-line-per-function hook the private helpers were re-partitioned: a
      `_ScoringParams` NamedTuple bundles the shared args; `_score_instances` became
      `_score_all_instances` (folding in the pool build); and `_locate_region`, `_score_class`, and a
      thin `_run_filter` wrapper were added. Public API and behaviour are unchanged.)*
- [x] Confirm the raw p-value is `(hits + 1) / (iterations + 1)` with a hit being a sample whose match
      count is `>=` observed, sampling is i.i.d. from the leave-one-out pooled composition, and length
      matches the motif's disordered region `D`.
- [x] Confirm background pooling is **per motif-bearing species** with leave-one-out.
- [x] Confirm the two-phase flow: unique instances `(regex, protein_id, region bounds)` are scored and
      de-duplicated, BH (`false_discovery_control(..., method="bh")`) runs once over the unique-instance
      p-vector, and q-values are broadcast back so fan-out duplicates share one q-value. Pass iff
      q-value `<= alpha` (target FDR); multi-class rows take the minimum-q class.
- [x] Confirm a missing-sequence / motif-not-in-disorder / unresolvable-motif row is dropped, while a
      wholesale-missing sequences dict for a present direction and an invalid `method` both raise
      `ValueError`.
- [x] Confirm output columns = input columns + the four Monte Carlo columns (`monte_carlo_hits`,
      `monte_carlo_pvalue`, `monte_carlo_qvalue`, `passes_monte_carlo`).

### 4. CLI wiring
- [x] Add `_add_monte_carlo_source_arguments`, `_add_monte_carlo_disorder_arguments`,
      `_add_monte_carlo_test_arguments`, `_build_monte_carlo_parser`, `monte_carlo_filter()` to
      `cli.py`.

### 5. Delete old code (per Q11)
- [x] None. No Monte Carlo file exists in the `refactoring` working tree (the `development` script
      implements a different null and is not ported).

### 6. Verification — no committed fixture (manual sign-off, as Modules 5/6 and Module 7's IUPred path)
- [x] Confirm there is no `case_study_output` Monte Carlo fixture (checked: none exists).
- [x] Deterministic seeded unit check on synthetic sequences with a synthetic pool: a motif repeated
      far more than the pooled composition predicts passes (low p-value); a motif occurring about as
      often as chance fails (high p-value); confirm same-`seed` runs are identical and a different
      `seed` perturbs `monte_carlo_hits` but not the pass/fail conclusion.
- [x] Unit-check the BH step: on a hand-built p-vector, `monte_carlo_qvalue` matches
      `scipy.stats.false_discovery_control(pvals, method="bh")`, q-values are monotone non-decreasing in
      p rank and `>=` the raw p-values, fan-out duplicate rows of one instance share a single q-value,
      and lowering `alpha` only shrinks (never grows) the passing set. Sanity: an all-null synthetic
      table (no true motifs) passes roughly `alpha`-fraction or fewer instances.
- [x] Unit-check the disorder-region helpers: `disordered_mask` thresholds correctly and
      `disordered_region` returns the run containing a position (and `None` for an ordered position),
      including region-edge cases.
- [x] Run end-to-end on a real Module 7 output (the AIUPred-filtered table) with
      `microbiolink-monte-carlo`: confirm the appended-4-column schema, output rows ⊆ input rows,
      `passes_monte_carlo` all-`True`, plausible pass counts. *(Ran the real case study: 615 human
      proteins → forward DMI (92,099) → Module 7 AIUPred IDR filter, cutoffs 0.5/0.5 (7,904) → Module 8
      (iterations=1000, alpha=0.05, seed=0) = **1,011 passing rows / 558 unique instances**. Schema,
      subset, and all-`True` all confirmed. At iterations=20000 → 1,161 rows, illustrating the p-value
      floor's effect on BH as expected.)*
- [x] Verify `dmi_type` routing and per-species pooling directly with a synthetic 2-row table (one
      `forward`, one `reverse`) by swapping motif-rich vs. motif-poor synthetic sequences between the
      human and bacterial dicts.
- [x] Shown to the user for confirmation.

## Verification

- `uv run ruff format` / `uv run ruff check` scoped to the new files (`utils/dmi_reader.py`,
  `workflow/monte_carlo.py`); `fasta.py`, `dmi.py`, `idr_filter.py`, `cli.py`, and
  `pyproject.toml` edits reviewed by eye against surrounding style (modified, not new).
- `uv run ty check`
- The Module 7 regression re-run (step 1) plus the seeded/synthetic and disorder-helper checks above.
