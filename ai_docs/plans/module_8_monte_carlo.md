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
   `>=` the observed count. A motif **passes** if its p-value is `<= alpha`.

This is a composition-matched *sampling* null, not a permutation shuffle: sampling from the pooled
frequency vector is the efficient canonical form of the identical composition-only null (see the
decision document's rejected options 3–4).

**Disorder dependency is required.** Unlike the earlier whole-protein-shuffle design, option (a)
**needs a disorder profile** — both to locate each motif's disordered region and to build the pooled
background. Module 8 computes the profile itself (Module 7 emits per-window *scores*, not per-residue
profiles or region boundaries), reusing the shared disorder helpers, and therefore takes the same
`method` / `disorder_cutoff` / device arguments as Module 7 and requires the same optional `idr`
dependency. This reverses the previous "Module 8 needs no disorder install" property; it is
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

## Shared code moves to `utils` (prerequisite refactor of Module 7)

Module 7 (`idr_filter.py`) owns helpers that Module 8 now needs identically. Rather than have Module 8
reach into `idr_filter`'s privates, the shared pieces move into `utils` and both modules import them
from there (the Q5/Q6 "shared internal utility, single implementation" rule). **This is a change from
the previous Module 8 plan**, where the disorder-profile helpers stayed in `idr_filter` because the
old whole-protein null had no use for them; option (a) uses them, so they move.

What moves, and where:

- **`microbiolink/utils/fasta.py` (extend)** — add `build_sequence_lookup(sequences)` (moved from
  `idr_filter._build_sequence_lookup`): reindex a header-keyed FASTA dict to
  `{uniprot_accession: sequence}` via `extract_uniprot_id`, `{}` if `sequences is None`.
- **`microbiolink/utils/dmi_table.py` (new)** — helpers that operate on a Module 6/7 interaction row,
  generic to any consuming module:
  - `motif_side(row)` — `(motif_uniprot_id, start, end)`, routed by `dmi_type` (moved from
    `idr_filter._motif_side`). Both Module 7 and Module 8 use all three fields.
  - `motif_class_ids(row)` — the motif class id(s) on the motif-bearing side (`human_annotation` if
    `dmi_type == 'forward'`, else `bacterial_annotation`), split on `'|'` into a list (Module 6
    merges co-located classes with `'|'`). New; used by Module 8.
  - `validate_direction_sequences(table, human_sequences, bacterial_sequences)` — raise `ValueError`
    if the table has `'forward'` rows but no `human_sequences` (symmetric for `'reverse'`). The
    sequence-presence half of `idr_filter._validate_inputs`.
- **`microbiolink/utils/disorder.py` (new)** — the disorder-profile helpers, moved out of
  `idr_filter.py` and made public so both modules share one implementation:
  - `iupred_profile(sequence)`, `aiupred_profile(sequence, force_cpu, gpu_num)` — moved verbatim from
    `idr_filter._iupred_profile` / `_aiupred_profile`.
  - `cached_profile(method, sequence, force_cpu, gpu_num) -> (disorder, binding)` — moved from
    `idr_filter._cached_profile`, keeping its module-level `lru_cache` so a protein referenced by many
    rows (or by both Module 7 and Module 8) is profiled once.
  - `validate_method(method)` — raise `ValueError` unless `method in {'iupred', 'aiupred'}` (the
    method half of `idr_filter._validate_inputs`), now shared since both modules take `method`.
  - `disordered_mask(disorder_profile, cutoff) -> np.ndarray` (new) — boolean per-residue mask
    `disorder_profile >= cutoff`. Only Module 8 uses it, but it belongs with the disorder helpers.
  - `disordered_region(mask, position) -> tuple[int, int] | None` (new) — the `(start, end)` half-open
    bounds of the maximal `True` run containing `position`, or `None` if `mask[position]` is `False`.

**`idr_filter.py` is then updated** to import `fasta.build_sequence_lookup`, `dmi_table.motif_side` /
`dmi_table.validate_direction_sequences`, and `disorder.cached_profile` / `disorder.validate_method`,
and to delete its now-duplicated `_build_sequence_lookup`, `_motif_side`, `_iupred_profile`,
`_aiupred_profile`, `_cached_profile`, and `_validate_inputs` (its two halves now come from
`disorder.validate_method` + `dmi_table.validate_direction_sequences`). Everything else in Module 7
(`OUTPUT_COLUMNS`, `_score_row`, `_to_output_row`, `_score_all_rows`, `filter_by_disorder`) keeps its
current behaviour — it now calls `disorder.cached_profile` instead of the local one. The refactor must
be behaviour-preserving — re-run Module 7's AIUPred case-study regression as a check.

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
- **Drop failing rows (do not annotate-and-keep).** Like `filter_by_disorder`, Module 8 returns only
  passing rows. The three new columns (`monte_carlo_hits`, `monte_carlo_pvalue`, `passes_monte_carlo`)
  are still added, with `passes_monte_carlo` `True` on every returned row by construction.
- **Rows with several `'|'`-joined motif classes take the minimum p-value.** Each class's regex is
  tested independently against the same region and the same background; the row's reported
  `monte_carlo_hits`/`monte_carlo_pvalue` are those of the class with the smallest p-value, and the row
  passes iff that minimum p-value `<= alpha` (over-represented under at least its best-supported
  class). Single-class rows — the common case — are unaffected. (Alternative considered: require *all*
  listed classes to pass; rejected as over-strict since co-located classes are usually variants of one
  pattern.)

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
2. **Output columns are the input's columns plus the three Monte Carlo columns** —
   `[*interaction_table.columns, 'monte_carlo_hits', 'monte_carlo_pvalue', 'passes_monte_carlo']`.
   Whatever came in is passed through; three columns are appended. Passing rows are built as
   `row.to_dict()` updated with the three new keys, then assembled with
   `pd.DataFrame(records, columns=output_columns)` (explicit columns so an empty result keeps the
   schema, matching the other modules).
3. **Routing, lookup, disorder, and validation come from `utils`**; regexes come from
   `dmi.load_motif_regexes`. Module 8 imports `utils.fasta`, `utils.dmi_table`, `utils.disorder`, and
   `workflow.dmi` — never `idr_filter`.
4. **A row is dropped, not raised, when it cannot be tested**: its motif protein has no sequence in
   the corresponding lookup, its motif does not lie in a disordered region, or none of its motif class
   ids resolve to a known regex (mirrors Module 7's per-row skip). A wholesale-missing
   `human_sequences`/`bacterial_sequences` for a *present* direction is a hard `ValueError`
   (`dmi_table.validate_direction_sequences`); an invalid `method` is a hard `ValueError`
   (`disorder.validate_method`).
5. **The sampling null is cached per `(regex, region-length, leave-one-out composition, iterations,
   seed)`.** Module 6's cross-join means the same `(motif class, protein, region)` recurs across many
   partner rows; a `functools.lru_cache`'d sampler produces each unique combination's null match-count
   distribution only once. The cache key uses the *integer* leave-one-out disordered-residue count
   vector (exact, stable, and fully determining), so results are reproducible **across** runs, not
   just within one — a different input table that yields a different pool gets different keys.
   Reproducibility of the sampling itself comes from deriving the RNG deterministically from `seed`
   plus stable hashes of the regex, the composition, and `region_length`
   (`np.random.default_rng([seed, crc32(regex), crc32(composition), region_length])`).
6. **Requires the `idr` optional dependency.** Option (a) computes disorder, so it needs the same
   IUPred/AIUPred install as Module 7. `numpy` (base) supplies the sampling; `re` is stdlib. This is a
   deliberate change from the previous Module 8 plan, which needed no `idr` extra.

## Target implementation

### `microbiolink/utils/fasta.py` (extend)

```python
def build_sequence_lookup(sequences: dict[str, str] | None) -> dict[str, str]:
    """Reindex a header-keyed FASTA dict by UniProt accession. Empty dict if sequences is None."""
```

### `microbiolink/utils/dmi_table.py` (new)

```python
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

### `microbiolink/utils/disorder.py` (new — moved profile helpers + new region helpers)

```python
def iupred_profile(sequence: str) -> tuple[np.ndarray, np.ndarray]:
    """Per-residue IUPred2 disorder and ANCHOR2 binding scores. (moved from idr_filter)"""


def aiupred_profile(sequence: str, force_cpu: bool, gpu_num: int) -> tuple[np.ndarray, np.ndarray]:
    """Per-residue AIUPred disorder and binding scores. (moved from idr_filter)"""


@functools.lru_cache(maxsize=None)
def cached_profile(method: str, sequence: str, force_cpu: bool, gpu_num: int) -> tuple[np.ndarray, np.ndarray]:
    """Cached (disorder, binding) profile, keyed by method/sequence/device. (moved from idr_filter)"""


def validate_method(method: str) -> None:
    """Raise ValueError unless method in {'iupred', 'aiupred'}."""


def disordered_mask(disorder_profile: np.ndarray, cutoff: float) -> np.ndarray:
    """Boolean per-residue mask: disorder_profile >= cutoff."""


def disordered_region(mask: np.ndarray, position: int) -> tuple[int, int] | None:
    """Half-open (start, end) of the maximal True run containing position, or None if mask[position] is False."""
```

### `microbiolink/workflow/dmi.py` (edit — add public accessor)

- Add `load_motif_regexes() -> dict[str, str]` (shown above). No other changes.

### `microbiolink/workflow/idr_filter.py` (edit — adopt the utils, delete duplicates)

- Delete `_build_sequence_lookup`, `_motif_side`, `_iupred_profile`, `_aiupred_profile`,
  `_cached_profile`, and `_validate_inputs`.
- Import `from ..utils import dmi_table, disorder, fasta`; use `fasta.build_sequence_lookup`,
  `dmi_table.motif_side`, `dmi_table.validate_direction_sequences`, `disorder.cached_profile`,
  `disorder.validate_method`.
- In `filter_by_disorder`, replace `_validate_inputs(...)` with `disorder.validate_method(method)` +
  `dmi_table.validate_direction_sequences(dmi_table, human_sequences, bacterial_sequences)`; `_score_row`
  calls `disorder.cached_profile(...)`.
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

from ..utils import dmi_table, disorder, fasta
from . import dmi

MONTE_CARLO_COLUMNS = ["monte_carlo_hits", "monte_carlo_pvalue", "passes_monte_carlo"]
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"  # canonical order for the composition count vector


def _count_matches(regex: str, sequence: str) -> int:
    """Number of re.finditer matches of regex in sequence (same matching Module 6 uses)."""


def _disordered_residue_counts(
    sequence: str, method: str, disorder_cutoff: float, force_cpu: bool, gpu_num: int
) -> np.ndarray:
    """Count vector over AMINO_ACIDS of the sequence's disordered residues.

    Profiles the sequence (disorder.cached_profile), masks residues at disorder_cutoff
    (disorder.disordered_mask), and tallies the masked residues into an AMINO_ACIDS-length vector.
    """


def _build_species_pools(
    interaction_table, human_lookup, bacterial_lookup, method, disorder_cutoff, force_cpu, gpu_num
) -> dict:
    """Precompute, per motif-bearing species, the pooled disordered-residue counts and per-protein
    counts (for leave-one-out).

    Iterates the unique motif proteins referenced by 'forward' rows (human) and 'reverse' rows
    (bacterial), accumulating a per-species global AMINO_ACIDS count vector and a
    {protein_id: count_vector} map. Uses disorder.cached_profile, so each protein is profiled once
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


def _score_row(
    row, sequence_lookup, species_pool, regexes, method, disorder_cutoff,
    iterations, alpha, seed, force_cpu, gpu_num
) -> tuple[bool, int, float] | None:
    """Run the over-representation test for one row's motif class(es).

    Resolves the motif protein/position (dmi_table.motif_side) and class ids
    (dmi_table.motif_class_ids); locates the disordered region containing the motif start
    (disorder.cached_profile + disordered_mask + disordered_region); forms the leave-one-out count
    vector (species global counts minus this protein's counts); for each resolvable class regex,
    computes the observed count in the region and the p-value from _cached_null_counts; keeps the
    smallest-p-value class. Returns (passes, monte_carlo_hits, monte_carlo_pvalue), or None (row
    dropped) if the sequence is missing, the motif is not in a disordered region, or no class id
    resolves. passes is (pvalue <= alpha); hits = number of null counts >= observed;
    pvalue = (hits + 1) / (iterations + 1).
    """


def _filter_all_rows(
    interaction_table, human_lookup, bacterial_lookup, species_pools, regexes, method,
    disorder_cutoff, iterations, alpha, seed, force_cpu, gpu_num
) -> list[dict]:
    """Run the test on every row; collect passing rows as row.to_dict() records with the three
    Monte Carlo keys added. Picks human_lookup/pool for 'forward' rows, bacterial otherwise."""


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

    For each row, the motif class's regex is counted in the motif's disordered region, then re-counted
    over `iterations` length-matched samples drawn from the pooled disordered-region composition of the
    motif's species (leave-one-out); the motif passes if its p-value (= (hits + 1) / (iterations + 1))
    is <= alpha. Tests whether the motif is more frequent than the disordered proteome's composition
    predicts.

    Args:
        interaction_table: Module 7's IDR-filtered table (extra columns are passed through).
        human_sequences: Human FASTA header -> sequence mapping. Required for 'forward' rows.
        bacterial_sequences: Bacterial FASTA header -> sequence mapping. Required for 'reverse'.
        method: 'iupred' or 'aiupred' (disorder predictor).
        disorder_cutoff: Per-residue disorder score at/above which a residue is disordered.
        iterations: Number of samples per motif.
        alpha: Maximum p-value a motif may have and still pass.
        seed: Random seed for reproducible sampling.
        force_cpu: Force CPU inference (aiupred only; ignored for iupred).
        gpu_num: GPU index (aiupred only; ignored for iupred).

    Returns:
        interaction_table restricted to passing rows, with monte_carlo_hits, monte_carlo_pvalue, and
        passes_monte_carlo columns appended.

    Raises:
        ValueError: If method is invalid, or a present direction's sequences are missing.
    """
```

- `filter_by_monte_carlo`: `disorder.validate_method(method)`;
  `dmi_table.validate_direction_sequences(interaction_table, human_sequences, bacterial_sequences)`;
  `human_lookup = fasta.build_sequence_lookup(human_sequences)` (and bacterial);
  `regexes = dmi.load_motif_regexes()`;
  `species_pools = _build_species_pools(...)`; `records = _filter_all_rows(...)`;
  `output_columns = [*interaction_table.columns, *MONTE_CARLO_COLUMNS]`;
  `return pd.DataFrame(records, columns=output_columns)`.
- `_score_row`: `protein_id, start, _ = dmi_table.motif_side(row)`;
  `sequence = sequence_lookup.get(protein_id)` (`None` → drop);
  `disorder_profile, _ = disorder.cached_profile(method, sequence, force_cpu, gpu_num)`;
  `region = disorder.disordered_region(disorder.disordered_mask(disorder_profile, disorder_cutoff),
  start)` (`None` → drop); `region_seq = sequence[region_start:region_end]`, `D = region_end -
  region_start`; `loo_counts = tuple(species_global - species_per_protein[protein_id])`;
  `resolved = [(mid, regexes[mid]) for mid in dmi_table.motif_class_ids(row) if mid in regexes]`
  (empty → drop); for each `regex`, `observed = _count_matches(regex, region_seq)`,
  `null = _cached_null_counts(regex, D, loo_counts, iterations, seed)`,
  `hits = sum(count >= observed for count in null)`, `pvalue = (hits + 1) / (iterations + 1)`; keep the
  class with the smallest `pvalue`; `passes = pvalue <= alpha`.
- `_filter_all_rows`: `for _, row in interaction_table.iterrows()`, pick human vs bacterial
  lookup/pool by `row['dmi_type'] == 'forward'`, call `_score_row`; skip `None`; when `passes`, append
  `{**row.to_dict(), 'monte_carlo_hits': hits, 'monte_carlo_pvalue': pvalue,
  'passes_monte_carlo': True}`.

### `microbiolink/cli.py` (extend, argparse only)

- `_add_monte_carlo_source_arguments(parser)`: `--interaction_file` (required — Module 7 output CSV),
  `--human_fasta_file` (optional), `--bacterial_fasta_file` (optional).
- `_add_monte_carlo_disorder_arguments(parser)`: `--method` (`iupred`/`aiupred`), `--disorder_cutoff`
  (`float`), `--force_cpu` (flag), `--gpu_num` (`int`, `0`) — the same disorder arguments Module 7's
  CLI exposes, minus `--binding_cutoff`.
- `_add_monte_carlo_test_arguments(parser)`: `--iterations` (`int`, `1000`), `--alpha` (`float`,
  `0.05`), `--seed` (`int`, `0`).
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
- [ ] Add `build_sequence_lookup` to `microbiolink/utils/fasta.py`.
- [ ] Create `microbiolink/utils/dmi_table.py` with `motif_side`, `motif_class_ids`,
      `validate_direction_sequences`.
- [ ] Create `microbiolink/utils/disorder.py` with `iupred_profile`, `aiupred_profile`,
      `cached_profile`, `validate_method` (moved from `idr_filter`) plus new `disordered_mask` and
      `disordered_region`.
- [ ] Update `idr_filter.py` to import from utils and delete its now-duplicated helpers
      (`_build_sequence_lookup`, `_motif_side`, `_iupred_profile`, `_aiupred_profile`,
      `_cached_profile`, `_validate_inputs`).
- [ ] Re-run the Module 7 AIUPred case-study regression from the Module 7 plan to confirm the refactor
      is behaviour-preserving (same passing-row set / scores as before).

### 2. Resource accessor + packaging
- [ ] Add `load_motif_regexes()` to `microbiolink/workflow/dmi.py`.
- [ ] Add `microbiolink-monte-carlo = "microbiolink.cli:monte_carlo_filter"` to `[project.scripts]` and
      note the `idr` extra requirement.

### 3. Core module
- [ ] Create `microbiolink/workflow/monte_carlo.py` with `_count_matches`,
      `_disordered_residue_counts`, `_build_species_pools`, `_stable_hash`, `_cached_null_counts`,
      `_score_row`, `_filter_all_rows`, `filter_by_monte_carlo`, importing only from `utils` and
      `workflow.dmi`.
- [ ] Confirm the p-value is `(hits + 1) / (iterations + 1)` with a hit being a sample whose match
      count is `>=` observed, sampling is i.i.d. from the leave-one-out pooled composition, and length
      matches the motif's disordered region `D`.
- [ ] Confirm background pooling is **per motif-bearing species** with leave-one-out.
- [ ] Confirm a missing-sequence / motif-not-in-disorder / unresolvable-motif row is dropped, while a
      wholesale-missing sequences dict for a present direction and an invalid `method` both raise
      `ValueError`.
- [ ] Confirm output columns = input columns + the three Monte Carlo columns.

### 4. CLI wiring
- [ ] Add `_add_monte_carlo_source_arguments`, `_add_monte_carlo_disorder_arguments`,
      `_add_monte_carlo_test_arguments`, `_build_monte_carlo_parser`, `monte_carlo_filter()` to
      `cli.py`.

### 5. Delete old code (per Q11)
- [ ] None. No Monte Carlo file exists in the `refactoring` working tree (the `development` script
      implements a different null and is not ported).

### 6. Verification — no committed fixture (manual sign-off, as Modules 5/6 and Module 7's IUPred path)
- [ ] Confirm there is no `case_study_output` Monte Carlo fixture (checked: none exists).
- [ ] Deterministic seeded unit check on synthetic sequences with a synthetic pool: a motif repeated
      far more than the pooled composition predicts passes (low p-value); a motif occurring about as
      often as chance fails (high p-value); confirm same-`seed` runs are identical and a different
      `seed` perturbs `monte_carlo_hits` but not the pass/fail conclusion.
- [ ] Unit-check the disorder-region helpers: `disordered_mask` thresholds correctly and
      `disordered_region` returns the run containing a position (and `None` for an ordered position),
      including region-edge cases.
- [ ] Run end-to-end on a real Module 7 output (the AIUPred-filtered table) with
      `microbiolink-monte-carlo`: confirm the appended-3-column schema, output rows ⊆ input rows,
      `passes_monte_carlo` all-`True`, plausible pass counts.
- [ ] Verify `dmi_type` routing and per-species pooling directly with a synthetic 2-row table (one
      `forward`, one `reverse`) by swapping motif-rich vs. motif-poor synthetic sequences between the
      human and bacterial dicts.
- [ ] Shown to the user for confirmation.

## Verification

- `uv run ruff format` / `uv run ruff check` scoped to the new files (`utils/dmi_table.py`,
  `utils/disorder.py`, `workflow/monte_carlo.py`); `fasta.py`, `dmi.py`, `idr_filter.py`, `cli.py`, and
  `pyproject.toml` edits reviewed by eye against surrounding style (modified, not new).
- `uv run ty check`
- The Module 7 regression re-run (step 1) plus the seeded/synthetic and disorder-helper checks above.
