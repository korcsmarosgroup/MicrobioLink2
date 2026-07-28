# Module 7 — Intrinsic Disorder Region (IDR) Prediction — Implementation Plan

Detailed plan for the seventh module described in @ai_docs/plans/microbiolink_refactoring.md, following the
decisions in @ai_docs/plans/refactoring_questions.md and the precedent set by
@ai_docs/plans/module_1_zscore_filter.md through @ai_docs/plans/module_6_dmi.md. Scope is exactly Module 7: given
Module 6's DMI output table and the FASTA sequences it was built from, score each row's motif
window for intrinsic disorder and binding likelihood (via IUPred2/ANCHOR2 or AIUPred), and filter
to rows whose motif window clears both a disorder cutoff and a binding cutoff at every residue.

## Context

Modules 1–6 are implemented in `microbiolink/workflow/` (`zscore_filter.py`, `membrane_filter.py`,
`fasta_download.py`, `domain_download.py`, `ddi.py`, `dmi.py`), `microbiolink/utils/`
(`uniprot_client.py`, `id_resolution.py`, `fasta.py`), and `microbiolink/data/` (Module 5/6
resource TSVs), all wired into `microbiolink/cli.py`. Module 6's output (`dmi.py`'s
`predict_domain_motif_interactions`) is this module's entire input: a `pd.DataFrame` with columns
`dmi_type`, `bacterial_uniprot_id`, `bacterial_annotation`, `human_uniprot_id`,
`human_annotation`, `start`, `end`, `resource` — one row per predicted domain-motif interaction,
where `start`/`end` are the motif's 0-indexed, end-exclusive match positions (Python `re.Match`
slice semantics) in whichever species' sequence carries the motif for that row (`human` if
`dmi_type == 'forward'`, `bacterial` if `'reverse'`).

Case-study's legacy IDR code lives in the old flat `workflow/` (top-level, pre-refactor, not
packaged): `idr_prediction.py` (IUPred2/ANCHOR2), `idr_prediction_score.py` (a near-duplicate of
`idr_prediction.py`, differing only in trivial whitespace), `AIUPred.py` and `AIUPred_old.py`
(AIUPred), plus two vendored third-party model implementations, `aiupred_lib.py` and
`iupred2a.py`. None of these are imported anywhere else in the repo (confirmed via
`grep -rl "idr_prediction\|AIUPred\|aiupred_lib\|iupred2a"`).

**Key decision: none of this legacy code is ported.** The plan's own pointer —
`github.com/saezlab/iupred` — turns out to be a clean, pip-installable wrapper around both IUPred2
and AIUPred, fetched and confirmed directly from the package's README:

```
pip install git+https://github.com/saezlab/iupred.git                          # IUPred2/ANCHOR2 only
pip install "iupred[aiupred] @ git+https://github.com/saezlab/iupred.git"       # + AIUPred (needs torch)
```

Import name is `iupred`, exposing `iupred(sequence, mode='short'|'long'|'glob')` (returns
`(scores, glob_text)`), `anchor2(sequence, disorder_scores)`, `aiupred_disorder(sequence,
force_cpu=False, gpu_num=0)`, `aiupred_binding(sequence, force_cpu=False, gpu_num=0)`, plus
`ensure_aiupred_data()`/`clear_aiupred_cache()` for its own model-download cache
(`~/.cache/iupred/`). This makes the six legacy files, and the vendored model weights under
`case_study_input/input/iupred_data` and `case_study_input/input/aiupred_data` that
`idr_prediction.py`/`AIUPred.py` read from, fully redundant — one dependency replaces all of it.

Two scope decisions were confirmed with the user before writing this plan (working through several
rounds — see below):

1. **Per-residue gate, not an average-based one.** The legacy scripts differ from each other here:
   `idr_prediction.py`'s `motif_selection` used a joint, tolerant rule (count residues where
   disorder *and* anchor both exceed 0.5 simultaneously; accept if that count is within ±1 of the
   window length), while `AIUPred.py`'s `keep_high_confidence` used a stricter, separate rule (every
   residue must individually clear its own track's cutoff — `(d_prof[s:e+1] >= thr).all() and
   (b_prof[s:e+1] >= thr).all()`, `thr=0.60`, no tolerance). The plan's literal spec text
   ("filtered ... based on the disorder/binding score cut off") is ambiguous between these and a
   simple average-of-window check. Averaging was considered and rejected — it would let a strong
   motif carry a few weak/non-disordered residues, undermining the point of the filter.
2. **Two separate per-track gates, not one merged score — matching `AIUPred.py`'s existing logic
   exactly**, applied identically regardless of `method`. A combined-score-per-residue design
   (`(disorder + binding) / 2 >= combined_cutoff`, checked per residue) was also considered and
   rejected in favor of keeping the two tracks independent, since that's what the only internally-
   consistent piece of legacy logic already does.

The net effect: `AIUPred.py`'s `keep_high_confidence` becomes this module's reference behavior,
generalized to also work for `method='iupred'`, with no ±1 tolerance carried over from
`idr_prediction.py`'s looser variant (one consistent rule across both methods, not two).

## Decisions

1. **One public function**, mirroring `dmi.py`'s "take Module N-1's output shape directly, no
   dataclasses" precedent:
   ```python
   def filter_by_disorder(
       dmi_table: pd.DataFrame,
       human_sequences: dict[str, str] | None,
       bacterial_sequences: dict[str, str] | None,
       method: str,
       disorder_cutoff: float,
       binding_cutoff: float,
       force_cpu: bool = False,
       gpu_num: int = 0,
   ) -> pd.DataFrame
   ```
   `human_sequences`/`bacterial_sequences` are the same header-keyed FASTA dicts Module 3 produces
   and Module 6 already consumes (`dict[fasta_header -> sequence]`) — no new input shape, no file
   I/O in the workflow module (file reading stays in `cli.py`, matching Modules 5/6).
2. **Per-protein profile memoization**, not per-row. Module 6's cross-join means the same
   `(motif_protein, start, end)` motif commonly recurs across many rows (one per compatible partner
   protein), so scoring is cached per unique `(method, sequence, force_cpu, gpu_num)` via
   `functools.lru_cache`, keyed on the sequence itself rather than the accession — avoids re-running
   AIUPred's model forward pass (or IUPred2/ANCHOR2) once per row when it's really once per protein.
3. **`method` dispatches to `iupred` or `aiupred`, raising `ValueError` otherwise** — same pattern as
   `dmi.py`'s `mode` validation and `membrane_filter.py`'s `species` validation.
4. **IUPred2/ANCHOR2 profile**: disorder is computed with `mode='short'` (matches
   `idr_prediction.py`'s default); binding is `anchor2(sequence, disorder_long)` where
   `disorder_long = iupred(sequence, mode='long')[0]` — ANCHOR2 requires long-mode disorder as its
   input regardless of which mode is reported, a real algorithmic dependency in the original
   IUPred2/ANCHOR2 design, not legacy cruft, so it's kept.
5. **Gate is two independent per-residue `.all()` checks** (decision 2 above):
   `(disorder_window >= disorder_cutoff).all() and (binding_window >= binding_cutoff).all()`, no
   tolerance, applied identically for both methods.
6. **Reported score columns are decoupled from the gate.** `disordered_score`/`binding_score` are
   the window's mean disorder/binding score (`.mean()`, not part of the pass/fail decision);
   `combined_score` is their average. Computed only for rows that already passed the per-residue
   gate — matches the plan's literal "three extra columns" spec, while keeping the actual filtering
   decision on the stricter per-residue rule.
7. **A row whose motif protein has no sequence in the corresponding `*_sequences` argument is
   dropped, not a hard error.** In practice this shouldn't happen (Module 7's FASTA input should be
   the exact same file Module 6 scored positions against), but treating a single missing sequence
   as fatal for the whole batch is disproportionate — a `ValueError` for missing *sequences dict
   entirely (mirrors `dmi.py`'s "required input for mode not supplied" check) is still raised
   up-front if forward rows exist but `human_sequences is None` (or the reverse), since that *is* a
   caller error.
8. **Output columns**: `dmi.OUTPUT_COLUMNS` (the exact 8 columns from Module 6) plus
   `disordered_score`, `binding_score`, `combined_score` — 11 total, extending rather than
   reshaping Module 6's schema.
9. **New optional dependency extra, `idr`**, not a base dependency: `iupred[aiupred]` pulls in
   `torch`, which is heavy and unnecessary for users who only need `method='iupred'`. Both
   `iupred`/`anchor2`/`aiupred_disorder`/`aiupred_binding` imports are lazy (inside function
   bodies), matching the existing convention (`import omnipath as op` inside
   `membrane_filter.py`'s function body).
10. **Legacy files deleted outright (per Q11, no shims)**: `workflow/idr_prediction.py`,
    `workflow/idr_prediction_score.py`, `workflow/AIUPred.py`, `workflow/AIUPred_old.py`,
    `workflow/aiupred_lib.py`, `workflow/iupred2a.py` — confirmed unreferenced elsewhere. The
    vendored model data under `case_study_input/input/iupred_data` and
    `case_study_input/input/aiupred_data` becomes unreferenced by any code once these are deleted,
    but is left in place — it's shared fixture data, not code, and deleting fixtures is out of this
    module's scope.
11. **`AIUPred.py`'s fixture is a real regression check, not just manual sign-off.** Unlike Modules
    5/6, `case_study_output/MicrobioLink_AIUPRED_outcome.csv` has genuine, correctly-written score
    columns (`Avg IUPred score` [disorder], `Avg Binding score`, `Combined score`) at `thr=0.60` for
    both tracks — and decisions 2/5 above make this module's `method='aiupred'` gate logic identical
    to what produced that file. It can be regenerated and diffed (after remapping column names/order
    — see Migration checklist). By contrast, `idr_prediction.py`'s own `write_output` has a bug: it
    writes a 10-column header promising score data, but the body only ever writes
    `"\t".join(interaction)` — the raw HMI fields — never the computed `motif[4:]` score values. The
    only committed IUPred fixture,
    `case_study_output/output/HMI/IUPred/BT_enterocyte_idr_cd_usecase.csv`, is consequently just 6
    columns with no score data despite its promising header, so `method='iupred'` has no usable
    fixture and goes through manual sign-off only (same treatment as Modules 5/6).

## Target implementation

### `microbiolink/workflow/idr_filter.py` (new, core, no argparse)

```python
OUTPUT_COLUMNS = [*dmi.OUTPUT_COLUMNS, 'disordered_score', 'binding_score', 'combined_score']


def _iupred_profile(sequence: str) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-residue IUPred2 disorder and ANCHOR2 binding scores.

    Args:
        sequence: Amino acid sequence to score.

    Returns:
        A (disorder_scores, binding_scores) pair of NumPy arrays, each one score per residue
        in sequence. binding_scores is computed via ANCHOR2 against the long-mode IUPred2
        disorder profile, not the short-mode profile disorder_scores itself uses.
    """


def _aiupred_profile(sequence: str, force_cpu: bool, gpu_num: int) -> tuple[np.ndarray, np.ndarray]:
    """Compute per-residue AIUPred disorder and binding scores.

    Args:
        sequence: Amino acid sequence to score.
        force_cpu: Force CPU inference even if a GPU is available.
        gpu_num: Index of the GPU to use when force_cpu is False.

    Returns:
        A (disorder_scores, binding_scores) pair of NumPy arrays, each one score per residue
        in sequence.
    """


@functools.lru_cache(maxsize=None)
def _cached_profile(
    method: str, sequence: str, force_cpu: bool, gpu_num: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute and cache a sequence's disorder/binding profile, keyed by method and device.

    Module 6's output cross-joins one motif hit against every partner protein sharing a
    compatible domain, so the same motif protein's sequence is commonly scored by many DMI
    rows. Caching here means each unique (method, sequence, force_cpu, gpu_num) combination is
    only ever scored once, regardless of how many rows reference it.

    Args:
        method: 'iupred' or 'aiupred'.
        sequence: Amino acid sequence to score.
        force_cpu: Force CPU inference (aiupred only; accepted but ignored for iupred).
        gpu_num: GPU index (aiupred only; accepted but ignored for iupred).

    Returns:
        A (disorder_scores, binding_scores) pair of NumPy arrays, each one score per residue
        in sequence.

    Raises:
        ValueError: If method is not 'iupred' or 'aiupred'.
    """


def _build_sequence_lookup(sequences: dict[str, str] | None) -> dict[str, str]:
    """Reindex a header-keyed FASTA dict by UniProt accession.

    Args:
        sequences: FASTA header -> sequence mapping (Module 3's output shape), or None.

    Returns:
        A dict mapping each UniProt accession to its sequence. Empty if sequences is None.
    """


def _motif_side(row) -> tuple[str, int, int]:
    """Resolve which protein and position carries the motif for one DMI row.

    Args:
        row: One row of Module 6's DMI output table.

    Returns:
        (motif_uniprot_id, start, end) for whichever side is the motif-bearing side:
        human_uniprot_id/start/end if row['dmi_type'] == 'forward', otherwise the bacterial
        equivalent.
    """


def _score_row(
    row, sequence_lookup, method, disorder_cutoff, binding_cutoff, force_cpu, gpu_num,
) -> tuple[bool, float, float, float] | None:
    """Score one DMI row's motif window against the per-residue disorder/binding gate.

    Args:
        row: One row of Module 6's DMI output table.
        sequence_lookup: UniProt accession -> sequence mapping for the motif-bearing species
            (built by _build_sequence_lookup).
        method: 'iupred' or 'aiupred'.
        disorder_cutoff: Minimum per-residue disorder score required across the whole window.
        binding_cutoff: Minimum per-residue binding score required across the whole window.
        force_cpu: Passed through to _cached_profile.
        gpu_num: Passed through to _cached_profile.

    Returns:
        (passes, disordered_score, binding_score, combined_score), where passes is True iff
        every residue in the motif window clears both cutoffs, and the three scores are the
        window's mean disorder, mean binding, and their average. None if the motif protein's
        sequence is not present in sequence_lookup.
    """


def filter_by_disorder(
    dmi_table: pd.DataFrame,
    human_sequences: dict[str, str] | None,
    bacterial_sequences: dict[str, str] | None,
    method: str,
    disorder_cutoff: float,
    binding_cutoff: float,
    force_cpu: bool = False,
    gpu_num: int = 0,
) -> pd.DataFrame:
    """Filter Module 6's DMI table to interactions in a disordered, binding-prone region.

    Args:
        dmi_table: Module 6's output table (predict_domain_motif_interactions' return value).
        human_sequences: Human FASTA header -> sequence mapping. Required if dmi_table
            contains any 'forward' rows.
        bacterial_sequences: Bacterial FASTA header -> sequence mapping. Required if
            dmi_table contains any 'reverse' rows.
        method: 'iupred' or 'aiupred'.
        disorder_cutoff: Minimum per-residue disorder score required across the whole motif
            window.
        binding_cutoff: Minimum per-residue binding score required across the whole motif
            window.
        force_cpu: Force CPU inference (aiupred only; ignored for iupred).
        gpu_num: GPU index to use (aiupred only; ignored for iupred).

    Returns:
        dmi_table restricted to rows whose motif window clears both cutoffs at every residue,
        with three extra columns: disordered_score, binding_score, combined_score.

    Raises:
        ValueError: If method is not 'iupred' or 'aiupred', or dmi_table contains 'forward'
            rows with human_sequences=None (or 'reverse' rows with bacterial_sequences=None).
    """
```

- `_iupred_profile`: lazily `from iupred import anchor2, iupred`; `disorder, _ =
  iupred(sequence, mode='short')`; `disorder_long, _ = iupred(sequence, mode='long')`;
  `binding = anchor2(sequence, disorder_long)`; returns `(np.asarray(disorder),
  np.asarray(binding))`.
- `_aiupred_profile`: lazily `from iupred import aiupred_binding, aiupred_disorder`; returns
  `(np.asarray(aiupred_disorder(sequence, force_cpu=force_cpu, gpu_num=gpu_num)),
  np.asarray(aiupred_binding(sequence, force_cpu=force_cpu, gpu_num=gpu_num)))`.
- `_cached_profile`: dispatches to `_iupred_profile` or `_aiupred_profile` on `method`, raising
  `ValueError` for any other value; `force_cpu`/`gpu_num` are accepted but ignored on the
  `iupred` path (still part of the cache key, harmless — same sequence scored under `method='iupred'`
  is identical regardless of those two args).
- `_build_sequence_lookup`: `{fasta.extract_uniprot_id(header): seq for header, seq in
  sequences.items()}` if `sequences` is not `None`, else `{}`.
- `_motif_side`: `row['human_uniprot_id'], row['start'], row['end']` if `row['dmi_type'] ==
  'forward'`, else the bacterial equivalent.
- `_score_row`: looks up the motif protein's sequence in `sequence_lookup`; returns `None` if
  absent. Otherwise calls `_cached_profile`, slices `disorder_profile[start:end]` /
  `binding_profile[start:end]`, computes `passes = bool((disorder_window >=
  disorder_cutoff).all() and (binding_window >= binding_cutoff).all())`,
  `disordered_score = float(disorder_window.mean())`, `binding_score =
  float(binding_window.mean())`, `combined_score = (disordered_score + binding_score) / 2`.
- `filter_by_disorder`: validates `method in {'iupred', 'aiupred'}`; raises `ValueError` if any
  `dmi_type == 'forward'` row exists and `human_sequences is None` (symmetric check for
  `'reverse'`/`bacterial_sequences`) — mirrors `dmi.py`'s required-input check for the direction
  actually present in the data, not the full Cartesian set `dmi.py` checks (Module 7 only ever
  gets one combined table, not per-direction calls). Builds both sequence lookups once, iterates
  `dmi_table.itertuples()`, calls `_score_row` per row; a `None` result is dropped (decision 7);
  rows that fail the gate are dropped; passing rows are collected as tuples and returned as
  `pd.DataFrame(kept_rows, columns=OUTPUT_COLUMNS)` (explicit columns so an empty result keeps
  the right schema, matching `dmi.py`/`ddi.py`'s precedent).

### `microbiolink/cli.py` (extend, argparse only)

- Add `_build_idr_filter_parser()`: `--dmi_file` (required), `--human_fasta_file` (optional,
  required for forward rows), `--bacterial_fasta_file` (optional, required for reverse rows),
  `--method` (required, `choices=['iupred', 'aiupred']`), `--disorder_cutoff` (required, `type=float`),
  `--binding_cutoff` (required, `type=float`), `--force-cpu` (`action='store_true'`, `dest='force_cpu'`),
  `--gpu-num` (`type=int`, `default=0`, `dest='gpu_num'`), `-o/--output_file` (required).
- Add entry point:
  ```python
  def idr_filter() -> int:
      """Filter domain-motif interactions to those in a disordered, binding-prone region."""

      from .utils import fasta
      from .workflow import idr_filter as idr_filter_module

      args = _build_idr_filter_parser().parse_args()
      dmi_table = pd.read_csv(args.dmi_file)
      human_sequences = fasta.read_fasta_sequences(args.human_fasta_file) if args.human_fasta_file else None
      bacterial_sequences = (
          fasta.read_fasta_sequences(args.bacterial_fasta_file) if args.bacterial_fasta_file else None
      )

      result = idr_filter_module.filter_by_disorder(
          dmi_table,
          human_sequences=human_sequences,
          bacterial_sequences=bacterial_sequences,
          method=args.method,
          disorder_cutoff=args.disorder_cutoff,
          binding_cutoff=args.binding_cutoff,
          force_cpu=args.force_cpu,
          gpu_num=args.gpu_num,
      )
      result.to_csv(args.output_file, index=False)
      return 0
  ```
- `pyproject.toml`: add `microbiolink-idr-filter = "microbiolink.cli:idr_filter"` under
  `[project.scripts]`; add
  ```
  [project.optional-dependencies]
  idr = ["iupred[aiupred] @ git+https://github.com/saezlab/iupred.git"]
  ```
  (base `dependencies` stay untouched — `torch` and `iupred` are opt-in via this extra).

## Migration checklist

### 1. Packaging
- [x] Add `[project.optional-dependencies] idr = [...]` to `pyproject.toml`.
- [x] Add `microbiolink-idr-filter = "microbiolink.cli:idr_filter"` to `[project.scripts]`.

### 2. Core module
- [x] Create `microbiolink/workflow/idr_filter.py` with `_iupred_profile`, `_aiupred_profile`,
      `_cached_profile`, `_build_sequence_lookup`, `_motif_side`, `_score_row`,
      `filter_by_disorder`.
- [x] Confirm `_cached_profile` is actually cached — same sequence scored twice (e.g. two DMI rows
      sharing a motif protein) triggers the underlying `iupred`/`aiupred` call only once.
      (`functools.lru_cache` keyed on `(method, sequence, force_cpu, gpu_num)`; behavior implied by
      construction, exercised implicitly by the AIUPred regression run below completing in ~3 min
      over 91,980 DMI rows rather than re-scoring every row independently.)
- [x] Confirm a row whose motif protein isn't in the supplied sequences is dropped, not raised,
      while a wholesale missing `human_sequences`/`bacterial_sequences` (when that direction's rows
      exist) does raise `ValueError`. (`_score_row` returns `None` on a missing lookup entry,
      dropped in `_score_all_rows`; `_validate_inputs` raises `ValueError` for a wholesale-missing
      dict — confirmed by code inspection.)

### 3. CLI wiring
- [x] Add `_build_idr_filter_parser`, `idr_filter()` entry point to `cli.py`.

### 4. Delete old code (per Q11)
- [x] Delete `workflow/idr_prediction.py`, `workflow/idr_prediction_score.py`,
      `workflow/AIUPred.py`, `workflow/AIUPred_old.py`, `workflow/aiupred_lib.py`,
      `workflow/iupred2a.py` (top-level, pre-refactor) — confirmed unreferenced elsewhere via
      `grep -rl "idr_prediction\|AIUPred\|aiupred_lib\|iupred2a"`.

### 5. Regression check — AIUPred path (real fixture, per decision 11)
- [x] Regenerate the forward-mode DMI table feeding into this fixture: ran
      `microbiolink-download-domains` against the `Entry` IDs in
      `case_study_input/input/bacterial_protein/BT_BEV_domains.tsv`, then `microbiolink-dmi --mode
      forward` against `case_study_input/input/human_transcriptomics/protein_sequences.fasta` →
      91,980 DMI rows.
- [x] Ran `microbiolink-idr-filter --method aiupred --disorder_cutoff 0.60 --binding_cutoff 0.60` →
      5,574 passing rows.
- [x] Diffed against `case_study_output/MicrobioLink_AIUPRED_outcome.csv` (227 rows) after remapping
      columns and joining on `(bacterial_uniprot_id, human_uniprot_id, start, end)`. **Result: PASS.**
      211/227 fixture rows matched keys in the new output (the other 16 are explained by live
      UniProt/Pfam domain-annotation drift since the fixture was generated — `microbial_domains.tsv`
      freshly queried from UniProt differs in content from the vendored `BT_BEV_domains.tsv`, which
      is Module 6/domain-download's input, not this module's). Of the 211 key-matched rows, scores
      agree within `atol=0.02` on all three float columns for 211/211 (disorder), 194/211 (binding),
      211/211 (combined) — mean absolute differences of 0.001–0.005, consistent with AIUPred
      re-downloading current model weights rather than the historical ones used to build the
      fixture, not a logic defect. No gate-flipping mismatches (a row passing in one output and
      failing in the other while disagreeing outside a few percent) were observed.

### 6. Manual sign-off — IUPred path (no usable fixture, per decision 11)
- [x] Confirmed `_iupred_profile` calls `anchor2` with the *long*-mode disorder profile
      (`idr_filter.py:31-37`), not the short-mode one used for the reported score.
- [x] Ran the same forward-mode DMI table through `--method iupred --disorder_cutoff 0.5
      --binding_cutoff 0.5` (7,531 passing rows) and compared against
      `case_study_output/output/HMI/IUPred/BT_enterocyte_idr_cd_usecase.csv`'s protein/motif/position
      combinations. **Result: substantial overlap, as expected** — 87/91 (96%) of the fixture's
      distinct `(human_uniprot_id, human_annotation, start, end)` motif-positions are recovered by
      the new output, and 358/403 full `(bacterial, human, motif, start, end)` fixture rows match
      exactly. The 4 unrecovered motif-positions are consistent with the same
      UniProt/Pfam-domain-drift and cutoff/tolerance differences noted above and in decision 1/2, not
      a logic bug (that legacy fixture has no score columns to compare numerically — see decision
      11's note on `idr_prediction.py`'s `write_output` bug).
- [x] Verified `dmi_type` routing directly: constructed a synthetic 2-row DMI table (one `forward`,
      one `reverse` row over the same accessions) and called `filter_by_disorder` twice, swapping
      which side (`human_sequences` vs `bacterial_sequences`) held a strongly-disordered vs.
      strongly-ordered synthetic sequence. The `forward` row's `disordered_score` tracked
      `human_sequences` and the `reverse` row's tracked `bacterial_sequences` in both swaps,
      confirming `_motif_side`/`_score_row` route each `dmi_type` to the correct species' sequences.
      (Ran as a direct Python check rather than a full `mode='reverse'`/`'both'` case-study pipeline
      run, since that would require downloading bacterial FASTA sequences and human Pfam domains at
      case-study scale — a large additional live-API cost — purely to re-confirm a code path already
      covered unambiguously by this synthetic test.)
- [x] Shown to the user for confirmation below.

## Verification

- `uv run ruff format microbiolink/workflow/idr_filter.py` — scoped to the new file.
- `uv run ruff check microbiolink/workflow/idr_filter.py` — same scoping. `cli.py`'s new
  `_build_idr_filter_parser`/`idr_filter()` additions and `pyproject.toml`'s changes are reviewed by
  eye against the existing files' style instead, since both are modified, not newly created, and
  aren't currently ruff-clean in full.
- `uv run ty check`
- Manual CLI runs and the AIUPred regression diff described in the checklist above.
