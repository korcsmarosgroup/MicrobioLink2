# Module 1 — Z-score Filter — Implementation Plan

Detailed plan for the first module described in @plans/microbiolink_refactoring.md, following
the decisions in @plans/refactoring_questions.md. Scope is exactly Module 1: filtering a gene
count matrix by a z-score cutoff.

## Source of truth

Three existing copies of this logic were compared:

- `case-study:workflow/z-score_filter_terminal.py` — **ground truth** (ties to
  `plans/refactoring_questions.md` Q11: case-study always wins over beta on overlapping logic).
- `MicrobioLink-2.1-beta:microbiolink/z-score_filter_terminal.py` — byte-identical copy of the
  above (stale duplicate, not a source of new logic).
- `MicrobioLink-2.1-beta:microbiolink/z_score_filter_terminal.py` — same algorithm, cleaned up
  into `main()` / `parse_args()` / `cli_main()`, no behavior change.
- `MicrobioLink-2.1-beta:microbiolink_api/expression.py` — same core algorithm, restructured as
  `filter_counts_by_zscore()` (DataFrame → DataFrame) / `filter_count_matrix_file()` (file I/O
  wrapper), with added edge-case guards.

All four compute the identical statistic per column: fit a Gaussian KDE to the non-NaN values,
take the mode (`mu`) as the peak of the KDE, take the mean of values above `mu` (`U`), derive
`sigma = (U - mu) * sqrt(pi/2)`, then `z = (value - mu) / sigma`. A value is kept as-is if
`z > cutoff`, otherwise it becomes NaN.

The legacy script's `count[list(zcount).index(x)] if x > args.zscore else 'NaN'` construct looks
like a bug (index-based lookup instead of positional), but is confirmed harmless: ties in
`zcount` only occur for ties in the original `count` value (the z transform is affine per
column), so the "wrong" index still returns a numerically identical value. NaN inputs never
reach `.index()` because the ternary's condition (`x > cutoff`) short-circuits to the else
branch first. **No behavior change here** — just replace it with a direct positional/vectorized
computation in the rewrite, since the index-lookup is O(n²) and confusing.

## Decisions from grilling

1. **Output NaN representation**: real `NaN` (empty CSV cell via pandas' default `to_csv`), not
   the legacy literal string `'NaN'`. The legacy string format exists only because downstream
   `get_proteins()` in `get_human_fasta.py` string-compares `line[1] != 'NaN'` — that consumer is
   itself being rewritten onto pandas as part of this same refactor (Module 2/3), so there is no
   remaining reason to keep the string convention. Note this as a deliberate, confirmed
   deviation from the byte-for-byte legacy CSV output.
2. **No new edge-case guards.** `microbiolink_api/expression.py` added guards for <2 non-NaN
   values, zero variance, empty upper tail, and `LinAlgError` from a degenerate KDE fit — none of
   which exist in the legacy script. **Decision: do not port these.** Module 1 matches the legacy
   script's behavior exactly, including the fact that a degenerate column (e.g. constant values,
   or fewer than 2 non-NaN entries) will raise rather than silently pass through unfiltered. If
   this becomes a real problem against actual data it should be handled as a deliberate follow-up
   change, not folded silently into this port.
3. **Boundary semantics**: strict `z > cutoff` to keep a value (matches legacy exactly). A value
   with `z == cutoff` is filtered to NaN. This narrows the ambiguous prose in the top-level plan
   ("NaN unless below cut-off") to match actual ground-truth behavior.
4. **File naming**: new file `microbiolink/zscore_filter.py`, functions `filter_counts_by_zscore()`
   and `filter_count_matrix_file()` (names carried over from `microbiolink_api/expression.py` —
   only the guards are dropped, not the naming/structure). This retires both of beta's ambiguous
   `z-score_filter_terminal.py` / `z_score_filter_terminal.py` names.

## Target implementation

### `microbiolink/zscore_filter.py` (core, no argparse)

```python
def filter_counts_by_zscore(
    count_matrix: pd.DataFrame,
    zscore_threshold: float = -3,
) -> pd.DataFrame:
    """Filter a count matrix by the MicrobioLink z-score rule."""
```

- Operates column-by-column exactly as the legacy script does (each sample/condition column
  gets its own KDE fit, independent of other columns).
- Per column: coerce to float, fit `gaussian_kde` on non-NaN values, locate `mu` via
  `np.linspace(min, max, 100)` + `kernel.evaluate` + `argmax` (same 100-point grid as legacy —
  do not change the resolution, since it's an unstated legacy parameter that could shift boundary
  values), compute `sigma`, compute `z`, keep original value where `z > zscore_threshold`, else
  `NaN`.
- Returns a new DataFrame, same shape/index/columns as input, same dtype rules as legacy (numeric
  in, numeric-or-NaN out).
- No file I/O, no argparse, no `print`. Pure function suitable for library/API use per the
  top-level plan's "public functions should be enough to run microbiolink without the cli."

```python
def read_count_matrix(filename: PathLike, index_col: int | str = 0) -> pd.DataFrame:
    """Read a gene/protein count matrix from disk."""


def filter_count_matrix_file(
    input_file: PathLike,
    zscore_threshold: float = -3,
    output_file: PathLike | None = None,
) -> pd.DataFrame:
    """Read, filter, and optionally write a count matrix."""
```

- `filter_count_matrix_file` is the file-based convenience wrapper the CLI calls into; also
  usable directly as a library function without ever touching argparse.
- `zscore_threshold` is typed/accepted as `float`, not `int` — the legacy argparse forced
  `type=int`, which is almost certainly an oversight (z-score cutoffs are not naturally integral,
  e.g. `-2.5`). Accepting `float` is a strict superset of accepting `int` (whole-number cutoffs
  like `-3` still work identically), so this isn't a behavior change for any input that worked
  before, and it removes a silent-truncation gotcha (e.g. `-2.5` used to silently truncate to
  `-2` in the legacy script — this rewrite requires no truncation).

### `microbiolink/cli.py` (argparse only, per Q2)

Add:

```python
def zscore_filter() -> int:
    from . import zscore_filter as zscore_filter_module

    parser = argparse.ArgumentParser(
        description='Gene expression filtration based on individual cell count and z-score.',
    )
    parser.add_argument('-i', '--input_file', required=True, help='Input CSV file with gene expression data.')
    parser.add_argument('-zscore', '--zscore', required=True, type=float, help='Z-score cut-off to filter lowly expressed genes.')
    parser.add_argument('-o', '--output_file', required=True, help='Output CSV file for filtered results.')
    args = parser.parse_args()

    zscore_filter_module.filter_count_matrix_file(
        args.input_file,
        zscore_threshold=args.zscore,
        output_file=args.output_file,
    )
    return 0
```

- No `default=-3` on `--zscore` in the CLI parser (the legacy script declared a default but also
  `required=True`, which makes the default dead code — drop the unreachable default, keep
  `required=True`).
- `pyproject.toml` entry point: `microbiolink-zscore-filter = "microbiolink.cli:zscore_filter"`
  (renamed from beta's `microbiolink-zscore-filter = "microbiolink.cli:z_score_filter_terminal"`
  to match the new function name — the CLI command name itself is unchanged for users).

## Inputs/Outputs recap (from top-level plan)

- **Input**: gene count matrix (CSV, genes as rows, samples/conditions as columns) + a
  user-defined z-score cutoff.
- **Output**: same-shape count matrix; each value is either the original count or `NaN` if its
  column-wise z-score does not exceed the cutoff.

## Migration checklist

### 1. Core module

- [ ] Create `microbiolink/zscore_filter.py`.
- [ ] Port `read_count_matrix()` from `microbiolink_api/expression.py`.
- [ ] Port `filter_counts_by_zscore()`, replacing the legacy `list(zcount).index(x)` lookup with
      direct positional computation (no behavior change, see "Source of truth").
- [ ] Confirm no edge-case guards are added — degenerate columns (< 2 non-NaN values, zero
      variance, empty upper tail, singular KDE) must raise, matching legacy behavior (decision 2).
- [ ] Port `filter_count_matrix_file()`, with `zscore_threshold: float` (not `int`).
- [ ] Confirm output NaN cells are real `NaN` (pandas default `to_csv` empty-cell behavior), not
      the literal string `'NaN'` (decision 1).

### 2. CLI wiring

- [ ] Add `zscore_filter()` to `microbiolink/cli.py` (argparse only — no `parse_args()` left in
      the core module, per Q2).
- [ ] `--zscore` argument is `type=float`, `required=True`, no `default` (the legacy default was
      dead code alongside `required=True`).
- [ ] Update `pyproject.toml` entry point to
      `microbiolink-zscore-filter = "microbiolink.cli:zscore_filter"`.

### 3. Regression check (per Q11)

Must run **before** section 4 deletes `workflow/z-score_filter_terminal.py` — the baseline needs
the old script to exist. Pull the baseline script explicitly from the `case-study` branch (the
documented ground truth), rather than relying on the `refactoring` branch's working-tree copy,
so this step stays correct even after that copy is deleted or if the two branches ever diverge:

- [ ] `git show case-study:workflow/z-score_filter_terminal.py > /tmp/baseline_zscore_filter.py`
      (or equivalent) to get an unambiguous ground-truth copy, independent of what's currently
      checked out on `refactoring`.
- [ ] Run that baseline script on `case_study_input/input/enterocyte_colon_CD_zscore.csv` at a
      chosen cutoff (e.g. `-3`) to produce a baseline output. (No pre-built fixture exists for
      this module — the file that looks like one,
      `case_study_input/input/human_transcriptomics/enterocyte_colon_CD_zscore.csv`, has an
      unknown/unrecorded cutoff baked in, so it can't be used directly.)
- [ ] Run the **new** `filter_count_matrix_file()` on the same raw input at the same cutoff.
- [ ] Diff old vs. new: same cells are NaN/`'NaN'` in both (accounting for the intentional
      representation difference from decision 1).
- [ ] Diff old vs. new: non-NaN numeric values match within float tolerance (KDE fit is
      floating-point, not expected to be bit-exact).

### 4. Delete old code (no shims, per Q11)

This refactor happens on the `refactoring` branch, which only contains case-study's files —
the beta files (`microbiolink/z-score_filter_terminal.py`, `microbiolink/z_score_filter_terminal.py`,
`microbiolink_api/expression.py`) were read via `git show origin/MicrobioLink-2.1-beta:...` purely
as reference for comparing implementations; they were never checked out here, so there is nothing
to delete for them on this branch.

- [ ] Delete `workflow/z-score_filter_terminal.py` (the only old copy of this logic present on
      the `refactoring` branch) — only after section 3's regression check has passed.

### 5. Manual sign-off

- [ ] Run the new CLI once end-to-end and inspect the output file by eye.
- [ ] Confirm the `float`-typed cutoff (vs. legacy `int`) and the dropped unreachable default
      behave as expected — these are small, intentional deviations from the legacy argparse
      behavior, not covered by the regression diff.
