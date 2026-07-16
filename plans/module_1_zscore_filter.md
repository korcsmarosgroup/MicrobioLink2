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

## Migration steps

1. Create `microbiolink/zscore_filter.py` with `read_count_matrix`, `filter_counts_by_zscore`,
   `filter_count_matrix_file`, ported from `microbiolink_api/expression.py` minus the edge-case
   guards (decision 2), with the index-lookup replaced by direct positional computation.
2. Add `zscore_filter()` to `microbiolink/cli.py`, wire up `pyproject.toml` entry point.
3. Delete `workflow/z-score_filter_terminal.py` (case-study) and both
   `microbiolink/z-score_filter_terminal.py` / `microbiolink/z_score_filter_terminal.py` and
   `microbiolink_api/expression.py` (beta) — no shims, per Q11.
4. Regression check (per Q11): no pre-built "before/after" fixture exists for this module
   specifically — `case_study_input/input/enterocyte_colon_CD_zscore.csv` is raw unfiltered
   input, and `case_study_input/input/human_transcriptomics/enterocyte_colon_CD_zscore.csv`
   (semicolon-delimited, legacy string-`'NaN'` format) is *an* old filtered output but the
   z-score cutoff used to produce it is not recorded anywhere. So: run the **old**
   `workflow/z-score_filter_terminal.py` and the **new** `filter_count_matrix_file()` on the same
   raw input (`case_study_input/input/enterocyte_colon_CD_zscore.csv`) at a chosen cutoff (e.g.
   `-3`, the legacy script's stated-but-unreachable default), then diff: same cells NaN in both,
   identical non-NaN numeric values (allow float tolerance from the KDE fit, not exact string
   match, since output NaN representation intentionally differs per decision 1).
5. Manually confirm the `float` vs `int` cutoff type change and dropped-default behavior look
   correct by running the new CLI once and inspecting output, since these are small intentional
   deviations from legacy argparse behavior.
