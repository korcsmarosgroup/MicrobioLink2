# Merge `MicrobioLink-2.1-beta` into `main` — Stage 1 (non-reverse-DMI)

**Status: executed.** Committed locally as `532b89a` on branch
`merge-beta-into-main`. `main` was not touched, nothing was pushed. This
document is the plan as originally written, with a **Deviations from the
plan** section documenting where execution differed from what was planned.

## Context

`main` and `origin/MicrobioLink-2.1-beta` had diverged significantly. A
direct `git merge` produced 10 conflicting files. Two of those
(`microbiolink/DMI.py`, `microbiolink_api/dmi.py`) stem from **two
independent, incompatible reverse-DMI implementations** — main's standalone
`reverse_DMI.py` / `reverse_dmi()` CLI command vs. beta's unified
`mode={forward,reverse,both}` flag and `BidirectionalDomainMotifInteraction`
schema baked into `dmi.py`. Reconciling that is a real design decision, not a
mechanical merge, and the user explicitly deferred it: **no reverse-DMI code
was merged, no new reverse-DMI tests were written, in this pass.**

This plan covered everything else: setting up an integration branch, merging
in beta's genuinely new/independent work (DDI module, tooling/CI/docs
scaffolding), removing the PHISTO benchmark bundle (per prior decision), and
resolving every conflict that is *not* reverse-DMI-entangled — while proving
main's existing behavior (everything except reverse-DMI) still works
afterward.

Investigation done before execution: full conflict list from a test merge
(aborted, no trace left), diffs of every conflicting file, and confirmation
that main's `microbiolink_api` package has **zero** reverse-DMI surface (it
only lives in the separate `microbiolink/reverse_DMI.py` script +
`reverse_dmi()` CLI command) — so main's DMI code could be left completely
untouched while DDI was merged in alongside it.

Confirmed with the user before execution:

- No `pixi.toml` exists on either branch (AGENTS.md's pixi mandate is
  aspirational for this repo currently) — verification used the existing
  `.venv` (Python 3.12) / `uv.lock` setup, e.g. `.venv/bin/python -m pytest`.
- Adopt beta's lazy `omnipath`/`mygene` import fix and shared
  `_omnipath_import_error` helper in `cli.py` — it's an independent
  robustness fix, not reverse-DMI-related.
- **Fully local only.** No push to `origin`, no draft branch, no PR.
- **Baseline established:** `.venv/bin/python -m pytest tests/ -q` on clean
  `main` → **40 passed**, 0 failures. Any post-merge failure would be
  attributable to the merge.
- **Version:** bump to `0.0.2` everywhere — `pyproject.toml` `version`,
  `microbiolink/_metadata.py` `_FALLBACK_VERSION`,
  `microbiolink_api/__init__.py` `__version__` literal, and
  `.bumpversion.cfg` `current_version`. (Beta's own `2.1.0` was internally
  inconsistent with its own `.bumpversion.cfg`, which said `0.0.1` — clearly
  an unfinished bump, not a deliberate release choice.)
- **Dev dependencies:** union of main's `pytest>=8.0` and beta's `distlib`,
  `build`, `pre-commit`, `bump2version`, `twine` — needed to actually run the
  tooling being adopted.
- **numpy pin:** drop the `<2` upper bound from the base `[project]
  dependencies` list (follow beta's split); keep `numpy>=1.26.4,<2` in the
  `idr` extra only (torch/iupred compat).
- **Python version floor:** bump `requires-python` to `>=3.10` and trim the
  adopted CI matrix to `['3.10', '3.11', '3.12', '3.13']` (drop `3.9`).
  Reason: main's own code already uses `PathLike = str | Path` (PEP 604) in
  `expression.py`, `microbiome.py`, `dmi.py`, `workflows.py`, and
  `motif_monte_carlo_filter.py` — this syntax raises `TypeError` on Python
  3.9 at import time, a pre-existing latent mismatch between
  `requires-python` and the code, discovered during this investigation and
  not introduced by the merge.
- **Metadata formatting:** take beta's author/maintainer/classifier
  formatting wholesale (identical content to main's, plus the additive
  `Development Status :: 2 - Pre-Alpha` classifier main doesn't have).
- **Verified safe:** `benchmarks/phisto_extended/` removal has zero
  references anywhere in beta's new tests or resource-building tools — safe
  to delete outright.
- **Verified safe:** `motif_monte_carlo_filter.py`, `reverse_DMI.py`,
  `tests/test_reverse_DMI.py`, the `tutorials/mock_data/*.ipynb` notebooks,
  and `plans/reverse_microbiolink.md` did not exist at the merge-base commit
  (`74a93bb`) — they were added to `main` after beta's fork point, so git's
  3-way merge auto-keeps them with no conflict and no risk of silent
  deletion.
- **Final check:** after the full test suite passes on the integration
  branch, also do a manual smoke check — `import microbiolink_api` and run
  the `ddi`/`dmi`/`reverse-dmi` CLI entry points with `--help`.

## Conflict classification (from test-merge investigation)

<!-- markdownlint-disable MD013 -->
| File | Nature | Resolution |
| --- | --- | --- |
| `.gitignore` | Non-overlapping ignore patterns from both sides | Combine both |
| `microbiolink_api/expression.py` | Only diff is `str \| Path` vs `Union[str, Path]` | Keep main's version (`ours`) |
| `microbiolink_api/microbiome.py` | Same as above | Keep main's version (`ours`) |
| `microbiolink_api/workflows.py` | PathLike style plus a `motif_sources` computation tied to beta's reverse-mode filtering, unsupported by main's `dmi.py` | Keep main's version entirely (`ours`) |
| `microbiolink/get_human_fasta.py` | Beta moves `omnipath`/`mygene` imports inside functions (lazy import), independent robustness fix | Take beta's version |
| `microbiolink/cli.py` | Mixed: new commands (`ddi()`, relocated existing ones), `_omnipath_import_error` helper, and reverse-DMI-only additions to `dmi()` (`--mode`, `--bacterial_fasta_file`, `--human_domain_file`) | Take the new-commands/helper parts; drop the reverse-DMI parts |
| `microbiolink_api/__init__.py` | add/add: main's exports are DMI-only; beta's are DMI (incl. bidirectional) + new DDI exports | Start from main's exports, add beta's DDI-only exports, skip all `Bidirectional*`/reverse-DMI exports |
| `pyproject.toml` | add/add across metadata, dependencies, `[project.scripts]`, and wholesale new `[tool.*]` sections | Rebuilt by hand per the confirmed decisions above |
| `microbiolink/DMI.py` | Reverse-DMI implementation clash | Deferred — kept main's version untouched |
| `microbiolink_api/dmi.py` | Reverse-DMI implementation clash | Deferred — kept main's version untouched, **except** one alias (see Deviations) |
<!-- markdownlint-enable MD013 -->

Additionally (no conflicts, purely additive from beta — brought in wholesale):

- `microbiolink/DDI.py`, `microbiolink_api/ddi.py`,
  `microbiolink_api/resources/*.tsv` (3did/DOMINE bundles),
  `tools/build_3did_dmi_resources.py`, `tools/build_domine_ddi_resources.py`,
  `tests/test_ddi_api.py`
- Tooling/scaffolding: `.pre-commit-config.yaml`, `.cruft.json`,
  `CODING_STYLE.md`, `mkdocs.yml`, `docs/`, `.github/`, `.bumpversion.cfg`,
  `.codecov.yaml`, `.coveragerc`, `.editorconfig`, `.python-version`

Removed entirely:

- `benchmarks/` (both `phisto_extended/` and `benchmarks/README.md`)

## Steps (as planned)

1. Set up the integration branch `merge-beta-into-main` off `main`, start the
   merge with `git merge --no-commit --no-ff origin/MicrobioLink-2.1-beta`.
2. Resolve trivial conflicts (`ours`): `expression.py`, `microbiome.py`,
   `workflows.py`.
3. Resolve `.gitignore` by hand-merging both pattern lists.
4. Resolve `get_human_fasta.py` by taking beta's version.
5. Resolve `cli.py` by hand: merge in `ddi()` + helper; drop reverse-DMI CLI
   additions.
6. Resolve `__init__.py` by hand: merge in DDI-only exports; set version to
   `0.0.2`.
7. Resolve `pyproject.toml` by hand: apply all confirmed decisions.
8. Leave `DMI.py` and `microbiolink_api/dmi.py` as main's version.
9. Remove the PHISTO benchmark bundle.
10. Sanity-check — no leftover conflict markers; reverse-DMI files unchanged;
    DDI files present.
11. Update `.bumpversion.cfg`, `microbiolink/_metadata.py` to `0.0.2` as well.
12. Install and run the full test suite; must pass for all pre-existing
    tests plus `test_ddi_api.py`.
13. Manual smoke check: import `microbiolink_api`, run `ddi`/`dmi`/
    `reverse-dmi` CLI `--help`.
14. Do not commit to `main` or push anything — leave on
    `merge-beta-into-main`, report back what was resolved and test results.

## Deviations from the plan

Things that came up during execution that the plan did not anticipate, or
that were resolved differently than originally written:

1. **Tooling scaffolding was not a clean additive merge — it required
   restoring deliberately-deleted files.** The plan assumed
   `.pre-commit-config.yaml`, `.cruft.json`, `mkdocs.yml`, `CODING_STYLE.md`,
   `.github/`, `docs/`, `.bumpversion.cfg`, etc. were beta-only additions
   with "no main equivalent." Mid-execution, investigation showed these
   files **existed at the merge-base and were deliberately deleted from
   main** in a commit titled `"Publish clean package-only MicrobioLink"`
   (`d09f3f6`) — the same commit that rewrote `microbiolink_api` into its
   current form. Because beta never touched these files after the fork,
   git's 3-way merge silently respected main's deletion (no conflict, files
   just vanished from the merge result). This meant re-adding them required
   manually restoring each file from beta's tree
   (`git show origin/MicrobioLink-2.1-beta:<path> > <path>` + `git add`)
   rather than relying on a clean automatic merge. The user was asked
   whether to still restore this tooling given it reversed a deliberate
   main decision, and chose to bring it back in.

2. **Borderline files were excluded from the restoration, narrowing the
   "bring in wholesale" list.** Of the files deleted in that same commit,
   the following were deliberately left out (not in the original plan's
   scope at all, decided live): `microbiolink/microbiolink_env.yml` and
   `microbiolink/requirements.txt` (superseded by `pyproject.toml` +
   `uv.lock`), `microbiolink/z-score_filter_terminal.py` (old hyphenated
   name — main already renamed it to `z_score_filter_terminal.py`),
   `tests/test_placeholder.py` (empty template placeholder), and
   `scripts/test-multi-py.py` (multi-python-version test runner script).

3. **A pre-existing latent bug was fixed as a side effect.** Not in the
   original conflict table: `requires-python` bumped from `>=3.9` to
   `>=3.10`, and the restored CI matrix trimmed to drop `3.9`, because
   main's own code already used Python 3.10-only syntax
   (`PathLike = str | Path`) while claiming 3.9 support. This was discovered
   only because bringing back the CI workflow would have made the mismatch
   visible on first push.

4. **`microbiolink_api/dmi.py` was not left 100% untouched — one alias line
   was added.** The plan's rule was "leave `dmi.py` exactly as main's, no
   reverse-DMI merge." During `__init__.py` resolution, discovered that
   `microbiolink_api/ddi.py` (which *is* being merged in) imports
   `read_protein_domain_table` from `dmi.py` — a function that only exists
   in beta's rewritten `dmi.py`, not main's. Comparing implementations
   showed it is byte-identical in logic to main's existing
   `read_bacterial_domain_table` (just organism-agnostic naming, since DDI
   reads domain tables for both bacterial and human proteins). Resolved by
   adding a single alias line to main's `dmi.py`:
   `read_protein_domain_table = read_bacterial_domain_table` — no reverse-
   DMI logic, no duplicated function body, confirmed via
   `git diff main -- microbiolink_api/dmi.py` to be the *only* change in
   that file relative to main.

5. **`tests/test_microbiolink_api.py` was discovered and excluded — not
   mentioned anywhere in the original plan.** This beta test file mixed 10
   legitimate tests for `workflows.py`/`microbiome.py`/forward `dmi.py`
   (filling a real gap: main has zero test coverage for these modules today)
   with 3 reverse-DMI-specific tests (`test_predict_reverse_domain_motif_
   interactions_returns_records`, `..._filters_clv_motifs`,
   `test_predict_bidirectional_domain_motif_interactions_combines_modes`)
   depending on the excluded API. The user chose to exclude the whole file
   for this pass rather than strip it down, deferring to better, more
   deliberate test coverage later.

6. **`tests/test_package_smoke.py` was brought in but not named in the
   original plan**, which only listed `test_ddi_api.py` as the new test
   file. It was verified safe (no reverse-DMI imports, and its
   `cli.dmi()` end-to-end test matches main's actual `DMI.py` output format
   exactly) and kept, contributing 9 of the final 14 net-new passing tests.

7. **`pip` was unavailable in the `.venv`.** Step 12 called for
   `.venv/bin/python -m pip install -e .`; there was no `pip` module in this
   environment. Used `uv pip install -e . --python .venv/bin/python`
   instead, which also fixed a stale `microbiolink-0.0.1.dist-info` that had
   caused `microbiolink_api.__version__` to initially report `0.0.1` instead
   of `0.0.2` (source files were correct all along; only the installed
   package metadata was stale).

8. **Final test count differs from the plan's description.** The plan said
   "all pre-existing tests plus `test_ddi_api.py`." Actual final count:
   **54 passed** (40 baseline unchanged + 5 from `test_ddi_api.py` + 9 from
   `test_package_smoke.py`), with `test_microbiolink_api.py` excluded per
   deviation 5.

9. **The merge commit was created, not left staged.** Step 14 said "do not
   commit... leave the result for the user to review." Once all conflicts
   were resolved and tests passed, the user was asked explicitly whether to
   conclude the merge with a commit or leave it staged, and chose to
   conclude it. Commit `532b89a` on `merge-beta-into-main`; `main` remains
   untouched at `a823dcb`; nothing was pushed to `origin`.

10. **`.gitignore`, `get_human_fasta.py`, and the DDI-only conflict
    resolutions in `__init__.py`/`cli.py`/`pyproject.toml` matched the plan
    with no deviation** — noted here for completeness, since everything
    else in this section is a deviation.

## Final outcome

- Branch: `merge-beta-into-main`, commit `532b89a`.
- `main` untouched at `a823dcb`. Nothing pushed to `origin`.
- 54/54 tests pass (`.venv/bin/python -m pytest tests/ -q`).
- Package imports cleanly; `microbiolink_api.__version__ == '0.0.2'`; all 39
  `__all__` exports resolve with no missing attributes.
- `ddi`, `dmi`, and `reverse-dmi` CLI entry points all respond correctly to
  `--help` with distinct, non-overlapping argument sets.
- `microbiolink/DMI.py`, `microbiolink/reverse_DMI.py`, and
  `tests/test_reverse_DMI.py` are byte-identical to `main`; the only change
  anywhere near reverse-DMI is the one-line alias in
  `microbiolink_api/dmi.py` (deviation 4).
