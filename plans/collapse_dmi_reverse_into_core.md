# Collapse DMI/reverse-DMI into `microbiolink/core/`

All code changes in this plan follow [CODING_STYLE.md](../CODING_STYLE.md).

**Status:** planned, not yet executed.

## Context

Candidate #1 from the 2026-07-10 architecture review of `microbiolink/`:
`microbiolink/DMI.py` and `microbiolink/reverse_DMI.py` reimplement
parsing/matching logic that already exists, more completely, in
`microbiolink_api/dmi.py` — concretely, `DMI.py`'s `extract_uniprot_id` raises
a plain `ValueError` on a malformed FASTA header while
`microbiolink_api/dmi.py`'s equivalent raises its own `InputFormatError` for
the identical check. Same concept, two implementations, already diverging.

During design discussion (grilling), the user established a broader
principle: `microbiolink_api` should be developed public-function-first, and
`microbiolink`'s CLI scripts should be thin wrappers that call those public
functions — the same process must run identically whether invoked via the
CLI or a public function call, for reproducibility. `microbiolink/DDI.py`
already proves this pattern works (it's a thin shim over
`microbiolink_api.ddi`); `DMI.py`/`reverse_DMI.py` never followed it.

This raised the further question of whether `microbiolink_api` should exist
as a separate top-level package at all. Investigation found no hard
technical reason for two top-level packages (no eager heavy-dependency
imports force separation — `microbiolink/__init__.py` is minimal), only
migration cost (an already-documented public import path, packaging/docs
wired around two package names). The user chose to fold that migration into
this plan now rather than defer it, given the package is pre-alpha
(`0.0.2`, `Development Status :: 2 - Pre-Alpha`) with no stable-release
promise yet.

A related, larger problem was also found and explicitly deferred: the
download/fetch script cluster (`download_protein_domains.py`,
`download_bacterial_proteins.py`, `get_protein_fasta.py`,
`get_bacterial_fasta.py`, `download_human_domains.py`, `get_human_fasta.py`)
has the same shallow-CLI-mixed-with-core-logic shape, and
`microbiolink_api/microbiome.py` already reaches into it today (via
`microbiolink.download_bacterial_proteins`, itself a re-export shim over
`microbiolink.download_protein_domains`). This is **out of scope** for this
plan — it needs its own dedicated grilling session (naming for the new core
module(s), reconciling at least two duplicate `read_ids` implementations and
three separate "is this gene expressed" implementations). One accepted
consequence: after Commit 1 below, `microbiolink/core/microbiome.py` will
still import from the top-level `microbiolink.download_bacterial_proteins`
— a pre-existing condition, not introduced by this plan, left for the
follow-up to resolve.

## Confirmed decisions

- **Output format:** preserve the exact legacy CLI output format for both
  commands (6-column `;`-separated CSV, `#`-prefixed header row) — no change
  to `dmi_output.tsv`/reverse-DMI output shape for existing pipelines.
- **Error behavior:** no special-casing needed. Once `DMI.py` delegates to
  `microbiolink_api.dmi` (soon `microbiolink.core.dmi`), it naturally raises
  `InputFormatError` (soon `MicrobioLinkError`/`InputFormatError`) instead of
  `ValueError` — consistent with the CLI-wraps-API principle, not a
  compatibility shim to preserve the old exception type.
- **`reverse_DMI.py`'s narrow public surface dropped:** `ReverseDomainMotifInteraction`
  dataclass, `predict_reverse_domain_motif_interactions_from_data`, and
  `filter_cleavage_motifs` are all deleted rather than kept as translating
  wrappers. `filter_cleavage_motifs` is redundant: `microbiolink_api.dmi`
  already excludes `CLV_*` motifs automatically inside its reverse-mode path
  (`_filter_reverse_motif_resources`/`_is_reverse_compatible_motif`).
- **`DMI.py`'s local parsing functions dropped:** `read_fasta_sequences`,
  `parse_elm_regex`, `parse_motif_domain`, `parse_protein_domain`,
  `create_uniprot_motif_dict`, `extract_uniprot_id` are all deleted.
  `microbiolink/motif_monte_carlo_filter.py`'s import of `extract_uniprot_id`
  is repointed to `microbiolink.core.dmi`.
- **New public function:**
  `resolve_dmi_resource_bundle_by_name(resource_set: str) -> DMIResourceBundle`
  added to `microbiolink_api/dmi.py` (destined for
  `microbiolink/core/dmi.py`), mapping `'default'`/`'elm'`/`'3did'` to
  `load_default_dmi_resource_bundle`/`load_default_elm_dmi_resource_bundle`/
  `load_default_3did_dmi_resource_bundle`. Exported from the package
  `__init__.py`, unit-tested directly.
- **CLI gains `--resource_set`:** both `DMI.py` and `reverse_DMI.py`'s
  `parse_args()` add `--resource_set {default, elm, 3did}` (default
  `'default'`), mirroring `DDI.py`'s existing `--resource_set` pattern
  exactly. `--elm_regex_file`/`--motif_domain_file` become optional
  (`required=False`, default `None`) so the packaged bundle can be used
  without supplying either file.
- **`cli.py:dmi()` fixed:** stops re-deriving its own `ArgumentParser`;
  delegates to `dmi_module.parse_args()` + `dmi_module.main(args)`, matching
  how `ddi()` and `reverse_dmi()` already delegate to their modules.
- **Package rename, clean break:** `microbiolink_api` is deleted outright and
  its contents moved into `microbiolink/core/` (see file list below) — no
  backward-compatible re-export shim, justified by pre-alpha status.
  `MicrobioLinkAPIError` renamed to `MicrobioLinkError` (not referenced by
  name anywhere outside `exceptions.py`/`__init__.py`, so zero-risk).
- **Test file renames:** `tests/test_microbiolink_api.py` →
  `tests/test_microbiolink_core.py`, `tests/test_ddi_api.py` →
  `tests/test_ddi_core.py` (imports updated to `microbiolink.core`).
- **`test_reverse_DMI.py` slimmed, not replaced wholesale:** delete the 3
  `filter_cleavage_motifs` tests and 3 `predict_reverse_domain_motif_interactions_from_data`
  tests (functions no longer exist); keep and update the 2 `main()`
  output-format tests (constructing `argparse.Namespace` with the new
  `resource_set` field, still asserting the legacy CSV shape).
- **New `tests/test_DMI.py`:** mirrors `test_reverse_DMI.py`'s remaining
  shape — `main()` output-format test(s) plus `--resource_set` variants.
  (No `test_DMI.py` exists today at all.)
- **`tests/test_package_smoke.py`:** add a `microbiolink-reverse-dmi`
  end-to-end case mirroring the existing `microbiolink-dmi` one (no smoke
  test for the reverse command exists today).
- **No golden byte-for-byte regression test** — declined; existing
  output-format assertions plus `microbiolink_api`'s (soon `core`'s) own
  data-level tests are considered sufficient.
- **Commit staging:** two commits, not one — see below. (Originally
  recommended one combined commit, matching `merge_beta_to_main.md` Stage
  2's rationale of avoiding a broken intermediate state; the user chose
  staged commits instead. This works cleanly here because `DMI.py`/
  `reverse_DMI.py` don't import from `microbiolink_api` today, so a
  pure-rename Commit 1 doesn't leave anything broken.)

## Commit 1: rename `microbiolink_api` → `microbiolink/core/`

Pure move, no behavior change to DMI/reverse-DMI yet.

<!-- markdownlint-disable MD013 -->
| From | To |
| --- | --- |
| `microbiolink_api/__init__.py` | `microbiolink/core/__init__.py` |
| `microbiolink_api/dmi.py` | `microbiolink/core/dmi.py` (+ new `resolve_dmi_resource_bundle_by_name`, added here so Commit 2 can use it immediately) |
| `microbiolink_api/ddi.py` | `microbiolink/core/ddi.py` |
| `microbiolink_api/workflows.py` | `microbiolink/core/workflows.py` |
| `microbiolink_api/exceptions.py` | `microbiolink/core/exceptions.py` (`MicrobioLinkAPIError` → `MicrobioLinkError`) |
| `microbiolink_api/microbiome.py` | `microbiolink/core/microbiome.py` |
| `microbiolink_api/expression.py` | `microbiolink/core/expression.py` |
| `microbiolink_api/resources/*` (7 `.tsv` files + `__init__.py`) | `microbiolink/core/resources/*` |
<!-- markdownlint-enable MD013 -->

Also in this commit:

- Repoint every internal cross-reference: `microbiolink/DDI.py`'s
  `from microbiolink_api.ddi import ...` → `from microbiolink.core.ddi import ...`;
  `workflows.py`'s internal imports; `microbiolink/core/microbiome.py`'s
  existing import from `microbiolink.download_bacterial_proteins` is left as
  a normal intra-package import (see "out of scope" above).
- No backward-compat shim package.
- Update `pyproject.toml` (`packages = ["microbiolink"]`, resource `include`
  path), `README.md`'s example, and the 3 `tutorials/mock_data/*.ipynb`
  notebooks' imports.
- Rename and update the two test files (`test_microbiolink_core.py`,
  `test_ddi_core.py`).
- End state: full test suite green; `DMI.py`/`reverse_DMI.py` byte-identical
  to before this commit.

## Commit 2: collapse `DMI.py`/`reverse_DMI.py` into thin shims

- `DMI.py`: drop all local parsing functions; `parse_args()` gains
  `--resource_set` and makes `--elm_regex_file`/`--motif_domain_file`
  optional; `main()` calls `microbiolink.core.dmi.predict_domain_motif_interactions(...)`
  and writes the unchanged legacy 6-column CSV.
- `reverse_DMI.py`: drop `ReverseDomainMotifInteraction`,
  `predict_reverse_domain_motif_interactions_from_data`,
  `filter_cleavage_motifs`; `main()` calls
  `microbiolink.core.dmi.predict_reverse_domain_motif_interactions(...)`
  directly and writes the unchanged legacy 6-column CSV (mapping
  `microbial_protein`→bacterial_protein, `domain`→human_domain,
  `host_protein`→human_protein). Same `--resource_set` addition.
- `motif_monte_carlo_filter.py`'s `extract_uniprot_id` import repoints to
  `microbiolink.core.dmi`.
- `cli.py:dmi()` delegates to `DMI.parse_args()`/`DMI.main()` instead of its
  own `ArgumentParser`.
- Tests: slim `test_reverse_DMI.py`, add new `test_DMI.py`, add the
  `microbiolink-reverse-dmi` smoke-test case.

## Out of scope (deferred, separate plan)

- Splitting the download/fetch script cluster
  (`download_protein_domains.py`, `download_bacterial_proteins.py`,
  `get_protein_fasta.py`, `get_bacterial_fasta.py`,
  `download_human_domains.py`, `get_human_fasta.py`) into core primitives vs
  thin CLI wrappers, and reconciling their duplicate `read_ids`/expression-check
  implementations. Candidates #4 and #5 from the 2026-07-10 architecture
  review.
- Backward-compatible `microbiolink_api` re-export shim — explicitly
  rejected (clean break instead), noted here only so a future pass doesn't
  re-propose it without reason.

## Deviations from the plan (Commit 1)

Minimal — execution matched the plan. Two details not spelled out verbatim
in the original text:

1. `resolve_dmi_resource_bundle_by_name` validates its input and raises
   `InputFormatError` for any `resource_set` value outside
   `{'default', 'elm', '3did'}`, matching the existing validation style of
   `predict_bidirectional_domain_motif_interactions_from_data`'s `mode`
   argument. The plan specified the function's existence and mapping but not
   this explicit error behavior.
2. `pyproject.toml` needed two more `microbiolink_api` references updated
   beyond the ones named in the plan: `[tool.coverage.run] source` and
   `[tool.ruff.lint.isort] known-first-party`. Both simply drop
   `"microbiolink_api"`, consistent with the single-package rename.

## Final outcome (Commit 1)

- All `microbiolink_api/*` content moved to `microbiolink/core/*` via
  `git mv`; `microbiolink_api` no longer exists anywhere in the tree.
- `MicrobioLinkAPIError` renamed to `MicrobioLinkError`.
- `microbiolink/DDI.py` repointed to `microbiolink.core.ddi`.
- `microbiolink/DMI.py` and `microbiolink/reverse_DMI.py` confirmed
  byte-identical to before this commit (`git diff` empty) — untouched, as
  planned.
- `pyproject.toml`, `README.md`, and all 3 `tutorials/mock_data/*.ipynb`
  notebooks updated (notebooks verified to still be valid JSON after the
  edit).
- `tests/test_ddi_api.py` → `tests/test_ddi_core.py`,
  `tests/test_microbiolink_api.py` → `tests/test_microbiolink_core.py`,
  imports updated in both.
- Editable install refreshed (`uv pip install -e .`); manual smoke check
  confirms `microbiolink-ddi/-dmi/-reverse-dmi --help` all still work and
  `import microbiolink_api` now correctly fails.
- Full test suite: **67/67 passed**, matching the pre-refactor baseline
  exactly — no regressions.
- Not yet committed to git as of writing this section; Commit 2 (the
  `DMI.py`/`reverse_DMI.py` thin-shim collapse) has not started.
