# Plan: Building the MicrobioLink Python Package

Status: proposed
Date: 2026-08-17

## Goal

Produce a correct, installable distribution (sdist + wheel) for `microbiolink`, decide
whether to keep `hatchling` or move to the `uv_build` backend, and document a repeatable
build + verify workflow. Publishing to PyPI is treated as optional future work.

## TL;DR

1. **The current build is broken** — the wheel ships the packaged `.tsv` data but **zero
   Python source files**, so `import microbiolink` fails after install. This is caused by a
   single misconfigured line in `pyproject.toml` and must be fixed first.
2. **Keep `hatchling`.** For this project the `uv_build` backend offers no meaningful
   benefit and adds cost (uv upgrade, layout/config changes). See the comparison below.
3. After the fix, `uv build` already produces a valid sdist + wheel; the rest of this plan
   is verification, metadata hygiene, and an optional publishing path.

---

## 1. Critical finding: the current wheel contains no code

Building today with `uv build` succeeds *silently* but produces a broken wheel.

Evidence (verified 2026-08-17):

```
$ uv build            # succeeds
$ unzip -l dist/microbiolink-2.1.0-py3-none-any.whl
  microbiolink/data/3did_dmi_classes.tsv
  microbiolink/data/... (6 .tsv files only)
  microbiolink-2.1.0.dist-info/...
# .py files in wheel: 0

$ pip install microbiolink-2.1.0-py3-none-any.whl   # installs fine
$ python -c "import microbiolink.workflow.ddi"
ModuleNotFoundError: No module named 'microbiolink.workflow'
```

### Root cause

In `pyproject.toml`:

```toml
[tool.hatch.build]
include = ["microbiolink/data/*.tsv"]
```

A top-level `[tool.hatch.build] include` **replaces** hatchling's default file selection for
*every* build target. Once it is present, hatchling includes *only* files matching that
pattern — i.e. just the `.tsv` data — and drops all `.py` source. The
`[tool.hatch.build.targets.wheel] packages = ["microbiolink"]` setting cannot recover the
files because the global `include` has already filtered them out.

### The fix (verified)

Remove the narrowing top-level `include`. The `.tsv` files live inside the `microbiolink`
package and are **not** git-ignored, so hatchling ships them automatically once the package
itself is included. Delete these three lines:

```toml
[tool.hatch.build]
include = ["microbiolink/data/*.tsv"]
```

Verified result after removal:

```
.py files in wheel:   30
.tsv files in wheel:   6
__pycache__ in wheel:  0
import microbiolink.workflow.ddi, microbiolink.cli.ddi  ->  IMPORT OK
```

If we later want to be *explicit* about data inclusion (defensive against future changes to
`.gitignore`), the correct hatchling idiom is `artifacts`/`force-include` scoped to the wheel
target, not a global `include`:

```toml
[tool.hatch.build.targets.wheel]
packages = ["microbiolink"]
force-include = { "microbiolink/data" = "microbiolink/data" }
```

But the minimal fix (just deleting the bad block) is sufficient and preferred.

> Note: `microbiolink/workflow/` is still present and required at runtime (e.g.
> `cli/ddi.py` does `from ..workflow import ddi`). The commit "Deleted the workflow package"
> referred to the *old top-level* Snakemake `workflow/` directory, not this package
> sub-module. It must remain in the distribution.

---

## 2. Build backend decision: `hatchling` vs `uv_build`

| Factor | hatchling (current) | uv_build |
| --- | --- | --- |
| Maturity | Mature, PyPA-adopted, ubiquitous | Stable only since uv 0.7.19 (mid-2025) |
| Installed tooling | Works with current uv 0.5.26 | Requires upgrading uv (≥0.7.19) |
| Build speed | Fast enough (~1s here) | Faster (Rust), matters for large/frequent rebuilds — not our case |
| Project layout | Flat layout (`microbiolink/` at repo root) works as-is | Prefers `src/` layout; flat layout needs explicit `module-name` config |
| Data files (`*.tsv`) | Included via package (once bug fixed); `force-include` escape hatch available | Package data included by default; `[tool.uv.build-backend]` config is newer/less battle-tested |
| Direct-reference deps (git `iupred`/`tiedie`) | Handled via `allow-direct-references` | Same PyPI restriction applies regardless of backend |
| Config surface | Well-documented, lots of examples | Minimal config, fewer escape hatches, smaller doc corpus |

### Recommendation: **stay on `hatchling`.**

Reasoning:

- `uv build` (the *frontend*) is already what we use and is backend-agnostic — we keep the
  fast uv-driven build experience regardless of backend choice. The `uv_build` *backend* only
  changes the build engine, and its headline advantage (raw build speed) is irrelevant for a
  package built occasionally.
- Switching would require upgrading uv (0.5.26 → ≥0.7.19) and likely restructuring to a
  `src/` layout or adding `[tool.uv.build-backend]` module config — real work for no user-facing
  benefit.
- hatchling is more flexible for the escape hatches we may want (`force-include`,
  version hooks, richer file selection) and has broader documentation.

Keep `uv_build` on the table only if we later (a) move to a `src/` layout during the ongoing
refactor, and (b) want the whole toolchain to be uv-native. It is a low-stakes reversible
choice — revisit if those conditions arise.

---

## 3. Step-by-step build & verify plan

### Step 1 — Fix packaging config
- Edit `pyproject.toml`: delete the `[tool.hatch.build] include = ["microbiolink/data/*.tsv"]`
  block (Section 1).

### Step 2 — Metadata hygiene (quick review while in the file)
- `requires-python = ">=3.9"` — confirm this is still accurate given current deps
  (`pandas>=2.2`, `scipy>=1.13` support 3.9, but the dev interpreter here is 3.12). Keep 3.9
  only if we actually test it; otherwise raise the floor.
- `version = "2.1.0"` — confirm this is the intended release version for this build.
- Duplicate/legacy files: repo has both `claude.md` (33 bytes) and `CLAUDE.md`; not a
  packaging blocker but worth a cleanup pass. Neither ships in the wheel.
- README installation section is stale (references the old Snakemake `workflow/` scripts and
  conda env). Not required for a build, but should be updated before any public release —
  tracked separately under the documentation plan.

### Step 3 — Clean build
```bash
rm -rf dist/
uv build
```
Produces `dist/microbiolink-2.1.0.tar.gz` and `dist/microbiolink-2.1.0-py3-none-any.whl`.

### Step 4 — Inspect artifacts
```bash
unzip -l dist/*.whl        # expect: 30 .py + 6 .tsv, no __pycache__
tar tzf dist/*.tar.gz      # sdist: package sources + data + pyproject + README + LICENSE
```
Checklist:
- [ ] Wheel contains all sub-packages: `microbiolink/`, `cli/`, `utils/`, `workflow/`, `data/`.
- [ ] All 6 `.tsv` files present under `microbiolink/data/`.
- [ ] No `__pycache__` / `.pyc` leaked.
- [ ] sdist does **not** balloon with the root PDF/JPG/logo (~1.5 MB). Confirm sdist stays
      small; if the images sneak in, exclude them from the sdist target.
- [ ] `entry_points.txt` lists all 10 `microbiolink-*` console scripts.

### Step 5 — Install-and-import smoke test (clean venv)
```bash
python3 -m venv /tmp/mbl_verify && source /tmp/mbl_verify/bin/activate
pip install dist/*.whl
python -c "import microbiolink.cli.ddi, microbiolink.workflow.ddi; print('import ok')"
microbiolink-ddi --help          # console script resolves
python -c "import importlib.resources as r, microbiolink.data as d; \
  print(list(p.name for p in r.files(d).iterdir() if p.name.endswith('.tsv')))"
deactivate
```
Expect: import ok, `--help` prints usage, all 6 `.tsv` resolvable via `importlib.resources`
(this is how `workflow/ddi.py` loads them — the packaging bug would have broken it in the
field).

### Step 6 — Validate metadata
```bash
uvx twine check dist/*        # or: pipx run twine check dist/*
```
Confirms long-description/metadata render correctly.

### Step 7 — Test the optional extras resolve (best-effort)
The `enrichment`, `idr`, and `tiedie` extras pull heavy / git-based deps. At minimum dry-run:
```bash
pip install "dist/microbiolink-2.1.0-py3-none-any.whl[enrichment]"   # PyPI deps only
```
Note the `idr` and `tiedie` extras use **git direct references** — see Section 4.

---

## 4. Publishing (optional / future work)

Not required to "build the package", but flagged because it constrains choices:

- **Direct-reference dependencies block PyPI upload.** The `idr` and `tiedie` extras use
  `... @ git+https://...` URLs. PyPI **rejects** packages whose metadata contains direct
  URL references. So the wheel as-is can be built and installed locally / from a git checkout,
  but **cannot be uploaded to PyPI** while those extras carry git URLs.
  - Options if we ever publish: (a) get `iupred`/`tiedie` onto PyPI and pin normal version
    specifiers; (b) drop them from packaged extras and document manual install in the docs;
    (c) publish only to a private index / distribute the wheel directly.
- If/when publishable: dry-run to **TestPyPI** first
  (`uvx twine upload -r testpypi dist/*`), verify a clean install from TestPyPI, then upload
  to real PyPI. Use a scoped API token, not a password.
- This intersects the "Get Started / installation" documentation page — align the documented
  install method with whatever distribution channel we choose.

---

## 5. Open questions for the user

1. **Target Python versions** — keep `>=3.9`, or raise the floor to match what we actually
   test (e.g. 3.10/3.11/3.12)? Affects the `requires-python` and any CI matrix.
2. **Is publishing to PyPI a goal?** If yes, we must resolve the git-dependency constraint in
   Section 4 before release. If no, the plan ends at Step 6 (local/GitHub-installable wheel).
3. **Data-file inclusion style** — minimal fix (delete bad block, rely on default inclusion)
   vs. explicit `force-include`. Recommend minimal now; revisit if `.gitignore` changes.
4. **Backend** — confirm we keep `hatchling` (recommended). Only reconsider `uv_build` if the
   refactor moves us to a `src/` layout.

## Acceptance criteria

- `uv build` produces sdist + wheel containing all Python sources **and** the 6 `.tsv` data
  files, with no `__pycache__`.
- Fresh-venv `pip install dist/*.whl` succeeds; `import microbiolink...` works; all 10
  `microbiolink-*` console scripts are registered and `--help` runs; packaged `.tsv` files are
  resolvable via `importlib.resources`.
- `twine check dist/*` passes.
- `pyproject.toml` build config reviewed and `hatchling` confirmed as the backend.
