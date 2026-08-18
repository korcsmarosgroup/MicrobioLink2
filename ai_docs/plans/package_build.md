# Plan: Building the MicrobioLink Python Package

Status: proposed
Date: 2026-08-17

## Decisions (locked)

- **Python floor: `>=3.12`.** No dependency requires an older version, so we drop 3.9–3.11
  (3.12 will soon be the oldest supported line).
- **Backend: `hatchling`** (unchanged).
- **Data inclusion: minimal** — delete the bad `include` block and rely on default package
  inclusion; no explicit `force-include`.
- **PyPI publishing: future goal only** — captured in the Future Work section, not part of
  this build.

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
- Set `requires-python = ">=3.12"` (was `>=3.9`). All runtime deps (`pandas>=2.2`,
  `scipy>=1.13`, `numpy>=1.26`, `mygene`, `omnipath`, `requests`) and the git extras support
  3.12, so nothing forces a lower floor.
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

## 4. Distribution scope for now

This build targets a **locally / GitHub-installable** wheel + sdist. That means:

- `uv build` output can be installed directly (`pip install <wheel>`) or via
  `pip install git+https://github.com/.../MicrobioLink2.git` / `uv add` from the repo.
- The git-based extras (`idr`, `tiedie`) work in this mode because direct URL references are
  allowed for local/VCS installs (`allow-direct-references = true` is already set).

The plan ends at Step 6 (verified artifacts). No PyPI upload.

---

## 5. Future Work: publishing to PyPI

Deferred, but recorded here because it constrains dependency choices when we get to it:

- **Direct-reference dependencies block PyPI upload.** The `idr` and `tiedie` extras use
  `... @ git+https://...` URLs. PyPI **rejects** any package whose metadata contains direct
  URL references, so the current `pyproject.toml` cannot be uploaded as-is.
  - Availability checked (2026-08-17): `tiedie` **is** on PyPI (but our pin is the
    `saezlab/tiedie` git fork — confirm the PyPI release is equivalent before switching);
    `iupred`/`aiupred` are **not** on PyPI (404).
  - Resolution options: (a) switch `tiedie` to a normal PyPI version specifier if the fork's
    changes are upstreamed/unneeded; (b) for `iupred`, either get it published or drop it from
    packaged extras and document the manual `git+` install in the docs; (c) publish only to a
    private index / distribute the wheel via GitHub Releases.
### Upload tooling: `uv publish` is the primary path

The upload itself is uv-native — we do **not** need twine to publish:

- **Primary: `uv publish`.** uv both builds and uploads, and handles auth the same ways twine
  does. Preferred order:
  1. **Trusted Publishing (OIDC) from GitHub Actions** — no stored secret at all; the workflow
     mints a short-lived credential. Configure the PyPI project to trust the repo/workflow, then
     `uv publish --trusted-publishing always`.
  2. Fallback for a manual/local upload: a scoped PyPI **API token** via `UV_PUBLISH_TOKEN`
     (or `uv publish --token ...`). Never a password.
- **`twine check` — optional pre-flight only.** uv has no equivalent, so keep this one twine
  command to validate that package metadata and the README long-description render correctly on
  the PyPI page: `uvx twine check dist/*` (runs via uvx, nothing installed into the project).
  twine is **not** needed for the upload.

### Release sequence (when ready)

1. `uvx twine check dist/*` — metadata/README render check.
2. Dry-run to **TestPyPI**: `uv publish --publish-url https://test.pypi.org/legacy/ ...`.
3. Verify a clean install from TestPyPI in a fresh venv.
4. Publish to real PyPI via `uv publish` (Trusted Publishing preferred).
5. Align the documented install method on the "Get Started" docs page with the chosen channel.

---

## Acceptance criteria

- `uv build` produces sdist + wheel containing all Python sources **and** the 6 `.tsv` data
  files, with no `__pycache__`.
- Fresh-venv `pip install dist/*.whl` succeeds; `import microbiolink...` works; all 10
  `microbiolink-*` console scripts are registered and `--help` runs; packaged `.tsv` files are
  resolvable via `importlib.resources`.
- `twine check dist/*` passes.
- `pyproject.toml` build config reviewed: `hatchling` backend confirmed, minimal data
  inclusion, and `requires-python = ">=3.12"`.
