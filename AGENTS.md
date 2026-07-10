# Code Style

Follow PEP8 and Google style guide with the following additional rules:

## Docstrings

Napoleon (Google) style, opening and closing quotes on their own lines. Use
type hints instead of types in docstrings.

## Markdown standards

Always run markdownlint on any markdown files created or edited.
Install using: `pixi global install markdownlint-cli`
Fix all linting issues before completing the task.

## Mandatory Testing

Every component of the data pipeline (loaders, transformers, clustering
wrappers, evaluators) MUST have a corresponding test file in the `tests/`
directory. No pipeline code should be considered "complete" or "validated"
until its tests pass.

## Testing preferences

Write all Python tests as pytest style functions, not unittest classes.
Use descriptive function names starting with `test_`.
Prefer fixtures over setup/teardown methods.
Use assert statements directly, not `self.assertEqual`.

## Testing approach

Never create throwaway test scripts or ad hoc verification files.
If you need to test functionality, write a proper test in the test suite.
All tests go in the `tests/` directory following the project structure.
Tests should be runnable with the rest of the suite (`pixi run pytest`).
Even for quick verification, write it as a real test that provides ongoing
value.

## Package management

This project uses Pixi for all package management.
Never run commands directly (`python`, `pytest`, etc.).
Always prefix commands with `pixi run <command>`.
Example: `pixi run python script.py` not `python script.py`.
Example: `pixi run pytest` not `pytest`.

## Reproducibility & Determinism

All scripts using randomized algorithms (e.g., PCA, t-SNE, K-Means) must set
a fixed `random_state` or seed to ensure reproducible results.

## Data Integrity

The `peaklist_inputs/` directory must be treated as read-only. All
intermediate data (binned matrices, normalized features) and final results
(plots, metrics) must be saved to a `results/` or `output/` directory.

## Path Management

Always use project-root relative paths (utilizing `pathlib.Path`) for data
loading to ensure scripts remain portable and runnable via `uv run` from any
directory. **Use `pathlib.Path` for ALL filesystem-related code in Python
(no `os.path` or raw strings for paths).**

## Dependency Management

Before implementing new phases, verify that all necessary dependencies
(e.g., `scikit-learn`, `pandas`, `seaborn`) are explicitly defined in
`pixi.toml` and synchronized.

## Planning non-trivial changes

Before implementing a non-trivial multi-step change or refactor, write a
plan document to `plans/<short-name>.md` describing what will change and
why, and get it confirmed before writing any code — follow the style of
`plans/merge_beta_to_main.md` (status line, context, confirmed decisions,
steps/commits, files touched).
After execution, add a "Deviations from the plan" section to the same file
documenting anything that differed from what was written, plus a final
outcome summary (test results, files touched).
