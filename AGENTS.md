# AGENTS.md

## Project context
MicrobioLink predicts host–microbe protein–protein interactions and analyzes downstream effects on host signaling pathways using multi-omic data and network biology. The current goal is to convert existing scripts into a modern Python package set managed with UV.

## Current roadmap (high level)
- Split selected scripts into independent packages.
- Use a shared project template to scaffold each new package.
- Adopt UV for tooling (build, test, lint, and publishing workflows).

## Conventions
- Prefer small, focused packages with clear boundaries and minimal cross-dependencies.
- Keep APIs stable and document public entry points.
- Avoid introducing heavyweight dependencies unless justified by a pipeline step.

## Workflow guidance
- Inspect existing scripts to identify logical package boundaries.
- For each package, define: purpose, inputs/outputs, and CLI/API surface.
- Scaffold new packages from the chosen template, then migrate code incrementally.
- Add tests around critical data transformations and network logic.

## Notes for new contributors
- Start by reading `README.md` for the current pipeline layout.
- Keep changes incremental; avoid large refactors until packaging structure is in place.

