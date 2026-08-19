# Domain Docs

How the engineering skills should consume this repo's domain documentation when exploring the codebase.

## Before exploring, read these

- **`CONTEXT.md`** at the repo root, or
- **`CONTEXT-MAP.md`** at the repo root if it exists — it points at one `CONTEXT.md` per context. Read each one relevant to the topic.
- **`ai_docs/decisions/`** — read the decision records (ADRs) that touch the area you're about to work in.

> Note: this repo keeps internal agent/dev docs under `ai_docs/` (`ai_docs/plans/`, `ai_docs/decisions/`, `ai_docs/agents/`). `docs/` is reserved for the published Zensical documentation site — never put internal docs there.

If any of these files don't exist, **proceed silently**. Don't flag their absence; don't suggest creating them upfront. The `/domain-modeling` skill (reached via `/grill-with-docs` and `/improve-codebase-architecture`) creates them lazily when terms or decisions actually get resolved.

## File structure

Single-context repo (this repo):

```
/
├── CONTEXT.md                         ← domain glossary (repo root)
├── ai_docs/
│   ├── decisions/                     ← decision records (ADRs)
│   │   ├── monte_carlo_shuffling_strategy.md
│   │   └── tiedie_inputs.md
│   ├── plans/                         ← planning docs
│   └── agents/                        ← agent/skill config (this file)
├── docs/                              ← published Zensical site (NOT internal docs)
└── microbiolink/
```

## Use the glossary's vocabulary

When your output names a domain concept (in an issue title, a refactor proposal, a hypothesis, a test name), use the term as defined in `CONTEXT.md`. Don't drift to synonyms the glossary explicitly avoids.

If the concept you need isn't in the glossary yet, that's a signal — either you're inventing language the project doesn't use (reconsider) or there's a real gap (note it for `/domain-modeling`).

## Flag ADR conflicts

If your output contradicts an existing ADR, surface it explicitly rather than silently overriding:

> _Contradicts ADR-0007 (event-sourced orders) — but worth reopening because…_
