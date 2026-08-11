# Module 10 — Functional Enrichment Analysis — Implementation Plan

Detailed plan for the tenth and final module described in @ai_docs/plans/microbiolink_refactoring.md,
following the cross-cutting decisions in @ai_docs/plans/refactoring_questions.md and the precedent set
by @ai_docs/plans/module_1_zscore_filter.md through @ai_docs/plans/module_9_tiedie.md. Scope is exactly
Module 10: take either the whole Module 9 TieDie network **or** a host–microbe interaction (HMI) table,
reduce it to the set of human targets, run over-representation (Enrichr) functional enrichment against a
background gene universe, and return one results **table** and one enrichment **plot**.

The module answers the biological question "which pathways / ontology terms are over-represented among the
human proteins that microbial proteins engage — either the direct binding targets (HMI level) or the whole
downstream signalling subnetwork TieDie recovers (TieDIE level)?"

## Context

Modules 1–9 are implemented in `microbiolink/workflow/` (`zscore_filter.py`, `membrane_filter.py`,
`fasta_download.py`, `domain_download.py`, `ddi.py`, `dmi.py`, `idr_filter.py`, `monte_carlo.py`,
`tiedie.py`), `microbiolink/utils/` (`uniprot_client.py`, `id_resolution.py`, `fasta.py`,
`dmi_reader.py`, `gene_matrix.py`), and `microbiolink/data/`, all wired into `microbiolink/cli.py`.

Module 10 sits at the **end** of the pipeline. It has two entry points controlled by an
`analysis_level` argument (the ground truth's `--analysis_level`):

- **`HMI`** — enrich over the human targets of microbial proteins, taken directly from a Module 6/7/8
  DMI table (or the Module 5 DDI table — same pair columns). This path does not require TieDie.
- **`TieDIE`** — enrich over **every node** of the Module 9 final network (`network_output`), i.e. the
  binding proteins plus the whole recovered downstream signalling subnetwork.

Per the module-map, "the output is a table and an enrichment plot which will be the same for each input" —
so both levels converge on the **identical** output contract; only the target-set extraction differs.

## Ground-truth source

Module 10 **has** ground-truth code on the `case-study` branch, and it is the reference:
`workflow/enrichr_id_database_ranking.py`. It:

1. Reads a background gene-symbol list (first column of a user file, custom separator).
2. Extracts the target gene set, **positionally by tab index**, per analysis level:
   - `read_target_gene_list_HMI` — columns `[0]` and `[5]` (human protein + bacterial protein of the
     6-column HMI table), unioned and de-duplicated.
   - `read_target_gene_list_TieDIE` — columns `[0]` and `[2]`. **This is a latent bug** (see below):
     against the final-network schema, `[0]` is `Target.node` and `[2]` is `Relationship`, so it enriches
     over the Target column only.
3. Translates the target UniProt accessions → human gene symbols via MyGene (`scopes='uniprot'`,
   `fields='symbol'`, `species='human'`).
4. Runs `gget.enrichr(target_symbols, database=..., background_list=...)`, filters to `adj_p_val < 0.05`,
   writes the CSV, and draws a top-20 bar plot ranked either by `combined_score` (manual twin-axis plot
   with `-log10(adj_p)` overlaid) or by `adj_p` (gget's built-in `plot=True`).

The refactored Module 10 is a **port** of that transformation logic, re-pointed at (a) the refactored
table schemas (read by **column name**, not tab position), (b) the shared
`id_resolution.translate_uniprot_to_gene_symbols` helper, and (c) a self-contained object-oriented
matplotlib figure (no pyplot global state). `case_study_output/output/Enrichment_analysis/`
(`gget_enrichr_results_reactome_usecase_hmi.csv` / `.png` and `..._tiedie.csv` / `.png`) is the regression
fixture (per Q11).

The Enrichr results CSV schema is the contract:
`rank, path_name, p_val, z_score, combined_score, overlapping_genes, adj_p_val, database`.

## Input schema adaptation — the key gap

The ground truth reads its inputs **positionally** by tab index, which is invalid against the refactored
CSV tables (different column order, comma-separated). Module 10 reads both inputs **by column name**:

- **HMI level** — the DMI table's `human_uniprot_id` column is the human target of the microbial protein
  in **both** DMI directions (Module 9 established this). The refactored extractor takes the **distinct
  `human_uniprot_id` values only** (Decision 2). This is a deliberate, cleaner divergence from the ground
  truth, which also folded in the bacterial column (`line[5]`) — but those bacterial accessions **provably
  cannot** map to a human gene symbol: UniProt accessions are globally unique to one organism, and the
  shared translator hard-filters `species="human"` (`id_resolution.py:71`) and keeps only entries that
  return a symbol (`:77`), so a bacterial accession returns `notfound` and is dropped. The two approaches
  yield the **identical** target set; reading `human_uniprot_id` by name makes the intent explicit and
  skips the dead round-trip. *(A DDI table (Module 5) exposes the same `human_uniprot_id` column, so the
  HMI level accepts it unchanged — available, not a separate code path.)*
- **TieDIE level** — the Module 9 final network's columns are `Target.node`, `Source.node`,
  `Relationship`, `layer`. The extractor takes the **union of `Source.node` and `Target.node`**
  (Decision 7), matching the module-map's "every node on the TieDie network." **This fixes a latent
  ground-truth bug:** `read_target_gene_list_TieDIE` read columns `[0]` and `[2]` — i.e. `Target.node`
  (660 distinct nodes in the fixture) plus `Relationship` (`stimulates>`/`inhibits>`, which drops in
  translation) — so it enriched over the **Target column only**, silently missing the **137 Source-only
  nodes** (bacterial binding proteins, and TFs / PPI-sources that never appear as a target). The correct
  union is **797** distinct nodes. Bacterial source nodes are included and naturally drop out in the
  human-only translation, so enrichment runs over the human nodes of the whole network.

Both extractors return a list of distinct UniProt accessions, which then go through the **single** shared
translator.

## Shared code and reuse

- **`microbiolink/utils/id_resolution.py`** — `translate_uniprot_to_gene_symbols` **already exists**
  (added for Module 9; MyGene `scopes='uniprot'`, `fields='symbol'`, `species='human'`, batched, omits
  accessions with no symbol). Module 10 reuses it verbatim and the ground truth's local
  `translate_uniprot_to_symbols` is **not** ported (Q5: single implementation). Because the shared helper
  already **omits** unmapped accessions, the refactored target-symbol list is clean — no `None` values to
  filter, unlike the ground truth's `list(target_genesymbols.values())` which carried `None`s through.
  The symbol list is taken as `list(mapping.values())` **without de-duplication**, matching the ground
  truth (Enrichr de-duplicates internally).
- **background reader** — the background gene universe is a plain one-symbol-per-line human gene list (the
  case study's `enterocyte_colon_CD_expressed_genes.csv`), read as the first column with a user
  separator, mirroring the existing first-column symbol reads.
- **`gget`** — the enrichment engine, called **only for computation** (`gget.enrichr(..., plot=False)`);
  Module 10 never uses gget's built-in plotting (Decision 3). Lazy-imported inside the function body.
- **`matplotlib`** — used only to build the plot, via the **object-oriented `Figure()` API** (no
  `pyplot`, no global state, no backend selection); lazy-imported inside `plot_enrichment`.

Both `gget` and `matplotlib` are new optional dependencies, gated behind an `enrichment` extra and
lazy-imported inside function bodies — exactly as `idr` gates Modules 7/8 (`idr_filter.py:55`) and
`tiedie` gates Module 9 (`tiedie.py:257`), and consistent with the coding-style rule that heavy deps are
imported inside the functions that use them (`id_resolution.py:59`, `tiedie.py:51`). Lazy import is also
what keeps the package importable without the `enrichment` extra installed.

## Decisions

Resolved with the user during the plan grilling (see the running ledger). Ports of ground-truth
transformation logic are noted as such.

1. **New `enrichment` optional extra: `gget` + `matplotlib`.** `gget.enrichr` is the enrichment engine
   (a thin client over the live [Enrichr](https://maayanlab.cloud/Enrichr) API); `matplotlib` draws the
   plot. Neither is a base dependency and both are only needed here, so they gate behind an `enrichment`
   extra with a `microbiolink-enrichment` console script — mirroring `idr`/`tiedie`. `mygene` is already a
   base dep (reused via `id_resolution`). Both are lazy-imported.
2. **HMI target set = distinct `human_uniprot_id`.** Read by column name; bacterial IDs are provably
   dropped by the human-only translator, so this equals the ground-truth result (see Input schema
   adaptation). A Module 5 DDI table is accepted unchanged (same column).
3. **Both plots are built by Module 10 with the object-oriented `Figure()` API; gget is called with
   `plot=False`.** The core returns `(results_df, figure)` and the CLI writes both — all transformation
   stays pure and only the thin CLI layer touches disk. Using `matplotlib.figure.Figure()` directly
   avoids pyplot's global state entirely: no backend juggling, no `matplotlib.use()`/`MPLBACKEND`
   ordering, and `fig.savefig()` auto-selects the Agg file-canvas for the PNG regardless of environment.
   - `combined_score` ranking: a **faithful port** of the ground truth's twin-axis plot — horizontal bars
     by combined score (`inf` excluded) with a twin x-axis (`ax.twiny()`) scatter of `-log10(adj_p_val)`.
   - `adj_p` ranking: **our own** horizontal bar of `-log10(adj_p_val)`. This is a faithful
     re-expression, **not** byte-identical to gget's old built-in `adj_p` PNG. *(Reversed from an earlier
     "reuse gget's plot" decision: gget's plot lives in pyplot global state, which forced a fragile
     backend-ordering hack; building both plots ourselves removes it. The **results CSV is unaffected** —
     gget still computes the enrichment identically — only the `adj_p` image differs from the old output.)*
4. **`adj_p_val < 0.05` significance filter and top-20-by-rank plot are ported verbatim.** The results
   table written to disk is the full filtered set (gget's column order, `index=False`); the plot shows the
   top 20 after sorting by the chosen ranking, with `combined_score == inf` rows excluded from the
   `combined_score` plot. Both ranking modes (`combined_score`, `adj_p`) are kept.
5. **TieDIE target set = `Source.node ∪ Target.node` (every node), fixing the ground-truth Target-only
   bug** (see Input schema adaptation). Both extractors feed the one shared translator; unmapped
   accessions are dropped.
6. **Background gene universe is a user-supplied gene-symbol list** (`--background_gene_list` + `--sep`,
   first column) — the expressed-gene set from Module 1 / the transcriptomics output. Passed straight to
   `gget.enrichr(background_list=...)`. No expression filtering happens in Module 10; the user supplies an
   already-expressed background, matching the ground truth and the case study.
7. **Enrichment is a live Enrichr call against a named, versioned library.** Default database
   `Reactome_2022` (the ground-truth default), with the same shortcut aliases exposed in help
   (`pathway`, `ontology`, `transcription`, …) — gget resolves these natively. `species="human"` is
   passed explicitly. Unlike Module 9's OmniPath fetch, the library name is version-pinned, so for a
   **fixed target set** the enrichment is reproducible (see Verification). *(Packaging an offline library
   snapshot is Future Work.)*
8. **Separator auto-derived from `analysis_level`** (Q1): `HMI` ⇒ comma (Module 6/7/8 CSV), `TieDIE` ⇒
   tab (Module 9 network writes `to_csv(sep="\t")`). No `--target_sep` flag — both upstream formats are
   fixed and known, so the level fully determines the separator.
9. **Empty results are handled gracefully** (Q4/Q6): when no term passes `adj_p_val < 0.05`, the core
   returns `figure=None`; the CLI writes the (empty) CSV, logs a clear warning, **skips** the PNG (no
   placeholder), and exits 0. `head(20)` naturally handles fewer than 20 significant rows.

## Target implementation

### `microbiolink/workflow/enrichment.py` (new, core, argparse-free)

```python
"""Functional enrichment analysis: reduce a network/HMI table to human targets, run Enrichr, plot."""

import numpy as np
import pandas as pd

from ..utils import id_resolution

HMI_HUMAN_COLUMN = "human_uniprot_id"
NETWORK_SOURCE_COLUMN = "Source.node"
NETWORK_TARGET_COLUMN = "Target.node"
DEFAULT_DATABASE = "Reactome_2022"
SIGNIFICANCE_CUTOFF = 0.05
TOP_N = 20
RESULT_COLUMNS = [
    "rank", "path_name", "p_val", "z_score",
    "combined_score", "overlapping_genes", "adj_p_val", "database",
]


# --- Target extraction (one per analysis level) ----------------------------------------------------

def extract_hmi_targets(hmi_table: pd.DataFrame) -> list[str]:
    """Distinct human UniProt targets from a Module 5/6/7/8 table (the human_uniprot_id column)."""


def extract_tiedie_nodes(network: pd.DataFrame) -> list[str]:
    """Every node of a Module 9 final network: the union of Source.node and Target.node."""


# --- Enrichment ------------------------------------------------------------------------------------

def run_enrichment(
    target_symbols: list[str],
    background_symbols: list[str],
    database: str = DEFAULT_DATABASE,
) -> pd.DataFrame:
    """Run gget.enrichr against `database` with the background list; filter to adj_p_val < 0.05.

    Lazy `import gget`; called with plot=False, species="human". Returns the filtered results table
    (RESULT_COLUMNS schema); may be empty if nothing passes the cutoff.
    """


# --- Plot (object-oriented Figure API; no pyplot) --------------------------------------------------

def plot_enrichment(results: pd.DataFrame, ranking: str = "combined_score"):
    """Build the top-20 enrichment plot; return a matplotlib Figure, or None if `results` is empty.

    Lazy `from matplotlib.figure import Figure` (no pyplot / backend selection).
    ranking='combined_score': ground-truth twin-axis plot — horizontal bars by combined score (inf
    excluded) with a twin x-axis scatter of -log10(adj_p_val). ranking='adj_p': horizontal bars of
    -log10(adj_p_val).
    """


# --- Orchestrator (the public entry point) ---------------------------------------------------------

def run_enrichment_analysis(
    target_table: pd.DataFrame,
    background_symbols: list[str],
    analysis_level: str,
    database: str = DEFAULT_DATABASE,
    ranking: str = "combined_score",
) -> tuple[pd.DataFrame, "matplotlib.figure.Figure | None"]:
    """Extract targets for the given analysis_level, translate to symbols, enrich, and plot.

    Args:
        target_table: A Module 5/6/7/8 HMI table (analysis_level='HMI') or a Module 9 final
            network (analysis_level='TieDIE').
        background_symbols: Expressed-gene-symbol universe for the enrichment background.
        analysis_level: 'HMI' (distinct human_uniprot_id) or 'TieDIE' (all network nodes).
        database: Enrichr gene-set library (default Reactome_2022).
        ranking: 'combined_score' or 'adj_p' — how the returned plot is ranked.

    Returns:
        (results_table, figure) — the significant-enrichment table and its top-20 plot; figure is
        None when no term passes the significance cutoff.
    """
```

- `run_enrichment_analysis`: dispatch on `analysis_level` to `extract_hmi_targets` /
  `extract_tiedie_nodes`; `symbols = list(id_resolution.translate_uniprot_to_gene_symbols(targets).values())`
  (already `None`-free, undeduplicated per the ground truth); call `run_enrichment` then
  `plot_enrichment`; return both. No disk access. `run_enrichment` (data) and `plot_enrichment` (figure)
  are fully separated — `ranking` reaches only `plot_enrichment`.

### `microbiolink/cli.py` (extend, argparse only)

- `_build_enrichment_parser()`: `--analysis_level` (required, choices `HMI`/`TieDIE`), `--target_file`
  (required — the DMI/HMI table or the Module 9 final network), `--background_gene_list` (required),
  `--sep` (background separator, required), `--database` (default `Reactome_2022`, help lists the
  shortcut aliases), `--ranking` (default `combined_score`, choices `combined_score`/`adj_p`),
  `--output_file` (required — results CSV), `--output_image` (required — plot PNG). **No `--target_sep`**
  (Decision 8).
- `enrichment()` entry point: read the target table with the level-appropriate separator
  (`pd.read_csv(sep="," if HMI else "\t")`), read the background symbols (first-column reader, `--sep`),
  `from .workflow import enrichment as enrichment_workflow` (lazy, per the cli precedent), call
  `run_enrichment_analysis(...)`, write `results.to_csv(output_file, index=False)`; if the figure is not
  None, `figure.savefig(output_image, dpi=300, bbox_inches="tight", transparent=True)` (ground-truth
  savefig parameters), else log a warning that no significant terms were found and no plot was written;
  `return 0`.

### `pyproject.toml`

- Add `enrichment = ["gget>=<floor>", "matplotlib>=<floor>"]` under `[project.optional-dependencies]`
  (pin minimums against the lockfile once installed; `gget`'s `background_list` support sets the floor).
- Add `microbiolink-enrichment = "microbiolink.cli:enrichment"` under `[project.scripts]`.
- Note that Module 10 requires the `enrichment` extra (`gget` calls the live Enrichr API at runtime).

## Migration checklist

### 1. Core module
- [x] Create `microbiolink/workflow/enrichment.py` with `extract_hmi_targets`, `extract_tiedie_nodes`,
      `run_enrichment`, `plot_enrichment`, `run_enrichment_analysis`.
- [x] Confirm `extract_hmi_targets` reads distinct `human_uniprot_id` by name (M5/M6/M7/M8 tables)
      (Decision 2) and `extract_tiedie_nodes` unions `Source.node`/`Target.node` = 797 for the fixture,
      fixing the Target-only bug (Decision 5). *(Verified: 797 = 660 Target + 137 Source-only.)*
- [x] Confirm `run_enrichment` reuses `id_resolution.translate_uniprot_to_gene_symbols` (Q5), calls
      `gget.enrichr(plot=False, species="human")`, and applies `adj_p_val < 0.05` (Decision 4).
- [x] Confirm `plot_enrichment` uses the OO `Figure()` API (no pyplot), builds both ranking modes, and
      returns `None` on empty input (Decisions 3, 9). *(Note: `_plot_combined_score` now coerces
      `combined_score` numerically — live gget returns it as an object column mixing floats and the
      string "inf", which crashed sorting; the CSV fixture had masked this.)*

### 2. CLI wiring + packaging
- [x] Add `_build_enrichment_parser` and `enrichment()` to `cli.py` (argparse only; core stays plain
      Python), with separator auto-derived from `--analysis_level` and the empty-figure skip+warn path
      (Decisions 8, 9). *(Parser split into `_add_enrichment_source_arguments` to stay under the
      50-line function limit.)*
- [x] Add the `enrichment` extra (`gget>=0.29.0`, `matplotlib>=3.8.0`) and `microbiolink-enrichment`
      script to `pyproject.toml`; `uv.lock` updated.

### 3. Delete old code (per Q11)
- [x] Delete `workflow/enrichr_id_database_ranking.py` once the port is validated (nothing in
      `microbiolink/` references it).

### 4. Verification — case-study regression + manual sign-off
- [~] **HMI level (exact-ish, Q5)** — **could not run against the true fixture input.** The empty
      `gget_enrichr_results_reactome_usecase_hmi.csv` was built from the *final IDR/MC-filtered* HMI
      table, which is not present in the repo as a fixture. The HMI path was instead validated against
      the *raw* prediction output (513 distinct targets), which legitimately yields ~154 terms — more
      targets, more enrichment, so **not** a target-extraction/translation regression. To close this
      exactly, re-run with `--target_file` pointed at the refactored final-filtered HMI table.
- [x] **TieDIE level (Q7)** — ran against the frozen `usecase_final_network.txt`; `extract_tiedie_nodes`
      returned the **797**-node union. Enrichment CSV: **schema identical** to the fixture; **660-term
      overlap** (Jaccard 0.888), same top-10 term set, 18/20 top-20 overlap. The 30 gained / 53 dropped
      terms are the expected divergence from the buggy 660-node Target-only fixture.
- [~] Confirm `run_enrichment_analysis` works with a Module 6, 7, and 8 table at the HMI level and a
      Module 5 DDI table. *(Validated on one reconstructed `human_uniprot_id` table; all four formats
      share that single column, so the code path is identical, but each was not run separately.)*
- [x] Visually inspect the produced PNGs — both `combined_score` (twin-axis port) and `adj_p`
      (re-expression) render correctly. **User signed off on the figures (2026-08-11).**

## Verification

- `uv run ruff format` / `uv run ruff check` scoped to the new/edited files (`workflow/enrichment.py`,
  `cli.py`, `pyproject.toml` reviewed by eye).
- `uv run ty check`.
- The case-study regression (HMI exact-ish + TieDIE node-set/overlap) plus the manual sign-off above.

## Future Work

- **Update the default library from `Reactome_2022`** — bump the default database to the latest Enrichr
  Reactome release once one is available, so enrichment reflects current pathway annotations.
- **Offline enrichment snapshot** — package a pinned Enrichr/Reactome gene-set library for reproducible,
  offline runs, replacing the live `gget.enrichr` call (same rationale as Module 9's pinned-network item).
- **Multi-database / combined enrichment** — run several libraries (Reactome + GO + KEGG) in one call and
  emit a faceted plot, once single-database parity is confirmed.
- **Per-layer TieDIE enrichment** — enrich each Module 9 network layer separately (binding proteins vs.
  TFs vs. DEGs) rather than pooling every node, to separate direct-target from downstream signal.
```
