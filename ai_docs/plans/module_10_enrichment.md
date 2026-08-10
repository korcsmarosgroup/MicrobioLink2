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
2. Extracts the target gene set, positionally, per analysis level:
   - `read_target_gene_list_HMI` — tab-split, columns `[0]` and `[5]` (human protein + bacterial
     protein), unioned and de-duplicated.
   - `read_target_gene_list_TieDIE` — tab-split, columns `[0]` and `[2]` (the network's two node
     columns), unioned and de-duplicated.
3. Translates the target UniProt accessions → human gene symbols via MyGene (`scopes='uniprot'`,
   `fields='symbol'`, `species='human'`).
4. Runs `gget.enrichr(target_symbols, database=..., background_list=...)`, filters to `adj_p_val < 0.05`,
   writes the CSV, and draws a bar plot (top-20) ranked either by `combined_score` (manual twin-axis plot
   with `-log10(adj_p)` overlaid) or by `adj_p` (gget's built-in `plot=True`).

The refactored Module 10 is a **faithful port** of that transformation logic, re-pointed at (a) the
refactored table schemas (read by **column name**, not tab position) and (b) the shared
`id_resolution.translate_uniprot_to_gene_symbols` helper. `case_study_output/output/Enrichment_analysis/`
(`gget_enrichr_results_reactome_usecase_hmi.csv` / `.png` and `..._tiedie.csv` / `.png`) is the regression
fixture (per Q11).

The Enrichr results CSV schema is the contract:
`rank, path_name, p_val, z_score, combined_score, overlapping_genes, adj_p_val, database`.

## Input schema adaptation — the key gap

The ground truth reads its inputs **positionally** by tab index, which is invalid against the refactored
CSV tables (different column order, comma-separated). Module 10 reads both inputs **by column name**:

- **HMI level** — the DMI table's `human_uniprot_id` column is the human target of the microbial protein
  in **both** DMI directions (Module 9 established this). The refactored extractor takes the **distinct
  `human_uniprot_id` values only**. This is a deliberate, cleaner divergence from the ground truth, which
  also folded in the bacterial column (`line[5]`) — but those bacterial accessions cannot map to a human
  gene symbol (`species='human'`) and drop out in translation, so the **effective** target set is
  identical. Reading `human_uniprot_id` by name makes that intent explicit and drops the dead bacterial
  round-trip. *(A DDI table (Module 5) exposes the same `human_uniprot_id` column, so the HMI level
  accepts it unchanged — flagged as available, not a separate code path.)*
- **TieDIE level** — the Module 9 final network's columns are `Source.node`, `Relationship`,
  `Target.node`, `layer`. The extractor takes the **union of `Source.node` and `Target.node`**, matching
  the ground truth's "every node on the TieDie network." Bacterial binding-protein source nodes are
  included (as the ground truth does) and naturally drop out in the human-only translation, so enrichment
  runs over the human nodes of the network.

Both extractors return a list of distinct UniProt accessions, which then go through the **single** shared
translator.

## Shared code and reuse

- **`microbiolink/utils/id_resolution.py`** — `translate_uniprot_to_gene_symbols` **already exists**
  (added for Module 9, MyGene `scopes='uniprot'`, `fields='symbol'`, `species='human'`, batched, omits
  accessions with no symbol). Module 10 reuses it verbatim and the ground truth's local
  `translate_uniprot_to_symbols` is **not** ported (Q5: single implementation). Because the shared helper
  already **omits** unmapped accessions, the refactored target-symbol list is clean — no `None` values to
  filter, unlike the ground truth's `list(target_genesymbols.values())` which carried `None`s through.
- **`microbiolink/utils/gene_matrix.py`** / the CLI's non-zero-expression reader — the background gene
  universe is the same expressed-gene set used elsewhere (the Module 1 z-score / transcriptomics output);
  the background reader mirrors the existing first-column symbol read.
- **`gget` + `matplotlib`** — new optional dependencies, gated behind an `enrichment` extra, exactly as
  `idr` gates Modules 7/8 and `tiedie` gates Module 9. Lazy-imported inside function bodies (coding-style
  rule for heavy deps, and consistent with `omnipath`/`mygene`/`iupred` usage across the package).

## Decisions

The decisions below follow the refactoring precedence rules and the ground-truth logic. The ones marked
**(confirm)** are flagged for the user before/while drafting the implementation, in the same spirit as the
Module 9 plan's flagged decisions; the rest are direct ports.

1. **New `enrichment` optional extra: `gget` + `matplotlib`.** `gget.enrichr` is the enrichment engine
   (a thin client over the live [Enrichr](https://maayanlab.cloud/Enrichr) API); `matplotlib` draws the
   plot. Neither is a base dependency and both are only needed for this module, so they gate behind an
   `enrichment` extra with a `microbiolink-enrichment` console script — mirroring `idr`/`tiedie`. `mygene`
   is already a base dep (reused via `id_resolution`).
2. **Enrichment is a live Enrichr call, against a *named, versioned* gene-set library.** Default database
   `Reactome_2022` (the ground-truth default), with the same shortcut aliases exposed in help
   (`pathway`, `ontology`, `transcription`, …). Unlike Module 9's OmniPath fetch, the library name is
   version-pinned (`_2022`), so for a **fixed target set** the enrichment is far more reproducible. The
   remaining non-determinism is the target set itself: the HMI target set is deterministic from the
   fixture DMI table, but the **TieDIE target set inherits Module 9's live-OmniPath drift**, so only the
   HMI-level fixture is expected to match closely (see Verification). *(Packaging an offline enrichment
   library snapshot is Future Work, same rationale as Module 9's pinned-network item.)*
3. **Core returns `(results_df, figure)`; the CLI writes the CSV and PNG.** Following the Module 9
   file-I/O rule, all transformation stays pure (table in, DataFrame + `matplotlib.figure.Figure` out) and
   only the thin CLI layer touches disk (`to_csv`, `figure.savefig`). The plotting is **ported verbatim
   from the ground truth** — no re-implementation: the `combined_score` ranking builds the manual
   twin-axis bar (`combined_score` bars + `-log10(adj_p_val)` scatter overlay) exactly as
   `enrichr_id_database_ranking.py` does, and the `adj_p` ranking uses **gget's own built-in plot**
   (`gget.enrichr(..., plot=True)`), staying identical to the previous implementation *(confirmed with the
   user — reuse gget for now)*. The one accommodation to the `(df, figure)` contract: since gget draws to
   the global pyplot state, `plot_enrichment` grabs the active figure with `plt.gcf()` after the gget call
   so the CLI can `savefig` it through the same return path as the manual plot. This keeps the CLI's disk
   write uniform without altering gget's chart.
4. **`adj_p_val < 0.05` significance filter and top-20-by-rank plot are ported verbatim.** The results
   table written to disk is the filtered set; the plot shows the top 20 after sorting by the chosen
   ranking, with `combined_score == inf` rows excluded from the `combined_score` plot (ground-truth
   behaviour). The two ranking modes (`combined_score`, `adj_p`) are both kept.
5. **HMI target set = distinct `human_uniprot_id`; TieDIE target set = union of `Source.node` /
   `Target.node`.** Read by column name (see "Input schema adaptation"). Both feed the one shared
   translator; unmapped accessions are dropped by the translator.
6. **Background gene universe is user-supplied gene symbols** (`--background_gene_list` + `--sep`, first
   column), the expressed-gene set from Module 1 / the transcriptomics file — the same universe Module 9
   contextualises against. Passed straight to `gget.enrichr(background_list=...)`. No z-score/expression
   filtering happens in Module 10; the user supplies an already-expressed background, matching the ground
   truth and the case study.

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


# --- Enrichment + plot -----------------------------------------------------------------------------

def run_enrichment(
    target_symbols: list[str],
    background_symbols: list[str],
    database: str = DEFAULT_DATABASE,
    ranking: str = "combined_score",
) -> pd.DataFrame:
    """Run gget.enrichr against `database` with the background list; filter to adj_p_val < 0.05.

    Lazy `import gget`. Passes plot=(ranking == "adj_p") so gget draws its built-in adj_p plot as a
    side effect (ground-truth behaviour). Returns the filtered results table (RESULT_COLUMNS schema).
    """


def plot_enrichment(results: pd.DataFrame, ranking: str = "combined_score"):
    """Build the top-20 enrichment plot; return a matplotlib Figure (lazy `import matplotlib`).

    ranking='combined_score': the ground-truth manual twin-axis plot — horizontal bars by combined
    score (inf excluded) with a twin x-axis scatter of -log10(adj_p_val). ranking='adj_p': gget's own
    built-in plot via gget.enrichr(..., plot=True), captured with plt.gcf() (ported unchanged).
    """


# --- Orchestrator (the public entry point) ---------------------------------------------------------

def run_enrichment_analysis(
    target_table: pd.DataFrame,
    background_symbols: list[str],
    analysis_level: str,
    database: str = DEFAULT_DATABASE,
    ranking: str = "combined_score",
) -> tuple[pd.DataFrame, "matplotlib.figure.Figure"]:
    """Extract targets for the given analysis_level, translate to symbols, enrich, and plot.

    Args:
        target_table: A Module 5/6/7/8 HMI table (analysis_level='HMI') or a Module 9 final
            network (analysis_level='TieDIE').
        background_symbols: Expressed-gene-symbol universe for the enrichment background.
        analysis_level: 'HMI' (distinct human_uniprot_id) or 'TieDIE' (all network nodes).
        database: Enrichr gene-set library (default Reactome_2022).
        ranking: 'combined_score' or 'adj_p' — how the returned plot is ranked.

    Returns:
        (results_table, figure) — the significant-enrichment table and its top-20 plot.
    """
```

- `run_enrichment_analysis`: dispatch on `analysis_level` to `extract_hmi_targets` /
  `extract_tiedie_nodes`; `symbols = id_resolution.translate_uniprot_to_gene_symbols(targets)` then take
  `list(symbols.values())` (already `None`-free); call `run_enrichment` then `plot_enrichment`; return
  both. No disk access.
- **`ranking` threading (ground-truth faithfulness):** because gget's built-in `adj_p` plot is a side
  effect of the `gget.enrichr(..., plot=True)` call, `run_enrichment` takes `ranking` and passes
  `plot=(ranking == "adj_p")` — matching the ground truth, where the enrichr call and the plot are made
  together in the `adj_p` branch. `plot_enrichment` then either builds the manual `combined_score` figure
  or captures gget's already-drawn figure via `plt.gcf()`. Filtering (`adj_p_val < 0.05`) is independent
  of ranking and always applied in `run_enrichment`.

### `microbiolink/cli.py` (extend, argparse only)

- `_build_enrichment_parser()`: `--analysis_level` (required, choices `HMI`/`TieDIE`), `--target_file`
  (required — the DMI/HMI table or the Module 9 final network), `--target_sep` (default `,` for HMI,
  `\t` for the Module 9 network — or a single `--target_sep` the user sets), `--background_gene_list`
  (required), `--sep` (background separator, required), `--database` (default `Reactome_2022`, help lists
  the shortcut aliases), `--ranking` (default `combined_score`, choices `combined_score`/`adj_p`),
  `--output_file` (required — results CSV), `--output_image` (required — plot PNG).
- `enrichment()` entry point: read the target table (`pd.read_csv` with the level-appropriate separator),
  read the background symbols (first-column reader, `--sep`), call
  `enrichment_workflow.run_enrichment_analysis(...)`, write `results.to_csv(output_file, index=False)` and
  `figure.savefig(output_image, dpi=300, bbox_inches="tight", transparent=True)` (ground-truth savefig
  parameters), `return 0`.

### `pyproject.toml`

- Add `enrichment = ["gget>=0.28", "matplotlib>=3.8"]` under `[project.optional-dependencies]` (pin
  minimums; confirm exact floors against the lockfile).
- Add `microbiolink-enrichment = "microbiolink.cli:enrichment"` under `[project.scripts]`.
- Note that Module 10 requires the `enrichment` extra (`gget` calls the live Enrichr API).

## Migration checklist

### 1. Core module
- [ ] Create `microbiolink/workflow/enrichment.py` with `extract_hmi_targets`, `extract_tiedie_nodes`,
      `run_enrichment`, `plot_enrichment`, `run_enrichment_analysis`.
- [ ] Confirm `extract_hmi_targets` reads distinct `human_uniprot_id` by name (M5/M6/M7/M8 tables) and
      `extract_tiedie_nodes` unions `Source.node`/`Target.node` (Decision 5).
- [ ] Confirm `run_enrichment` reuses `id_resolution.translate_uniprot_to_gene_symbols` (Q5) and applies
      the `adj_p_val < 0.05` filter (Decision 4).
- [ ] Confirm plotting is ported verbatim: manual twin-axis for `combined_score`, gget's built-in
      `plot=True` (captured via `plt.gcf()`) for `adj_p` (Decision 3).

### 2. CLI wiring + packaging
- [ ] Add `_build_enrichment_parser` and `enrichment()` to `cli.py` (argparse only; core stays plain-Python).
- [ ] Add the `enrichment` extra and `microbiolink-enrichment` script to `pyproject.toml`.

### 3. Delete old code (per Q11)
- [ ] Delete `workflow/enrichr_id_database_ranking.py` once the port is validated (nothing in
      `microbiolink/` references it).

### 4. Verification — case-study regression + manual sign-off
- [ ] **HMI level** — run against the refactored HMI/DMI table + the case-study background list, database
      `Reactome_2022`, and diff the results table against `gget_enrichr_results_reactome_usecase_hmi.csv`.
      This target set is deterministic, so a close match is expected (allow for Enrichr library revisions).
- [ ] **TieDIE level** — run against the Module 9 final network + background, diff against
      `..._tiedie.csv`. The target set inherits Module 9's live-OmniPath drift (Decision 2), so an exact
      match is **not** expected — verify schema identity and manually confirm the top enriched terms are
      biologically coherent (new-feature / drift sign-off, per the refactoring plan).
- [ ] Visually inspect both PNGs against the fixtures and show to the user for sign-off.
- [ ] Confirm `run_enrichment_analysis` works with a Module 6, 7, and 8 table at the HMI level (all expose
      `human_uniprot_id`) and a Module 5 DDI table (same column).

## Verification

- `uv run ruff format` / `uv run ruff check` scoped to the new/edited files (`workflow/enrichment.py`,
  `cli.py`, `pyproject.toml` reviewed by eye).
- `uv run ty check`.
- The case-study regression (HMI-level table diff) plus the end-to-end manual sign-off above.

## Future Work

- **Offline enrichment snapshot** — package a pinned Enrichr/Reactome gene-set library for reproducible,
  offline runs, replacing the live `gget.enrichr` call (same rationale as Module 9's pinned-network item).
- **Multi-database / combined enrichment** — run several libraries (Reactome + GO + KEGG) in one call and
  emit a faceted plot, once single-database parity is confirmed.
- **Per-layer TieDIE enrichment** — enrich each Module 9 network layer separately (binding proteins vs.
  TFs vs. DEGs) rather than pooling every node, to separate direct-target from downstream signal.
```
