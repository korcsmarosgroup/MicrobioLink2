# Module 9 — TieDie Network Propagation — Implementation Plan

Detailed plan for the ninth module described in @ai_docs/plans/microbiolink_refactoring.md, following the
cross-cutting decisions in @ai_docs/plans/refactoring_questions.md and the precedent set by
@ai_docs/plans/module_1_zscore_filter.md through @ai_docs/plans/module_8_monte_carlo.md. Scope is exactly
Module 9: take a MicrobioLink DMI table (optionally augmented with a Module 5 DDI table) plus a
user-supplied differentially-expressed-gene (DEG) list, build the three TieDie input files, run the
TieDie tied-diffusion algorithm, and reassemble its output into a single human-readable host–microbe
signalling network and node-annotation table.

TieDie ("Tied Diffusion for Subnetwork Discovery") connects an **upstream** node set (here: the human
proteins targeted by microbial proteins, predicted by DMI and/or DDI) to a **downstream** node set
(here: the transcription factors
inferred to drive the observed DEGs) by bidirectional heat diffusion over a human signalling network,
returning the linking subnetwork. The algorithm lives in the saez lab package
`github.com/saezlab/tiedie`.

The two heat-input modelling choices (upstream bacterial-target heat and downstream TF-activity heat,
i.e. Decisions 5, 6, and 11 below) are documented in depth — with pros/cons and better alternatives — in
@ai_docs/decisions/tiedie_inputs.md. Read that decision document for the *why* behind the heats, this
plan for the *how*.

## Context

Modules 1–8 are implemented in `microbiolink/workflow/` (`zscore_filter.py`, `membrane_filter.py`,
`fasta_download.py`, `domain_download.py`, `ddi.py`, `dmi.py`, `idr_filter.py`, `monte_carlo.py`),
`microbiolink/utils/` (`uniprot_client.py`, `id_resolution.py`, `fasta.py`, `dmi_reader.py`), and
`microbiolink/data/` (Module 5/6 resource TSVs), all wired into `microbiolink/cli.py`.

Module 9 is the first module that (a) depends on an external algorithm package rather than resources we
ship, (b) requires an on-disk working directory (TieDie is file-in / folder-out), and (c) fetches its
core network live from OmniPath. It sits at pipeline position **DMI → (Module 7 IDR) → (Module 8 Monte
Carlo) → Module 9**, but its only hard predecessor is a DMI-shaped table: per the module-map, the DMI
input may be the Module 8, Module 7, **or** Module 6 output depending on how much filtering the user
wants, because Module 9 reads only the protein-pair columns that all three share.

**The three steps.** The refactoring plan specifies TieDie as a three-step process, matching the two
ground-truth scripts plus the algorithm package between them:

1. **Build inputs** — convert the DMI table + DEG list + expressed-gene list into TieDie's three input
   files (`upstream.input`, `downstream.input`, `pathway.sif`), plus the contextualised TF–target
   network that step 3 needs.
2. **Run TieDie** — invoke the TieDie algorithm on those three files; it writes an output folder
   containing `tiedie.cn.sif` (the causal subnetwork), `heats.NA` (linker heats), and diagnostics.
3. **Format outputs** — reassemble `tiedie.cn.sif` + `heats.NA` + the host–microbe interactions + the
   contextualised TF–target network + the DEG values into one **final network file** and one **node
   annotation file**.

**Ground-truth source.** Unlike Module 8, Module 9 **has** ground-truth code on the `case-study` branch,
and it is the reference:
- Step 1: `workflow/tiedie_input_processing.py` — fetches OmniPath PPI + CollecTRI TF–target networks,
  contextualises them to expressed genes, and writes `upstream.input`, `downstream.input`,
  `pathway.sif`, `contextualised_regulator-target_network.txt`, `contextualised_regulators_of_targets.txt`.
- Step 2: `workflow/TieDie/` — a vendored, Python-2-era copy of TieDie. **Not** ported; superseded by
  the modern pip-installable `saezlab/tiedie` (see Decision 1).
- Step 3: `workflow/processing_tiedie_output.py` — combines the TieDie output into the final network
  table and node table.

The refactored Module 9 is a **faithful port** of steps 1 and 3 (their transformation logic is the
contract), re-pointed at (a) the refactored DMI table schema and (b) the modern `tiedie` package for
step 2. `case_study_output/output/TieDIE/` (`usecase_upstream.input`, `usecase_downstream.input`,
`usecase_pathway.sif`, `usecase_final_network.txt`, `usecase_node_table.txt`) is the regression fixture
(per Q11).

**The TieDie package interface (verified against `saezlab/tiedie`, v2).** Modern Python package,
`requires-python >=3.9`, deps `numpy`/`scipy`/`networkx`, installable via
`pip install git+https://github.com/saezlab/tiedie.git`.
- **CLI:** `tiedie -n pathway.sif -u upstream.input -d downstream.input --output_folder OUT`, with
  `-s/--size` (size control, default 1.0), `-a/--alpha` (linker cutoff override), `-c/--depth` (causal
  path depth, default 3), `-p/--permute` (permutations, default 1000), `--pagerank`, `--pcst`,
  `-k/--kernel`. Programmatic entry point: `tiedie.cli.main(args_list)`.
- **Input `.sif`:** tab-separated `source <interaction> target` (e.g. `-a>`, `-t|`; MicrobioLink emits
  `stimulates>` / `inhibits>`, which TieDie's `classify_interaction` accepts).
- **Input `.input` (heats):** tab-separated `gene <heat> <sign(+/-)>`.
- **Output folder** contains `tiedie.sif`, `tiedie.cn.sif` (+`.txt`), `heats.NA`, `node_types.NA`,
  `node.stats`, `report.txt`, `score.txt`, `permuted_scores.txt`. Step 3 consumes `tiedie.cn.sif` and
  `heats.NA`.
- **Note the API gap that drives Decision 1:** the public `tiedie` API (`ScipyKernel`, `parse_heats`,
  `parse_net`, `normalize_heats`, `filter_linkers`, `find_linker_cutoff`, `connected_subnets`, …)
  covers diffusion and linker extraction, but the **causal-path filtering that produces
  `tiedie.cn.sif` lives only in `tiedie/cli.py`** (`find_consistent_paths`, `extract_subnetwork`) and is
  not exported. Reproducing it via the public API would mean re-porting that logic; driving the CLI
  gets it for free.

## Input schema adaptation — the key gap

The ground-truth step 1 reads a host–microbe interaction file with columns `# Human Protein`,
`Bacterial protein` (and an optional `sign`). The refactored DMI table (Module 6 `OUTPUT_COLUMNS`, extended
by Modules 7/8) instead has `dmi_type`, `bacterial_uniprot_id`, `human_uniprot_id`, … and **no `sign`
column**. Module 9 therefore reads the DMI table **as a DataFrame** and projects the two columns it
needs — `human_uniprot_id` (the upstream human target) and `bacterial_uniprot_id` (its microbial
partner) — rather than re-reading a differently-named on-disk file. `human_uniprot_id` is the human
target of the bacterial protein in **both** DMI directions (in `forward` rows it is the motif side, in
`reverse` rows the domain side), so no `dmi_type` branching is needed for the upstream set.

**DDI is the same shape.** The optional Module 5 DDI table (`ddi.OUTPUT_COLUMNS`) carries
`bacterial_uniprot_id`, `bacterial_pfam_domain`, `human_uniprot_id`, `human_pfam_domain`, `resource` —
so it too reduces to `(human_uniprot_id, bacterial_uniprot_id)` host–microbe pairs, exactly like the
DMI table. Module 9 therefore treats "host–microbe interactions" (HMI) as the **union of the DMI and DDI
protein-pair projections** and threads that single combined pair set through both the upstream heat
(Decision 5) and the step-3 network assembly. No other DDI column is used (the Pfam domains do not enter
the TieDie inputs). DDI is optional: when no DDI table is supplied, HMI is just the DMI pairs, identical
to before.

With no `sign` column present on either table, upstream signs follow the ground-truth **else-branch**:
every upstream node gets sign `-` and the heat is the interaction count per human protein (see Decision 5
for the exact count semantics). This is deterministic and matches the ground-truth behaviour when the HMI
file lacks a `sign` column.

## Shared code and reuse

- **`microbiolink/utils/fasta.py`** — not needed (Module 9 works on the DMI table, not sequences).
- **`microbiolink/utils/dmi_reader.py`** — Module 9 needs only `human_uniprot_id` / `bacterial_uniprot_id`,
  which are direction-independent, so it reads those columns directly and does **not** use
  `dmi_reader.motif_side` (that helper resolves the motif-bearing side, which is not what "human target"
  means here). No change to `dmi_reader`.
- **`microbiolink/utils/id_resolution.py`** — step 3 maps node UniProt accessions → human gene symbols
  (the reverse of `_translate_gene_symbols_to_uniprot`). Add a public `translate_uniprot_to_gene_symbols`
  helper here (single implementation, Q5-style), porting the ground truth's `batch_uniprot_to_gene_symbol`
  (MyGene, `scopes='uniprot'`, `fields='symbol'`, `species='human'`, batched), so both directions live in
  one module.
- **OmniPath access** — the human PPI and TF–target networks are fetched live (Decision 4). The lazy
  `import omnipath as op` inside the function body follows the Module 2 (`membrane_filter.py`) precedent.
- **`tiedie` package** — a new optional dependency, driven through its CLI (Decision 1).

## Decisions

The four pivotal decisions were resolved with the user before drafting:

1. **Step 2 runs the pip-installable `saezlab/tiedie`, driven through its CLI.** Add `tiedie` as an
   optional dependency (mirroring the `idr` extra) and invoke `tiedie.cli.main([...])` **in-process**
   (not a subprocess — no PATH/interpreter coupling, and `SystemExit` from its internal `sys.exit` is
   caught and re-raised as a `RuntimeError` with the captured stderr). This gets the causal-path
   `tiedie.cn.sif` — which the public API does not expose (see the API gap above) — for free, and tracks
   upstream fixes. The old vendored `workflow/TieDie/` copy is **not** ported. *(Subprocess invocation is
   the fallback if in-process `main()` proves to have unwanted global state; noted, not chosen.)*
2. **One orchestrator CLI command, `microbiolink-tiedie`, runs all three steps end-to-end.** The
   module-map says "the output should be the output from step 3", so the user-facing command produces the
   final network + node table directly. The three steps remain **separate public functions** in the
   workflow module (so a library caller can run them independently or swap the intermediate files), but
   only the orchestrator is wired to a console script. *(Per-step CLI commands are deferred; the split
   functions make adding them later trivial if wanted.)*
3. **DDI (Module 5) is supported as an optional upstream input, merged with the DMI table.** *(Changed
   from deferred at the user's request.)* Both tables reduce to `(human_uniprot_id, bacterial_uniprot_id)`
   host–microbe pairs, so Module 9 forms the combined HMI set as their **union** and threads it through
   the upstream heat (Decision 5) and the step-3 network assembly. DDI is **optional**: omit it and the
   behaviour is exactly the DMI-only path (which is the only path the case study exercises). Because the
   ground-truth scripts never handled DDI, the DDI-augmented path is **net-new functionality with no
   case-study fixture** and its output must be manually confirmed (per the refactoring plan's new-feature
   rule). DDI and DMI evidence are weighted **equally** — a distinct `(human, bacterial)` pair counts once
   regardless of which table(s) it came from, and a pair present in both is not double-counted; see the
   DDI note in @ai_docs/decisions/tiedie_inputs.md (Decision 1) for the equal-weight rationale and the
   confidence-weighting alternative.
4. **The human PPI and TF–target networks are fetched live from OmniPath each run**, exactly as the
   ground truth does: `op.interactions.OmniPath.get(genesymbols=1)` for the PPI and
   `op.interactions.Transcriptional.get(databases='CollecTRI', genesymbols=1)` for the TF–target network.
   `omnipath` is already a base dependency. *(Packaging a pinned OmniPath/CollecTRI snapshot for
   reproducibility and offline runs is recorded under Future Work.)*

The remaining decisions follow the refactoring precedence rules and the ground-truth logic:

5. **Upstream heat = count of distinct `(human_uniprot_id, bacterial_uniprot_id)` pairs per human
   protein; sign `-`.** *(Confirmed with the user over the raw-row-count and distinct-partner
   alternatives.)*

   How the ground truth actually computes this: `read_hmi_file` keeps only `# Human Protein` /
   `Bacterial protein` (no dedup), then `process_upstream_input` drops the bacterial column and does
   `groupby('# Human Protein').size()` — a **raw row count**. It is not a weighted sum: every DMI
   instance row contributes `+1`, and because the bacterial column is dropped before grouping, nothing is
   de-duplicated. A human protein hit by one bacterium through three motif/domain matches scores `3` from
   that bacterium alone. `normalize_heats` then rescales the absolute-value sum of all heats to 1000
   before diffusion, and the sign does not affect diffusion (see Decision 6 / the upstream-sign note).

   Module 9 **deliberately diverges**: it de-duplicates the **combined HMI set (DMI ∪ DDI, Decision 3)**
   to distinct `(human_uniprot_id, bacterial_uniprot_id)` pairs before counting, so the heat reads as
   "how many distinct microbial proteins target this human protein" rather than "how many predicted
   interaction instances land on it." The rationale: the refactored Module 6 table fans out per (motif hit
   × compatible domain × partner), so a raw count would let motif/domain multiplicity — an artefact of the
   prediction granularity, not biological signal — dominate the heat. Distinct-partner counting is the
   defensible semantics and is robust to that fan-out. A pair appearing in both the DMI and DDI tables is
   counted once (union, not sum), so the two evidence sources are weighted equally. The cost is that the
   upstream heats will **not** match the case-study
   `usecase_upstream.input` byte-for-byte (they match only up to the dedup); this is expected and
   accepted. Because `normalize_heats` rescales anyway, the downstream effect on network topology is
   second-order. *(If exact case-study fidelity is later preferred, switching to the raw-row count is a
   one-line change in `build_upstream_heats`; distinct-partner counting is the confirmed default.)*
6. **Downstream heats are the CollecTRI TF activity scores, ported verbatim from the ground truth.**
   For each contextualised TF, merge its DEG targets' log2FC, sign each target's contribution by the
   TF→target `consensus_stimulation` (stimulation keeps the sign, inhibition flips it), average across
   targets, and emit `TF <mean_signed_log2FC> <+/->`. This transformation is the contract; it is ported
   line-for-line, only re-parameterised on the refactored inputs.
7. **The `pathway.sif` is the expressed-gene-contextualised OmniPath PPI**, with `consensus_stimulation`
   mapped `1 → stimulates>`, `0 → inhibits>`, de-duplicated, written as headerless tab-separated
   `source <direction> target` (UniProt accessions). Ported verbatim.
8. **"Expressed genes" contextualise the PPI and TF–target networks and come from a transcriptomics
   file** (`--transcriptomics_file` + value column + separator), keeping a gene iff its value is present,
   non-NaN, and non-zero — the same non-zero expression filter Module 3 owns (Q4). This is the
   Module 1 z-score output (or any count matrix) supplied by the user. The DEG list
   (`--deg_file`) is separate: it is the *endpoint* set, used **as supplied** — see Decision 11.
11. **The DEG list is taken as-is; no p-value column or filtering is exposed.** *(Confirmed with the
    user.)* The ground-truth step 1 accepts an optional `endpoint_pvalue_column` and, when given, filters
    the endpoint genes to a hardcoded `< 0.05`; step 3 reads the same column but its filter line is
    commented out. Module 9 **drops this input entirely**: the user is expected to supply an
    already-significant DEG list (as the case-study `..._degs_fc05.csv` already is), so the tool does no
    p-value filtering of its own. This removes the `--deg_pvalue_column` CLI argument, the
    `endpoint_pvalue_column` parameter, and the `< 0.05` filter from the port. The `endpoint_value_column`
    (the log2FC column, Decision 6) is unaffected and still required.
9. **TieDie is a file-in / folder-out algorithm, so step 2 needs a working directory — the one
   `workflow`-module exception to "no file I/O in core".** The orchestrator writes the three input files
   into a working directory (a `tempfile.TemporaryDirectory` by default, or a user-given `--work_dir`
   for inspection), runs `tiedie.cli.main`, and reads the output folder back. This is unavoidable given
   TieDie's interface; it is confined to the two thin functions (`run_tiedie`, `run_tiedie_pipeline`) and
   documented. All *transformation* logic stays pure (DataFrame in, DataFrame out); only the TieDie
   hand-off touches disk.
10. **Requires the `tiedie` optional dependency; no new base dependency.** `omnipath`, `pandas`,
    `numpy`, `mygene` are already base deps. `tiedie` (and its `networkx`) is the one new install, gated
    behind the `tiedie` extra, exactly as `idr` gates Modules 7/8.

## Target implementation

### `microbiolink/utils/id_resolution.py` (extend)

```python
def translate_uniprot_to_gene_symbols(uniprot_ids: list[str]) -> dict[str, str]:
    """Map human UniProt accessions to gene symbols via MyGene (batched).

    Reverse of _translate_gene_symbols_to_uniprot. Accessions with no symbol are omitted; the caller
    falls back to the accession itself (matching the ground truth).
    """
```

- Lazy `from mygene import MyGeneInfo`; `querymany(..., scopes='uniprot', fields='symbol',
  species='human', returnall=True)`; batch in chunks of 100 (ground-truth batch size).

### `microbiolink/workflow/tiedie.py` (new, core, argparse-free)

```python
"""TieDie tied-diffusion network propagation: build inputs, run TieDie, format the output network."""

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from ..utils import id_resolution

UPSTREAM_SIGN_DEFAULT = "-"          # DMI/DDI tables carry no sign column (Decisions 2, 5)
HMI_PAIR_COLUMNS = ["human_uniprot_id", "bacterial_uniprot_id"]
STIMULATION_INTERACTION = "stimulates>"
INHIBITION_INTERACTION = "inhibits>"
NETWORK_COLUMNS = ["Source.node", "Relationship", "Target.node", "layer"]
NODE_TABLE_COLUMNS = [
    "node", "all_nodes", "bacteria_layer", "bindingprot_layer",
    "ppi_layer", "tf_layer", "deg_layer", "gene_symbol",
]


# --- OmniPath network fetch (live; Decision 4) -----------------------------------------------------

def fetch_omnipath_networks() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fetch the OmniPath PPI and CollecTRI TF-target networks (lazy omnipath import)."""


# --- Step 1: build TieDie inputs -------------------------------------------------------------------

def combine_host_microbe_interactions(
    dmi_table: pd.DataFrame,
    ddi_table: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Union the DMI and (optional) DDI tables into distinct host-microbe protein pairs (Decision 3).

    Projects each table to HMI_PAIR_COLUMNS (human_uniprot_id, bacterial_uniprot_id), concatenates, and
    drops duplicates — so a pair in both tables is kept once (equal-weight union, not sum). With
    ddi_table=None this is just the DMI table's distinct pairs.
    """


def contextualise_networks(
    ppi: pd.DataFrame,
    tf_target: pd.DataFrame,
    expressed_genes: list[str],
    endpoint_genes: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Restrict the PPI to expressed↔expressed edges and the TF-target network to expressed TFs of DEGs.

    Returns (contextualised_ppi, contextualised_tf_target). The TF-target table is also the step-3
    'contextualised regulator-target network'.
    """


def build_pathway_sif(contextualised_ppi: pd.DataFrame) -> pd.DataFrame:
    """Build the headerless source/<direction>/target SIF (consensus_stimulation -> stimulates>/inhibits>)."""


def build_upstream_heats(hmi_pairs: pd.DataFrame) -> pd.DataFrame:
    """Build the upstream heat table from the combined HMI pairs (combine_host_microbe_interactions):
    one row per human_uniprot_id, heat = distinct bacterial partner count (Decision 5), sign = '-'.
    Columns: [protein, heat, sign], no header."""


def build_downstream_heats(
    contextualised_tf_target: pd.DataFrame,
    endpoint_genes: pd.DataFrame,
    value_column: int,
) -> pd.DataFrame:
    """Build the downstream heat table: per TF, the sign-corrected mean DEG log2FC and its sign
    (Decision 6). Columns: [tf, heat, sign], no header."""


# --- Step 2: run TieDie (file-in / folder-out; Decision 1 + 9) -------------------------------------

def run_tiedie(
    upstream: pd.DataFrame,
    downstream: pd.DataFrame,
    pathway_sif: pd.DataFrame,
    work_dir: Path,
    size: float = 1.0,
    alpha: float | None = None,
    depth: int = 3,
    permute: int = 1000,
    use_pagerank: bool = False,
) -> Path:
    """Write the three inputs into work_dir, run tiedie.cli.main in-process, return the output folder.

    Raises RuntimeError if TieDie fails (SystemExit from its internal sys.exit is captured with stderr).
    """


# --- Step 3: format the output network -------------------------------------------------------------

def assemble_network(
    causal_sif: pd.DataFrame,
    hmi_pairs: pd.DataFrame,
    contextualised_tf_target: pd.DataFrame,
) -> pd.DataFrame:
    """Combine the three layers (bacteria-bindingprot, bindingprot-tf, tf-deg) into one edge table.

    The bacteria-bindingprot layer is built from the combined HMI pairs (DMI ∪ DDI), so DDI-derived
    host-microbe edges appear alongside DMI ones (Decision 3)."""


def assemble_node_table(
    whole_net: pd.DataFrame,
    heats: pd.DataFrame,
    endpoint_genes: pd.DataFrame,
    value_column: int,
) -> pd.DataFrame:
    """Build the per-node annotation table: layer membership, gene symbol (uniprot->symbol), linker heat,
    and DEG log2FC."""


# --- Orchestrator (the public entry point) ---------------------------------------------------------

def run_tiedie_pipeline(
    dmi_table: pd.DataFrame,
    endpoint_genes: pd.DataFrame,
    expressed_genes: list[str],
    endpoint_value_column: int,
    ddi_table: pd.DataFrame | None = None,
    work_dir: str | None = None,
    size: float = 1.0,
    alpha: float | None = None,
    depth: int = 3,
    permute: int = 1000,
    use_pagerank: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run all three TieDie steps and return (final_network, node_table).

    Combines the DMI table with the optional DDI table into the host-microbe pair set, fetches the
    OmniPath networks, contextualises them, builds the three TieDie inputs, runs TieDie in a working
    directory (a TemporaryDirectory if work_dir is None), reads tiedie.cn.sif + heats.NA back, and
    assembles the final network and node tables. The DEG table is used as supplied (no p-value filtering;
    Decision 11).

    Args:
        dmi_table: A Module 6/7/8 DMI table (only human_uniprot_id/bacterial_uniprot_id are read).
        endpoint_genes: The DEG table, already filtered to significant genes (first column = gene symbol).
        expressed_genes: Gene symbols with non-zero expression (contextualisation universe).
        endpoint_value_column: 1-based column of the log2FC/expression value in endpoint_genes.
        ddi_table: An optional Module 5 DDI table, merged into the host-microbe pair set (Decision 3);
            None runs the DMI-only path.
        work_dir: Working directory for TieDie I/O; a TemporaryDirectory if None.
        size, alpha, depth, permute, use_pagerank: TieDie parameters (passed through).

    Returns:
        (final_network, node_table) as data frames.

    Raises:
        RuntimeError: If TieDie fails.
    """
```

- `run_tiedie`: writes `upstream.input`, `downstream.input`, `pathway.sif` into `work_dir` (tab-separated,
  headerless, `index=False`); builds the arg list (`-n`, `-u`, `-d`, `--output_folder`, `-s`, `-c`,
  `-p`, optional `-a`, `--pagerank`); `import tiedie.cli` lazily; calls `tiedie.cli.main(arglist)` inside
  a `try/except SystemExit` that re-raises as `RuntimeError`; returns `work_dir / "TieDIE"` (or the
  `--output_folder` it passed).
- `run_tiedie_pipeline`: first `hmi_pairs = combine_host_microbe_interactions(dmi_table, ddi_table)`,
  then thread `hmi_pairs` into `build_upstream_heats` and `assemble_network` (everything downstream sees
  the combined DMI ∪ DDI pair set, Decision 3).
- `assemble_network` / `assemble_node_table`: ported from `processing_tiedie_output.py`
  (`process_edges` / `process_node_table` + `fetch_gene_symbols`), re-pointed at the combined HMI pair
  set (`human_uniprot_id`/`bacterial_uniprot_id` in place of `# Human Protein`/`Bacterial protein`)
  and calling `id_resolution.translate_uniprot_to_gene_symbols`.

### `microbiolink/cli.py` (extend, argparse only)

- `_add_tiedie_source_arguments(parser)`: `--dmi_file` (required — a Module 6/7/8 DMI CSV),
  `--ddi_file` (optional — a Module 5 DDI CSV, merged into the host-microbe pair set; Decision 3),
  `--deg_file` (required — the endpoint DEG table, pre-filtered to significant genes; Decision 11),
  `--deg_sep` (default `,`), `--deg_value_column` (`int`, required — 1-based log2FC column),
  `--transcriptomics_file` (required — expressed-gene source), `--transcriptomics_sep` (default `,`),
  `--transcriptomics_value_column` (`int`, required).
- `_add_tiedie_algorithm_arguments(parser)`: `--size` (`float`, `1.0`), `--alpha` (`float`, optional),
  `--depth` (`int`, `3`), `--permute` (`int`, `1000`), `--pagerank` (flag), `--work_dir` (optional —
  keep TieDie's intermediate files instead of a temp dir).
- `_build_tiedie_parser()`: the two helpers plus `--network_output` (required) and `--node_output`
  (required).
- `tiedie()` entry point: read the DMI CSV (`pd.read_csv`), the optional DDI CSV (`pd.read_csv` if
  `--ddi_file` given, else `None`), the DEG table, and the expressed-gene list (reuse the
  non-zero-expression read: value present, non-NaN, non-zero); call
  `tiedie_workflow.run_tiedie_pipeline(..., ddi_table=ddi_table)`; write the two returned tables with
  `to_csv(sep="\t", index=False)`; `return 0`.

### `pyproject.toml`

- Add `tiedie = ["tiedie @ git+https://github.com/saezlab/tiedie.git"]` under
  `[project.optional-dependencies]` (`allow-direct-references` is already `true`).
- Add `microbiolink-tiedie = "microbiolink.cli:tiedie"` under `[project.scripts]`.
- Note in install requirements that Module 9 requires the `tiedie` extra (Python `>=3.9`, pulls in
  `networkx`).

## Migration checklist

### 1. Shared-code extraction
- [x] Add `translate_uniprot_to_gene_symbols` to `microbiolink/utils/id_resolution.py`.

### 2. Core module
- [x] Create `microbiolink/workflow/tiedie.py` with `combine_host_microbe_interactions`,
      `fetch_omnipath_networks`, `contextualise_networks`, `build_pathway_sif`, `build_upstream_heats`,
      `build_downstream_heats`, `run_tiedie`, `assemble_network`, `assemble_node_table`,
      `run_tiedie_pipeline`.
- [x] Confirm `combine_host_microbe_interactions` unions the DMI and optional DDI tables to distinct
      `(human, bacterial)` pairs (equal-weight, a pair in both counted once), and passes through unchanged
      when `ddi_table=None` (Decision 3).
- [x] Confirm `build_upstream_heats` counts distinct (human, bacterial) pairs over the combined HMI set,
      sign `-` (Decision 5).
- [x] Confirm `build_downstream_heats` reproduces the ground-truth sign-corrected mean-log2FC TF score
      (Decision 6). *(Byte-exact 292/292 against the frozen `contextualised_regulator-deg_network.txt`.)*
- [x] Confirm `build_pathway_sif` maps `consensus_stimulation` 1→`stimulates>`, 0→`inhibits>`, headerless
      (Decision 7). *(Ported via `np.where` — pandas 2.2 rejects `.loc` string-into-int assignment.)*
- [x] Confirm `run_tiedie` writes the three inputs, runs `tiedie.cli.main` in-process, catches
      `SystemExit`, returns the output folder (Decisions 1, 9). *(Verified end-to-end with tiedie 2.1.0.)*
- [x] Confirm `assemble_network` / `assemble_node_table` reproduce the ground-truth three-layer network
      and node table, re-pointed at the refactored DMI columns and `translate_uniprot_to_gene_symbols`.
      *(Node-table `nan,nan` `ppi_layer` artifact deliberately cleaned to `NA` — see Decision appendix in
      @ai_docs/decisions/tiedie_inputs.md; a confirmed divergence from `usecase_node_table.txt`.)*

### 3. CLI wiring + packaging
- [x] Add `_add_tiedie_source_arguments`, `_add_tiedie_algorithm_arguments`, `_build_tiedie_parser`,
      `tiedie()` to `cli.py`.
- [x] Add the `tiedie` extra and `microbiolink-tiedie` script to `pyproject.toml`.

### 4. Delete old code (per Q11)
- [x] Delete `workflow/tiedie_input_processing.py`, `workflow/processing_tiedie_output.py`, and the
      vendored `workflow/TieDie/` directory once the port is validated. *(They were tracked in the
      `refactoring` tree after all; deleted via `git rm` after the end-to-end sign-off — nothing in
      `microbiolink/` referenced them.)*

### 5. Verification — case-study regression + manual sign-off
- [x] Run step 1 on the case-study inputs (the DMI/IDR table, `Enterocytes BEST4_degs_fc05.csv`, the
      transcriptomics file) and diff `upstream.input`, `downstream.input`, `pathway.sif` against
      `case_study_output/output/TieDIE/usecase_*.input` / `usecase_pathway.sif`.
      **Result:** `upstream` — same 76 proteins, heats differ only by the distinct-pair dedup (Decision 5),
      sign `-` vs fixture `+` (Decision 2). `downstream`/`pathway` do **not** match byte-exact against a
      *live* OmniPath fetch (downstream: all 292 fixture TFs present as a subset of 620, heats shifted by
      CollecTRI growth; pathway: 93% edge overlap, rest is OmniPath drift) — the "match exactly"
      expectation only holds against a *frozen* snapshot, which the downstream port does (292/292). This
      is the live-fetch caveat and motivates the pinned-snapshot future work (Decision 4).
- [x] Run the full `microbiolink-tiedie` orchestrator and diff the final network + node table against
      `usecase_final_network.txt` / `usecase_node_table.txt`. **Result:** exit 0; schemas identical;
      node table 0 duplicates, bacteria layer 36/36 exact, `nan` cleanup verified (0 vs 398 fixture rows);
      bacteria-bindingprot edges 89% of fixture. Interior (PPI/TF-DEG) is ~2× larger purely from live
      CollecTRI/OmniPath growth. Shown to the user for sign-off.
- [x] Confirm `run_tiedie_pipeline` works with a Module 6, Module 7, and Module 8 table as the DMI input
      (all three expose `human_uniprot_id`/`bacterial_uniprot_id`). *(End-to-end run used a Module-7-shaped
      IDR table; M6/M7/M8 `OUTPUT_COLUMNS` all verified to carry both pair columns.)*
- [x] **DDI path (net-new, no fixture):** run the orchestrator with `--ddi_file` set to a Module 5 DDI
      CSV and confirm the DDI human targets appear in the upstream heats and the DDI host-microbe edges in
      the final network's bacteria-bindingprot layer; that a `(human, bacterial)` pair present in both DMI
      and DDI is counted once; and that omitting `--ddi_file` reproduces the DMI-only result exactly.
      **Result:** all four properties confirmed with a synthetic DDI table against the real `tiedie.cn.sif`
      (new target heat +1; new edge present; dup pair counted once; DDI=None ≡ DMI-only). Shown to the user.

## Verification

- `uv run ruff format` / `uv run ruff check` scoped to the new/edited files (`workflow/tiedie.py`,
  `utils/id_resolution.py`, `cli.py`, `pyproject.toml` reviewed by eye).
- `uv run ty check`.
- The case-study regression (step 1 file diff) plus the end-to-end manual sign-off above.

## Future Work

- **Confidence-weighted DDI/DMI upstream heat** — weight each host-microbe pair by its Module 7/8
  evidence scores instead of the equal-weight distinct-pair count (see @ai_docs/decisions/tiedie_inputs.md
  Decision 1). The natural upgrade once DDI and DMI are both feeding the upstream set (Decision 3).
- **Pinned network snapshot** — package a versioned OmniPath PPI + CollecTRI snapshot (Decision 4) for
  reproducible, offline runs, replacing the live fetch.
- **Per-step CLI commands** — expose `microbiolink-tiedie-input` / `-run` / `-output` alongside the
  orchestrator (Decision 2) if users want to inspect or swap intermediate files.
