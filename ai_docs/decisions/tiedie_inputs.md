# Decision — TieDie input-heat strategies for Module 9

Records how Module 9 converts MicrobioLink data into the two TieDie **heat** inputs — the upstream
(bacterial-target) heats and the downstream (TF-activity) heats — the trade-offs of the chosen
strategies, and the alternatives worth revisiting. See @ai_docs/plans/module_9_tiedie.md for the
implementation (these correspond to that plan's Decisions 5, 6, and 11).

TieDie diffuses heat from an **upstream** node set and a **downstream** node set across a signalling
network and returns the linking subnetwork. Each heat file is tab-separated `gene <heat> <sign(+/-)>`.
Two modelling questions decide what those numbers are:

1. **Upstream:** how hot is each human protein that microbial proteins target, and what sign?
2. **Downstream:** how hot is each transcription factor inferred from the DEGs, and what sign?

Neither heat is dictated by biology alone — each is a modelling choice with consequences for which
subnetwork TieDie returns. `normalize_heats` rescales the **absolute-value sum** of each heat vector to
1000 before diffusion, so only *relative* magnitudes within a vector matter; the sign is used later, in
the causal-path consistency check, not in the diffusion itself.

---

## Decision 1 — Upstream (bacterial-target) heat = distinct `(human, bacterial)` pair count

**Chosen.** For each human protein, heat = the number of **distinct bacterial proteins** that target it
(equivalently, distinct `(human_uniprot_id, bacterial_uniprot_id)` pairs). Every human target gets sign
`-` (Decision 2). A human protein hit by 3 microbes = heat 3; hit by 1 microbe through 5 different
motif/domain matches = heat 1.

**Source is the combined DMI ∪ DDI set (equal weight).** The pairs come from the DMI table and,
when supplied, the optional Module 5 DDI table (both reduce to `(human, bacterial)` pairs). The two are
**unioned**, so a pair predicted by both DMI and DDI is counted **once**, not summed — DMI and DDI
evidence are weighted equally. When no DDI table is supplied the set is just the DMI pairs. This is the
module-map's "if DDIs are included" path (Module 9 plan Decision 3); it is net-new relative to the ground
truth, which only ever saw DMI/HMI pairs, so the DDI-augmented output has no case-study fixture and needs
manual confirmation.

**How the ground truth differs.** The original `tiedie_input_processing.py` does a **raw row count**
(`groupby('# Human Protein').size()`) with no de-duplication over a DMI/HMI-only file, so every
motif/domain-match row adds `+1`. Module 9 de-duplicates the combined DMI ∪ DDI set to distinct pairs
first.

### Why this one

- **Robust to Module 6 fan-out.** The refactored DMI table fans out one motif hit across every
  compatible domain × partner, so a raw count would be dominated by prediction granularity (how many
  ways a motif can match) rather than biological signal (how many microbes bind). Distinct-pair counting
  strips the fan-out.
- **Interpretable.** "How many distinct microbial proteins target this human protein" is a meaningful
  quantity; "how many predicted interaction rows mention it" is not.
- **Reproducible.** Independent of the exact fan-out multiplicity, so stable across DMI-table revisions.

### Pros and cons

| | |
|---|---|
| **Pro** | Immune to motif/domain fan-out; interpretable; reproducible across table changes. |
| **Pro** | Still rewards multiply-targeted hubs (a human protein hit by many microbes is legitimately central to the host–microbe interface). |
| **Con** | Diverges from the original — will not byte-match `usecase_upstream.input` (matches only after dedup). |
| **Con** | Remains a **degree-like** heat: a human protein targeted by many microbes gets very hot and can dominate the upstream diffusion, independent of interaction quality. |
| **Con** | Ignores **interaction confidence** entirely — a weak, borderline DMI counts the same as a strong, IDR- and Monte-Carlo-passing one. |
| **Con** | Treats all bacterial partners as equal (no weighting by how many microbes, how confidently, or which microbial proteins). |

### Alternatives that could be better

- **Confidence-weighted heat** *(most promising).* Weight each distinct interaction by its Module 7
  `combined_score` (IDR/binding) and/or Module 8 `monte_carlo_qvalue` (significance), summing weights per
  human protein instead of counting. This makes the upstream heat reflect *evidence strength*, not just
  target count — the natural way to let the earlier filtering modules inform TieDie. It also subsumes the
  DMI-vs-DDI weighting question: rather than the current equal-weight union, DMI and DDI interactions
  could carry different (or evidence-derived) weights. Cost: the tables must carry those score columns
  (only true if Module 7/8 ran; DDI has no comparable score today), and the scores need mapping onto a
  sensible heat scale.
- **Raw row count (original).** Keep the ground-truth behaviour. Pro: exact fixture fidelity. Con: lets
  motif/domain multiplicity — an artefact — drive the heat. Rejected as the default for that reason;
  trivially available as a one-line switch.
- **Degree-compressing transform** (`log1p(count)` or `sqrt(count)`). Keeps the ordering but compresses
  the hub disparity so a protein hit by 40 microbes doesn't outweigh the rest of the network 40:1. Cheap;
  worth pairing with either counting scheme.
- **Binary presence (heat = 1 for any targeted protein).** Removes upstream degree disparity entirely —
  pure topology, every entry point equal. Pro: no hub domination. Con: discards the "how many microbes"
  signal, which is real information.

---

## Decision 2 — Upstream sign = `-` for every human target

**Chosen.** Every upstream node is written with sign `-`, inherited from the ground truth's else-branch
(used whenever the interaction table has no `sign` column — which the refactored DMI table never does).

### Why this one

- **Deterministic and faithful.** Matches the original's fallback exactly; no per-interaction sign data
  exists in the DMI table to do otherwise.
- **Conservative reading.** Treats microbial binding as *perturbing / interfering* with host protein
  function, a defensible default when the true functional effect is unknown.

### Pros and cons

| | |
|---|---|
| **Pro** | Deterministic, reproducible, matches the ground truth. |
| **Con** | Biologically arbitrary — a bacterial protein binding a host protein can **activate or inhibit** it; a blanket `-` asserts inhibition for all. |
| **Con** | A uniform sign makes TieDie's causal-path consistency check **nearly non-discriminating on the microbial side** — every source pushes the same direction, so path validity is decided almost entirely by the downstream signs and the network edge signs. |

### Alternatives that could be better

- **Per-interaction functional sign** *(the real fix).* If a future module (or an external resource such
  as IntAct, or a mechanistic effect predictor) annotated each host–microbe interaction as
  activating/inhibitory, that sign would make the upstream causal reasoning meaningful. This is the
  natural upgrade the ground truth's `if 'sign' in hbps.columns` branch was already written to accept —
  the data just doesn't exist yet in the pipeline.
- **Leave upstream unsigned and diffuse magnitude only.** Skip the causal-path signing on the microbial
  side (which is currently near-vacuous anyway) and treat the upstream purely as diffusion sources. Fewer
  false claims of directionality. Con: loses the (weak) sign filtering entirely.

---

## Decision 3 — Downstream (TF) heat = mean sign-corrected log2FC of a TF's DEG targets

**Chosen.** For each contextualised TF, merge its DEG targets' log2FC, sign-correct each by the TF→target
regulation direction (activation keeps the sign, inhibition flips it), and **average** across targets.
That mean is the heat; its sign is the heat's sign. Ported verbatim from
`process_downstream_input`.

### Why this one

- **Controls for target count.** Because it is a **mean, not a sum**, a TF with 50 DEG targets and one
  with 3 are on the same scale — hub TFs do not automatically dominate (contrast the upstream count
  heat, which *is* degree-driven).
- **Coherent activity direction.** Sign-correcting by mode of regulation makes every target's
  contribution point the same way for a genuinely active TF, so the averaged sign reads as inferred
  activity direction.
- **Faithful and simple.** Matches the ground truth and the case-study fixture; no extra dependency.

### Pros and cons

| | |
|---|---|
| **Pro** | Mean-not-sum neutralises the biggest disparity source (regulon size). |
| **Pro** | Sign-correction yields an interpretable activity direction per TF. |
| **Con** | **log2FC is a log scale and the mean preserves magnitude** — a TF averaging log2FC ≈ 10 gets ~10× the heat of one averaging ≈ 1, and `normalize_heats` preserves that ratio. |
| **Con** | **Single-/few-target TFs are unstable** — one extreme target (the case study had log2FC ≈ 11) makes a TF very hot off a single noisy observation, with no averaging to temper it. |
| **Con** | **Signed cancellation** — a TF whose targets split up/down averages toward ≈ 0 and effectively drops out, even with large individual |log2FC|. |
| **Con** | **No confidence weighting** — ignores DEG p-values, regulon coverage, and how many targets support the estimate; a 2-target and a 200-target mean are treated identically. |

### Alternatives that could be better

- **A statistical TF-activity method (VIPER / decoupleR ULM/MLM)** *(most principled).* Replace the
  hand-rolled mean with an established TF-activity inference that weights targets by regulon confidence
  and mode of regulation and returns a normalised activity score with significance. Far more robust to
  the single-target and cancellation problems. Notably, **TieDie already ships a version of this**: its
  `--d_expr` + `-m/--min_hub` mode calls `master_reg.ActivityScores.findRegulators`, computing TF
  activity from the differential-expression file directly (with a minimum-regulon-size gate). Adopting
  that mode would hand the whole downstream-heat problem to the algorithm's own, better-tested routine.
  Cost: a heavier dependency surface / different CLI wiring.
- **Minimum-target-count gate.** Require a TF to have ≥ k DEG targets before it gets a heat, killing the
  single-target spikes. Cheap, directly targets the worst instability. (TieDie's `--min_hub` is exactly
  this idea for its own activity mode.)
- **Magnitude-compressing transform** (rank, `sign(x)·log1p(|x|)`, or winsorising log2FC). Compresses the
  log-scale disparity and caps extreme fold changes while preserving ordering and sign.
- **Weight by DEG significance.** Fold change alone ignores confidence; weighting each target by (or
  gating on) its adjusted p-value would downweight noisy calls. Pairs naturally with removing the input
  p-value filter (Decision 4) by moving significance from a hard cut to a soft weight.
- **Sign-only / ±1 heats.** Use only the inferred direction, all magnitudes equal — pure topology,
  removing every magnitude disparity. Con: discards the (real) signal that some TFs are more strongly
  perturbed than others.

---

## Decision 4 — DEG list used as supplied; no p-value filtering in the module

**Chosen.** Module 9 drops the ground truth's optional `endpoint_pvalue_column` / hardcoded `< 0.05`
filter and uses the DEG table as given. The user is expected to supply an already-significant DEG list.

### Why this one

- **Single responsibility.** Significance thresholding is the user's differential-expression step, not
  TieDie's; folding it in hardcodes a 0.05 cut the tool shouldn't own.
- **The original already half-abandoned it.** Step 1 applied the filter but step 3's copy of it was
  commented out — the pipeline was already inconsistent about enforcing it.

### Pros and cons

| | |
|---|---|
| **Pro** | Simpler interface; no hidden hardcoded threshold; user keeps control of significance. |
| **Con** | No safety net — an unfiltered table silently pulls non-significant genes into the endpoint set. |
| **Con** | Removes the one place fold change and significance could have been combined. |

### Alternatives that could be better

- **Keep an optional, *configurable* filter** (`--deg_pvalue_column` + `--deg_pvalue_cutoff`, default off).
  Restores the safety net without hardcoding 0.05, for users who pass a full DE table.
- **Soft significance weighting** instead of a hard cut — carry the p-value into the downstream heat
  (Decision 3 alternative) so weak DEGs contribute less rather than being included or excluded outright.

---

## Summary

| Input | Chosen strategy | Main weakness | Best alternative to revisit |
|---|---|---|---|
| Upstream heat | Distinct `(human, bacterial)` pair count over DMI ∪ DDI | Degree-like; ignores interaction confidence; DMI/DDI equal-weighted | Confidence-weighting by Module 7/8 scores |
| Upstream sign | Blanket `-` | Biologically arbitrary; near-vacuous causal filtering | Per-interaction functional sign from an effect resource |
| Downstream heat | Mean sign-corrected log2FC per TF | Log-scale disparity; single-target instability; cancellation | VIPER/decoupleR or TieDie's own `--d_expr` activity mode; min-target gate |
| DEG filtering | Use as supplied (no filter) | No safety net | Optional configurable p-value filter or soft weighting |

All four are **faithful-port defaults** chosen for fidelity, simplicity, and reproducibility against the
case study. The alternatives above are the paths to a more *statistically grounded* set of TieDie inputs
once the pipeline carries the extra signal (interaction confidence, functional effect, regulon
statistics) they need.

---

## Appendix — possible future development: boolean node-table layer columns

Out of scope for the heat inputs above — this concerns the step-3 **node-table output** formatting, not
the diffusion inputs — but recorded here for continuity.

**Current encoding.** The node-annotation table marks each node's layer membership with one column per
layer (`bacteria_layer`, `bindingprot_layer`, `ppi_layer`, `tf_layer`, `deg_layer`) whose **value is the
role's display name** (`bacteria_layer = "bacteria"`, `tf_layer = "tf"`, …) or `"NA"` when the node is
not in that layer. `all_nodes` is the comma-join of the non-`NA` values, so the label text is what makes
that roll-up self-describing. `ppi_layer` is the special case: a node can be a PPI **source**
(`"bindingprot and/or protein"`), a PPI **target** (`"protein and/or tf"`), or both, so it holds a
comma-joined string of whichever roles apply. Absent membership is `"NA"` (already cleaned up from the
ground truth's `nan,nan` artifact — a deliberate divergence from `usecase_node_table.txt`).

**Proposed simplification (deferred).** Replace the label-as-value scheme with **boolean** flag columns
(`True`/`False`, non-members `False`) and derive `all_nodes` from the column *names* of the `True` flags.
This removes the redundant column-name/value duplication (`bacteria_layer = "bacteria"`) and makes layer
membership a clean one-hot encoding.

**Open design question (why deferred).** The `ppi_layer` dual role has no single obvious mapping:
- **Single `ppi_layer` bool** — `True` if the node is in the PPI layer at all (source *or* target).
  Simplest, but drops the source/target distinction (which used vague "and/or" labels anyway).
- **Two bools `ppi_source_layer` / `ppi_target_layer`** — preserves the distinction losslessly, at the
  cost of an extra column and less friendly `all_nodes` tokens (`ppi_source,ppi_target`).

Secondary choices to settle if pursued: the exact `all_nodes` token names (short names like `bacteria`,
`bindingprot`, `tf`, `deg` vs. keeping the descriptive PPI labels), and whether the boolean columns keep
the `_layer` suffix. Since the node table is a manual-sign-off artifact (live OmniPath already makes it
non-byte-exact), this further divergence from the fixture is acceptable when the time comes.
