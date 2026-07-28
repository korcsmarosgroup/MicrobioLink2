# Decision — Monte Carlo shuffling strategy for Module 8

Records the two candidate Monte Carlo null hypotheses considered for Module 8, why they differ, and
which one was chosen. See @ai_docs/plans/module_8_monte_carlo.md for the resulting implementation
plan.

## Context

Module 8 tests whether a predicted domain-motif interaction's motif is "real" via a Monte Carlo
shuffle. Two distinct nulls were considered. They sound similar ("is this motif–disorder association
by chance?") but test different things, have opposite length biases, and depend on almost unrelated
variables. The disorder-relocation null (Method 1) is what an earlier standalone script on the
`development` branch implemented; the motif-shuffle null (Method 2) is what Module 8 now uses.

## Method 1 — Disorder relocation

- **Null:** the motif's *position* and the protein's disordered regions are independent — any overlap
  is chance.
- **Mechanic:** keep the motif fixed; take the disordered regions as length-preserving blocks and
  slide them to random positions (~1000×); each shuffle, check whether the motif window still lands
  in disorder.
- **p-value:** `(hits + 1) / (iterations + 1)`, where a hit is a shuffle whose motif-window overlap
  is `>=` the observed overlap.
- **Needs a disorder profile:** yes (IUPred / AIUPred).
- **What drives it:** the disorder *fraction* `D/L`. As an approximation for a single disordered
  block of length `D` in a protein of length `L` and a motif of length `m`:
  `p ≈ (D − m + 1) / (L − D + 1)`.
  - Longer protein → **more** significant (more empty sequence for the block to miss).
  - Smaller disordered region → **more** significant (strongest lever).
  - Larger motif → weakly more significant.
- **Bias:** dominated by **protein length / disorder fraction** — properties of the protein, not
  evidence the interaction is real. That is the confound.

## Method 2 — Motif shuffle (over-representation) — CHOSEN

- **Null:** the motif is no more frequent than expected given the protein's amino-acid composition.
- **Mechanic:** composition-preserving shuffle — randomly permute the motif-bearing protein's
  residues (the multiset of amino acids is fixed, only order is destroyed), then re-run the motif
  regex.
- **p-value:** `(hits + 1) / (iterations + 1)`, where a hit is a shuffle whose regex match count is
  `>=` the count observed in the real sequence.
- **Needs a disorder profile:** no — nothing calls IUPred / AIUPred; no `method`, threshold, binding,
  or device flags; no `idr` dependency.
- **What drives it:** `E[matches] ≈ (L − w + 1) × p_match`, where `w` is the motif window width and
  `p_match` the per-position match probability under the protein's composition.
  - Shorter protein → slightly **more** significant (fewer positions to hit by chance) — length
    enters only *linearly*.
  - More specific / longer motif → **strongly** more significant — `p_match` drops *geometrically*
    per constrained residue, so motif specificity dominates protein length.
- **Bias:** dominated by **motif specificity**; protein length is a gentle linear factor.

## Why the two are not interchangeable

| | Disorder relocation (Method 1) | Motif shuffle (Method 2) |
|---|---|---|
| Null | motif position vs. disorder is chance | motif over-represented given AA composition |
| Dominant driver | protein length / disorder fraction | motif specificity |
| Longer protein | **more** significant | **less** significant |
| Motif length/specificity | weakly more significant | strongly more significant |
| Needs disorder profile | yes | no |

Neither is length-neutral, but they are biased by almost unrelated things and answer different
questions: Method 1 asks *"is this motif's position special relative to disorder?"* (blind to whether
the motif sequence is common); Method 2 asks *"is this motif's sequence special for this protein?"*
(blind to where it sits or whether it is accessible). A motif can pass one and fail the other.

## Decision

**Module 8 implements Method 2 (motif-shuffle over-representation).**

Rationale:
- It targets the question of interest — whether a predicted motif is genuinely over-represented rather
  than an artefact of amino-acid composition — with the dominant variable (motif specificity) being a
  property of the motif itself, not a length confound of the protein.
- It removes the disorder-profile dependency from Module 8 entirely, so the module runs without any
  IUPred / AIUPred install and does not recompute disorder that Module 7 already handles.
- Method 1's dominant driver (protein length / disorder fraction) is largely a protein property rather
  than evidence the interaction is real, making it the weaker null for this filter.

The two nulls are near-orthogonal, so they remain complementary: a motif that clears both would be
significant for two nearly-uncorrelated reasons. Running both as independent filters is a possible
future extension, but Module 8's current scope is Method 2 only.
