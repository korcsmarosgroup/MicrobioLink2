# Decision — Monte Carlo null-model strategy for Module 8

Records the null-model options considered for Module 8's motif-significance filter, the option chosen
for the current implementation, why the alternatives were set aside, and which alternatives are worth
revisiting later. See @ai_docs/plans/module_8_monte_carlo.md for the resulting implementation plan.

Module 8 asks: is a predicted domain–motif interaction's motif *genuinely* present, or is it an
artefact of the protein's amino-acid composition? Every option below is a way of building the null
against which the observed motif is judged. The design evolved in steps: pick the over-representation
family → restrict to disordered regions → pool the background across proteins → sample from the
pooled composition rather than shuffle it. The end point is **option (a)** below.

## Final decision (current implementation) — Option (a)

**Monte Carlo over-representation against a pooled disordered-region composition.**

Mechanic:

1. Restrict attention to **disordered regions** (where SLiMs occur), not whole proteins.
2. Build **one background amino-acid frequency vector** by tallying residues across the disordered
   regions of all proteins in the pipeline (leave-one-out: exclude the protein under test from its
   own background).
3. For a motif instance in a protein whose disordered region has length `D`, Monte Carlo **sample**
   `x` synthetic length-`D` sequences by drawing each position i.i.d. from the pooled frequency
   vector, re-run the *same* motif regex from Module 6 (treated as a **black box**), and count
   matches per sample.
4. p-value = `(hits + 1) / (iterations + 1)`, where a hit is a sample whose regex match count is
   `>=` the count observed in the protein's real disordered region.

Why this one:

- **Over-representation, not disorder relocation.** It tests whether the motif *sequence* is special,
  with motif specificity — a property of the motif — as the dominant driver, rather than a protein
  length / disorder-fraction confound (see Method 1 below).
- **Disordered-region background, not whole-protein.** The disordered region is the composition the
  motif actually competes against; a whole-protein background is diluted by the ordered core and
  biases *toward* false significance — the wrong direction for a filter (see rejected option 2).
- **Pooled background, not single-region.** A large pooled background is a stable estimate and
  removes the unfair penalty a noisy single-region background places on small / non-specific motifs
  (see rejected option 3).
- **Sample the composition, don't shuffle the pool.** Sampling from the frequency vector is the
  efficient, canonical form of the identical composition-only null (see rejected option 4).
- **Black-box regex.** Re-running the Module 6 regex on synthetic sequences avoids any pattern
  introspection, making this the simplest option to implement correctly *now*.

## Background: the two families first considered

Two distinct nulls were considered at the top level. They sound similar but test different things,
have opposite length biases, and depend on almost unrelated variables.

### Method 1 — Disorder relocation

- **Null:** the motif's *position* and the protein's disordered regions are independent — any overlap
  is chance.
- **Mechanic:** keep the motif fixed; slide the disordered blocks (length-preserving) to random
  positions; each shuffle, check whether the motif window still lands in disorder.
- **Driver:** the disorder *fraction* `D/L`. Approximately `p ≈ (D − m + 1) / (L − D + 1)` for a
  motif of length `m`, disordered block `D`, protein length `L`. Longer protein → *more* significant;
  smaller disordered region → *more* significant (strongest lever).
- **Bias:** dominated by **protein length / disorder fraction** — a property of the protein, not
  evidence the interaction is real. That is the confound that rules it out.

### Method 2 — Motif over-representation (CHOSEN FAMILY)

- **Null:** the motif is no more frequent than expected given the relevant amino-acid composition.
- **Mechanic:** compare the observed motif count against a composition-matched null (by shuffling or,
  as chosen, by sampling from a composition).
- **Driver:** `E[matches] ≈ (D − w + 1) × p_match`, where `w` is the motif window width and `p_match`
  the per-position match probability under the background composition. Motif specificity dominates
  (`p_match` drops geometrically per constrained residue); region length enters only *linearly*.
- **Bias:** dominated by **motif specificity**; length is a gentle linear factor.

|                          | Disorder relocation (M1)              | Motif over-representation (M2)          |
|--------------------------|---------------------------------------|-----------------------------------------|
| Null                     | motif position vs. disorder is chance | motif over-represented given AA comp.   |
| Dominant driver          | protein length / disorder fraction    | motif specificity                       |
| Longer protein/region    | **more** significant                  | weakly **less** significant             |
| Motif length/specificity | weakly more significant               | strongly more significant               |

Method 2 answers the biologically relevant question ("is this motif's sequence special?"), so the
over-representation family was chosen. Everything below is a refinement *within* Method 2.

## Options considered and set aside (for now)

1. **Method 1 — disorder relocation.** Rejected: its dominant driver (protein length / disorder
   fraction) is a protein property, not evidence the interaction is real. It also requires a disorder
   profile without answering the sequence-specificity question we care about.

2. **Whole-protein shuffle.** Shuffle the whole protein's residues, re-run the regex. Rejected:
   disordered regions are enriched in SLiM-typical residues (P, E, S, K, R, Q) and depleted in bulky
   hydrophobics; the whole-protein composition is diluted by the ordered core, making those residues
   look rarer than they are in context. This *understates* `p_match` for SLiMs → *inflates*
   significance → lets spurious motifs pass. Wrong direction for a filter. (It was the original plan
   because it needs no disorder profile; that convenience does not justify a biased null.)

3. **Single-region shuffle / composition.** Use only the disordered region the motif sits in.
   Rejected: a single (often short) region gives a noisy composition estimate, and small /
   non-specific motifs — whose `p_match` is dominated by a few common-residue frequencies — are
   hypersensitive to that noise. A self-inclusion effect compounds it (the motif's own residues,
   preserved under shuffling, easily re-form a short motif), pushing small motifs toward
   non-significance. Both are artefacts of estimating the background from the one region the motif is
   in. Short regions also give coarse, granular p-values. Pooling (option a) fixes all three.

4. **Whole-pool shuffle.** Permute the entire concatenated pool of disordered residues and read
   length-`D` windows from it. Rejected as *redundant*, not wrong: a permutation destroys everything
   except composition, so a window from the shuffled pool is (for a large pool) indistinguishable
   from an i.i.d. draw from the frequency vector — the *same* composition-only null as option (a).
   The frequency vector already holds all the information the shuffle would materialise, so the
   permutation is wasted work. It is also expensive: building a null *distribution* would require
   re-shuffling the pool per iteration (O(pool × iterations)) or sliding correlated windows over a
   single realisation. Option (a) is the efficient canonical form of this same null.

5. **Analytic frequency method (no masking).** Compute `p_match` in closed form from the concrete
   match (class mass per constrained position, wildcards = 1) and use a Binomial/Poisson tail.
   Deferred, not rejected on merit: it is cheaper, exact, deterministic, and resolves the small
   p-values that multiple-testing correction needs — but it requires **regex introspection**
   (identifying constrained vs. wildcard positions and their class mass), which is the fiddly,
   error-prone part. Monte Carlo keeps the regex a black box, so it is simpler to ship now. This
   method returns as a future candidate, paired with masking — see below.

6. **Proteome-wide disordered background (pool option b).** Build the background from disordered
   regions across the whole proteome. Deferred: it is the unbiased, principled background, but needs
   whole-proteome disorder prediction (heavy precompute) or a shipped precomputed table. The
   pipeline-pool background (option a) is far cheaper, reuses disorder we already compute, and its
   only bias — being enriched for the motifs under test — makes `p̂` mildly *conservative*, which is
   the safe direction for a filter. Worth upgrading to (b) later; see below.

## Dependencies and implications of the chosen option

- **Reintroduces the disorder dependency.** To know which residues are disordered (for both the pool
  and each target's region), Module 8 reuses the existing disorder-profile helpers in
  `idr_filter.py` (promoted to public for sharing: `cached_profile`, `iupred_profile`,
  `aiupred_profile`). This reverses the earlier "Module 8 needs no IUPred/AIUPred install" property.
  Acceptable because the realistic pipeline is Module 7 (IDR) → Module 8, so the dependency is already
  present, and `cached_profile` means the profile is near-free if Module 7 ran in-process.
- **Module 8 becomes a post-IDR step.** Running it standalone on the 8-column Module 6 table is
  dropped: "sample the disordered region" is undefined without a disorder profile.

## Multiple-testing correction — Benjamini–Hochberg FDR (resolves open questions 1 and 2)

Module 8 tests many motif instances in one run, so raw per-test p-values are corrected for multiple
testing with **Benjamini–Hochberg (BH)**, controlling the **false discovery rate (FDR)** — the
expected fraction of reported interactions whose motif is spurious. FDR is the right criterion for a
discovery filter; Bonferroni / FWER (bounding *any* false positive) is far too strict over thousands
of instances and would gut sensitivity for a guarantee this step does not need.

Decisions:

- **`alpha` is repurposed as the target FDR** (default `0.05`), not a per-test p-value threshold. A
  motif passes iff its BH-adjusted q-value is `<= alpha`. Both the raw `monte_carlo_pvalue` and the
  adjusted `monte_carlo_qvalue` are reported; the pass/fail decision uses the q-value.
  - The q-value is computed **independently of `alpha`**: BH adjusts the ranked p-values, and `alpha`
    only sets the cutoff the q-value is compared against. Changing the FDR target changes which rows
    pass, never the q-values themselves — so q-values are computed once and can be re-thresholded
    without recomputation.
  - `0.05` is the familiar default and matches the previous per-test threshold. `0.10` is a common
    looser choice for discovery screens where sensitivity matters more; user-configurable either way.
- **The unit of testing is the unique motif instance `(motif class, protein, disordered region)`, not
  the DMI row (this resolves open question 2).** Module 6's fan-out repeats one instance across many
  (domain × partner) rows, all carrying an *identical* p-value; correcting over rows would inflate the
  test count `m` with exact duplicates and distort the BH ranking. BH runs over de-duplicated
  instances, and each instance's q-value is broadcast back to every row that shares it. A protein with
  several motif instances contributes several independent tests — the natural reading of "test each
  instance independently."
- **BH, not Benjamini–Yekutieli.** Instances are near-independent; the only coupling is the shared
  leave-one-out pool, a weak *positive* dependence within BH's PRDS validity. BY's
  arbitrary-dependence guarantee would be needlessly conservative here.
- **The MC p-value floor caps power, so `iterations` must scale with the test count.** The smallest
  p-value MC can produce is `1/(iterations + 1)`, while BH's threshold for the most significant of `m`
  instances is `≈ alpha/m`; resolving that tail requires roughly `iterations ≳ m/alpha`. The
  per-unique-null caching makes high iteration counts affordable — cost is per unique
  `(regex, region-length, composition)` null, not per instance. MC's discreteness (p-values are
  multiples of `1/(iterations + 1)`) also makes BH mildly conservative; accepted, discrete-aware BH
  variants are out of scope.

**Structural consequence.** BH needs the full p-value vector before it can threshold anything, so the
row-at-a-time "score → threshold → drop" flow is replaced by a two-phase **"score every instance →
BH-correct → filter"**. See @ai_docs/plans/module_8_monte_carlo.md for the implementation.

This is also the strongest argument for the deferred **frequency method with masking** below: its
continuous p-values have no floor, removing the `iterations ≳ m/alpha` blow-up that BH exposes here.

## Future options worth revisiting

- **Frequency method with low-complexity masking.** The analytic frequency method (option 5) is the
  natural cheap upgrade: exact, deterministic, continuous p-values that resolve small tails for
  multiple-testing correction. Its main statistical weakness is overlapping-window clumping, which is
  worst in low-complexity / repetitive disordered regions — exactly what the pooled MC also glosses
  over. Adding **low-complexity masking** (the standard SLiMSearch/QSLiMFinder-style mitigation)
  addresses that weakness directly. Cost: the regex introspection deferred in option 5, plus a
  masking pass. This is the most likely next iteration.

- **Un-shuffled pool (real contiguous-window sampling).** Instead of sampling i.i.d. from the pooled
  *composition*, sample real contiguous length-`D` windows from the *un-shuffled* pooled disordered
  regions. This is the one variant that carries more information than the frequency vector: it
  preserves local sequence structure (low-complexity runs, di-peptide correlations), giving a
  higher-order, more faithful null. Cost: hold the pool of real sequences and sample windows from it,
  and accept a stronger (less conservative) null. Good when composition-only is too crude.

- **Proteome-wide background (pool option b).** Replace the pipeline pool with an unbiased
  proteome-wide disordered background (precomputed table or whole-proteome disorder prediction).
  Removes the mild enrichment bias of option (a).

## Open questions

Both original open questions (multiple-testing correction; multiple motifs per protein) are resolved
above under "Multiple-testing correction". None remain open.
