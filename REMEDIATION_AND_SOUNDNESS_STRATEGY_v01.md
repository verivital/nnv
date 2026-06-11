# PR #290 — Soundness Remediation & Future-Detection Strategy

Companion to `CR_TRANSFORMER_CLAUDE.md` (the finding list). This file is the **how**: the
principled bound approach, the execution sequence, the automated detection mechanism, and the
precision/scalability balance. Formal-verification rule: **every reach() must be a sound
over-approximation** (computed set ⊇ all true outputs), verified empirically by MC-containment
and, where tractable, by analytic argument — *but not so loose as to be useless*.

## 1. The bound philosophy (precision ↔ scalability)

For each nonlinear layer we use, in order of preference:
1. **Exact** when the op is affine over the set rep (positional encodings, FC, conv, +const).
2. **Sound linear relaxation** (CROWN / auto_LiRPA style): two affine functions
   `g_lo(x) ≤ f(x) ≤ g_hi(x)` valid on the per-neuron interval `[l,u]`, encoded with one fresh
   error predicate so input↔output correlation is preserved (tight). For a C² activation the
   linearization error of a center line is bounded **soundly** by the chord + second-derivative
   bound: on a subinterval of width h, `|f(x) − chord(x)| ≤ M·h²/8` with `M ≥ max|f''|` — so a
   *sampled* error max plus the `M·h²/8` correction is a guaranteed enclosure (the bug in SiLU/
   GELU was sampling **without** this correction → unsound between samples).
3. **Sound interval box** (fresh decoupled predicate) when (2) is not yet derived — sound but
   looser; acceptable as an interim, never as a silent identity.
4. **Fail loud** (`error`) when no sound bound is implemented for that case — never return a
   set that might exclude a real output.

Never: sample-only error bounds, identity for an active op, exact-star labels on
over-approximations, or a guard that a default constructor disarms.

## 2. The detection strategy (so we never miss these again)

**`tests/soundness/soundness_harness.m`** — a generic property-based MC-containment driver:
for a `(layer, makeInputSet, evaluateFn)` triple it draws many random + corner + adversarial
concrete points from the input set, runs `layer.reach`, and asserts every `evaluateFn(point) ∈
reach`. Crucially it uses **adversarial input sets** (wide, nonuniform centers/radii, sign-
straddling, multi-position) — the existing tests missed bugs by using tiny uniform boxes.

- **`tests/soundness/test_all_layers_soundness.m`** runs the harness over every layer NNV
  supports, with a per-layer input-set spec. A reach that is unsound **or** that errors on a
  supported input fails the suite. This is the regression net + the dev-time check.
- Dev rule: **adding/❝fixing❞ any reach() requires a harness entry**; CI runs it (and it is NOT
  allow-listed — see the `ci_allowed_failures.txt` rule). The oracle must be correct first, so
  `evaluate()` bugs (MHA softmax axis/head-split) are fixed before trusting MC results.

## 3. Execution sequence (soundness-impact order)

1. **CI gate real** (done): skips≠pass, 0/0 fails, soundness tests un-allow-listed.
2. **evaluate() correctness** (the MC oracle): MHA softmax-axis/head-split, DynamicMatmul
   broadcast, ElementwiseAffine align.
3. **Build the harness** + seed it with every transformer/activation layer.
4. **Production-path soundness**: SiLU error-bound (+ SwiGLU/RMSNorm), Addition multi-input,
   exact-star labeling.
5. **Norm grouping + sequence parity**: LayerNorm per-position group + reachSequence; Softmax
   ImageStar group; PlaceholderLayer.evaluateSequence.
6. **Attention**: SDPA/MHA guards (ValueDim, empty-weights, ImageStar K/V, zono path),
   check_soundness oracle, then the **sound multi-token bound** (value-hull box + per-row
   softmax bound, Shi 2020 / Wei 2023 / DeepT) to replace fail-loud where tractable.
7. **GELU** erf-vs-tanh; **Concat/Reshape** layout; **matlab2nnv** topology/dead-branch;
   **Embedding** random-weights.
8. **VNN-COMP**: all benchmarks via the Python importer with **sound** op handling or principled
   refusal (cctsdb/collins/vit need Slice/Expand/Where/ScatterND/ArgMax/Split/Pow or multi-token
   attention); runner counterexample-order + multi-spec try/catch.
9. **Tests/CI**: assumeFail→verifyError, tolerances, figure pollution; investigate
   test9_NN_Segmentation.
10. **Re-review** (fresh code-review pass) to confirm the soundness problems are gone.

## 4. Status (updated as work proceeds)

- DONE: LayerNorm sound bound (`84cb90ea1`); PlaceholderLayer UnsupportedOp (`9978cef50`); CI
  gate real (`11a24df8f`); consolidated review (`aae9e3b14`).
- See the bottom of this file / git log for per-item progress.

## 5. Literature anchors for sound bounds

- Softmax: **Wei et al. 2023** (AISTATS, arXiv:2303.01713) convex/concave bounds.
- Self-attention: **Shi et al. 2020** (ICLR, arXiv:2002.06622) closed-form linear bounds;
  **Bonaert et al. 2021** DeepT (PLDI) abstract transformers; **Vertex-Softmax 2026**
  (arXiv:2605.10974) exact softmax.
- Activations (SiLU/GELU/Swish) linear relaxation: CROWN / auto_LiRPA (Zhang et al.), and the
  chord + |f''| enclosure for the sampled-error correction.
