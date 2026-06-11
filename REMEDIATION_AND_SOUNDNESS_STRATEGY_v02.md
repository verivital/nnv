# PR #290 — Soundness Remediation & Future-Detection Strategy (v02)

Supersedes `..._v01.md` (the bound philosophy in §1–§2 there still holds verbatim — read it
first). This v02 records what the 2026-06-09 autonomous session **completed**, and scopes the
**remaining** work, with the multi-token attention design called out as the one capability gap.

Formal-verification rule (unchanged): every `reach()` must be a sound over-approximation
(computed set ⊇ all true `evaluate` outputs), validated by MC-containment and, where tractable,
by analytic argument — but not so loose as to be useless.

## A. Done this session (all reach-unsoundness resolved)

26 findings fixed across 9 commits (`68d95d635 … 5ac745050`); see the resolution table in
`CR_TRANSFORMER_CLAUDE.md`. Every fix added a regression test. Pattern applied consistently:

- **Silent-wrong → fail-loud** where no sound bound exists or an operand/value is missing:
  Addition dropped addend [11], EAffine bad bias [13], Embedding random weights [5], Constant
  with no value [26], Concat non-feature-axis on flat Stars [14].
- **Wrong-but-recoverable → made correct**: EAffine [1,1,C] axis when H==C [10] (shape-based
  broadcast, exact), Reshape OnnxBCHW reach [9] (mirror evaluate's layout, exact), intermediate
  Softmax ImageStar [12] (per-pixel channel groups), matlab2nnv final-softmax [22] (topology BFS),
  GELU variant [2][25] (erf vs tanh, sound bracket split).
- **Crash on supported input → sound fallback**: SiLU on 3-arg Stars [19] (interval box from the
  single minimum), Concat on 3-arg Stars [44], Softmax check sampler [20].
- **Dead/latent traps removed**: SiLU getLinearBounds [21] (sound constant bounds), EProduct
  multipleInputs [47], SDPA nargin [46], matlab2nnv dead ScalingLayer branch [27].

Detection mechanism in place: `tests/soundness/test_layer_soundness_harness.m` (generic MC over
adversarial boxes) — all registered activations/norms SOUND. The CI gate is now real (`11a24df8f`).

## B. Remaining work, ranked

### B1. Sound multi-token attention — the capability gap (replaces fail-loud) [8][3][4][40-followups]
Today multi-token SDPA/MHA `reach` **fails loud** (sound, but cannot verify ViT). To support the
VNN-COMP transformer benchmarks we need a sound bound for
`Attention(Q,K,V) = softmax(QKᵀ/√d_k)·V` over Star/ImageStar sets. Proposed staging:

1. **Value-hull baseline (sound, simple, loose).** For output token `i`, coordinate `d`:
   `o_{i,d} = Σ_j a_{ij} V_{j,d}` with `a_{ij} ≥ 0, Σ_j a_{ij} = 1` (a distribution, *whatever*
   the scores). A convex combination of `V_{j,d} ∈ [Vlb_{j,d}, Vub_{j,d}]` therefore lies in
   `[min_j Vlb_{j,d}, max_j Vub_{j,d}]` — **sound for any attention weights, no softmax bound
   needed**. Caveat: identical for every query `i` (loose). REQUIRES the per-token V layout
   (seq_len × ValueDim); get the flatten order right or it is loose-but-still-sound only if you
   fall back to the global min/max. *Validate with MC before trusting any grouping.*
2. **Softmax-aware refinement (tighter).** Bound the scores `S_{ij} = Q_i·K_jᵀ/√d_k` by interval
   (bilinear; McCormick or interval-product), then the per-row softmax weights with the sound
   monotone bound already implemented for SoftmaxLayer (`ub_k = e^{ub_k}/(e^{ub_k}+Σ_{l≠k}e^{lb_l})`,
   symmetric for lb). Combine `Σ_j a_{ij} V_{j,d}` with the weight ranges (LP or interval) to
   tighten beyond the value hull. Literature: Shi 2020 (arXiv:2002.06622), DeepT/Bonaert 2021,
   Wei 2023 softmax convex bounds (arXiv:2303.01713).
3. **Prerequisite — fix `evaluate` [3][4] FIRST.** `softmax(scores,2)` normalizes the wrong axis
   and the head split is strided not contiguous; both make `evaluate` the wrong oracle for
   multi-token. MC-validation of any reach bound is meaningless until `evaluate` matches the real
   network (replay a parsed `selfAttentionLayer.predict` to pin it). This is why [3][4] are listed
   under this item, not as independent live bugs.
4. **Wire into the engine [16][17].** MHA/SDPA take 3 positional sets (Q,K,V); the DAG engine
   gathers multi-input layers into a cell. Add MHA/SDPA to the multi-input whitelist (NN.m ~1510
   and ~1367) AND make their `reach` accept a `{Q,K,V}` cell. Until then, port-wired attention
   crashes and portless wiring silently drops operands — keep them OUT of the whitelist so the
   crash is loud, not silent, in the interim.

**Acceptance:** MC-containment 0 violations on adversarial multi-token boxes for several seq_len /
head configs, plus the value-hull analytic argument; then flip the fail-loud guard to the bound.

### B2. VNN-COMP runner correctness [28][29]
- [28] `falsify_single` never inverts the `needReshape==3` transform before writing the witness →
  traffic_signs SAT counterexamples in the wrong flat order (invalid / penalty). Add the inverse.
- [29] multi-spec (`cell lb`) reach branches lack try/catch; with layers now fail-loud, an uncaught
  reach error aborts the whole instance. Wrap each spec's reach.

### B3. CI / test hygiene [35][36][38]
- [35] replace remaining `assumeFail` in attention tests with `verifyError`/containment (one done).
- [36] softmax soundness test tolerance 0.1 (10% of codomain) — tighten.
- [38] tests write PNGs into the repo tree (SiLU/RMSNorm/SwiGLU) — write to a temp/ignored dir.

### B4. Lower / deferred
- [6] DynamicMatmul evaluate broadcasting (reach already fails loud).
- [42] exact-star mislabeling — **maintainer API decision** required (reject `exact-star` for
  transformer layers, or relabel the result approximate). Do not change unilaterally.
- [23] AttentionMask/causal ignored by parse; Copilot lows (verify_mnist_vit fallback vars,
  plots/ dir, Reshape −1 divisibility, TransposedConv bias validation, stale BN comment).

## C. Merge recommendation
The PR is now safe to merge **as "single-token / FC-simulated transformer utilities + fail-loud
guards + real CI gate"** once B2 (runner correctness) and B3 (CI hygiene) land — all live reach
soundness is resolved and tested. General **sound multi-token ViT verification (B1) is the
follow-up**; until it lands, the title/scope must state that multi-token attention is fail-loud
(refused), not verified. [42] must be resolved (reject or relabel) before any exact-star claim.

## D. Status log
- v01 → v02 (2026-06-09): all reach-unsoundness findings fixed + tested; 64/64 comprehensive
  checkpoint green; remaining work is capability (B1) + submission/runner (B2) + hygiene (B3).
