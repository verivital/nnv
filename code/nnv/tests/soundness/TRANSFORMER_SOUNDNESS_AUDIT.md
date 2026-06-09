# Transformer / Attention Soundness Audit — 2026-06-08

Autonomous audit of NNV's transformer-verification layers for **soundness** (does the computed
reachable set contain every actual output?), run against **R2026a** via the MATLAB MCP, with
empirical Monte-Carlo (MC) containment checks. Fixes the silent-unsoundness holes found; the
remaining gap (multi-token attention) is now **fail-loud** instead of silently wrong. The user's
actual verified models use FC-simulated attention and are unaffected.

All changes are on the `transformer` branch (pushed to `ttj`), **not** merged upstream.

---

## TL;DR

| Layer / path | Soundness status | Action |
|---|---|---|
| `BatchNormalizationLayer.reach` (ImageStar + Star) | **Sound** (0/3000 MC) | none — Nov-2025 "too-tight" bug already fixed by recent channel-aware changes |
| `Softmax.reach_star_approx_bounds` (the bounds path) | **Sound** (0/3000 MC, ∈[0,1]) | none |
| `SoftmaxLayer.reach`, intermediate (attention) | **WAS unsound** (silent identity) | **FIXED** — now routes to sound bounds, fails loud on unsupported types |
| ONNX loader mid-network softmax | **WAS unsound** (default IsFinalLayer=true) | **FIXED** — `load_nnv_from_mat` auto-marks non-final softmax |
| `SDPA` / `MHA` reach, seq_len > 1 | **Unsound** (returns V-bounds) | **GUARDED** — now errors `multiTokenUnsound` instead of silently degrading; sound multi-token bound deferred |
| Gelu / LayerNorm / RoPE / Sinusoidal / Learned PE | Sound by construction | none (positional are exact affine) |

---

## Fixes committed this session

1. **`SoftmaxLayer.reach` silent `catch → identity`** (commit `d5f3aa05e`)
   - The intermediate (`IsFinalLayer=false`) path called `Softmax.reach_star_approx(IS, method)`
     with `method=''` (no caller passes one); that dispatch errored on the empty method and a
     `try/catch` silently returned the **unsound identity passthrough**. So even an explicitly
     non-final softmax produced wrong (input-domain) sets.
   - Now calls `Softmax.reach_star_approx_bounds(IS)` directly (no fragile method string), and an
     unsupported set type **raises** rather than silently degrading.
   - Verified: `final=false` → bounds in **[0,1]**, MC-containment **sound=1** (3000 samples);
     `final=true` unchanged (identity, sound for argmax-on-logits specs). All 28 vnncomp25
     regression tests pass.

2. **ONNX loader did not mark non-final softmax** (commit `d5f3aa05e`)
   - `load_nnv_from_mat` built `SoftmaxLayer(name)` with the default `IsFinalLayer=true`, so a
     decomposed-attention softmax (`MatMul → Softmax → MatMul`) loaded as an unsound identity.
   - Added a position-aware post-pass marking every non-last `SoftmaxLayer` non-final. Erring
     toward non-final is **safe** (bounds are sound, just looser); only a mislabeled *non-final*
     softmax left as identity would be unsound, which this prevents.

3. **Attention multi-token fail-loud** (commit `27311daf1`)
   - `SDPA.compute_attention_bounds` and `MHA.compute_attention_bounds_single_head` return **V's
     bounds**, exact only for `seq_len=1` (softmax over one score == 1 ⇒ output == V). For
     `seq_len > 1` the output is a convex combination **across tokens** that can leave any single
     token's value range, so returning V's bounds is **unsound** — but the layers flatten and
     silently proceed.
   - Added guards in `reach_star_approx` / `compute_mha_bounds`: when configured with a per-token
     `ValueDim`/`EmbedDim` and the input is an exact multiple (>1 tokens), **raise**
     (`multiTokenUnsound`). No-op when single-token or dims unknown, so all 19 SDPA + 24 MHA
     single-token unit tests still pass.

4. **CI soundness regression** (commit `79fa2480e`)
   - `tests/soundness/test_transformer_soundness_regression.m` — 6 MC/fail-loud checks, <5s,
     7/7 pass, rides the matrix CI so soundness is **continuously verified in parallel**.

---

## What remains (flagged, not fixed — needs review / is research-grade)

### A. Sound multi-token attention bound (the real fix for seq_len > 1) — **HIGH**
The output for token *i* is `out_i = Σ_j softmax(scores_i)_j · V_j`, a convex combination of the
value rows. A **sound** (loose) bound, ignoring Q/K entirely, is per output-coordinate *d*:
`out_(i,d) ∈ [ min_j V_lb_(j,d) , max_j V_ub_(j,d) ]` (the bounding box of the value hull, valid
because attention weights live in the simplex). Implementing it needs the reach API to keep the
`[seq_len, dim]` structure (currently it flattens to one vector). Tighter results come from a
proper per-row softmax bound (refs below). Until implemented, the guard keeps it sound-by-refusal.

### B. Tighter softmax bounds (optional) — **MEDIUM**
Current `Softmax.compute_softmax_bounds` is a sound but loose **interval/box** over-approximation
(decoupled output predicates). Wei et al. 2023 give tighter **linear** bounds (exp-reciprocal /
log-sum-exp decomposition) that would reduce UNKNOWN verdicts on transformer specs. Drop-in target:
`Softmax.compute_softmax_bounds_tight` (currently a sampling placeholder).

### C. `matlab2nnv` softmax → PlaceholderLayer (identity) — **LOW**
The MATLAB import path treats `SoftmaxLayer` as a pass-through `PlaceholderLayer` (identity). Sound
for a final argmax-on-logits spec; **unsound** for a standalone *mid-network* MATLAB softmax (rare —
MATLAB self-attention is self-contained in `selfAttentionLayer → MHA`). Flag only.

### D. Untracked soundness tests not in CI — **LOW/MEDIUM**
`tests/soundness/test_SLM_layers_soundness.m` and `test_SoftmaxLayer_reach_soundness.m` are
untracked (not committed) so CI never runs them; they're also demonstrations (record unsoundness)
rather than clean pass/fail. The new `test_transformer_soundness_regression.m` supersedes them for
CI. Decide whether to commit/convert the SLM ones.

---

## End-to-end validation (2026-06-08, MCP / R2026a)

- **`examples/Transformer/MNIST/verify_mnist_vit.m` (FC-simulated attention): WORKS** — 4/5
  verified at ε=0.5 (matches the prior result), 0 errors. Confirms the softmax/loader/attention
  changes don't regress the real pipeline. Per-layer activation soundness also re-confirmed:
  Gelu 25/25, Sigmoid 2/2, Tanh 2/2, LayerNorm 3/3.

- **`verify_mnist_vit_attention.m` (REAL `MultiHeadAttentionLayer` at index 5): BROKEN —
  pre-existing, HIGH priority.** Symptoms: `Could not extract attention weights. Using identity
  projections` → `LayerNormalizationLayer.reach` receives an **empty ImageStar**
  (`estimateRanges` warns) → garbage output bounds (e.g. class lb −134.9, max-other ub 155.9 for
  ε=0.02) → **0/5 verified**. The user's own `debug_soundness.m` /
  `diagnose_attention_mismatch.m` / `analyze_bound_looseness.m` confirm this is a known-hard area.
  - **Not caused by this session's changes**: `matlab2nnv`, `LayerNormalizationLayer`, and the
    `ImageStar` path are untouched; the multi-token guard provably didn't alter flow (it didn't
    fire — see next point).
  - **Why the multi-token guard missed it**: weight extraction failed, so `compute_mha_bounds`
    took the `isempty(W_Q)` branch and set `EmbedDim = n`; the guard's `n > EmbedDim` test is then
    false, so a genuinely multi-token ViT looks single-token and the unsound bounds flowed on. The
    guard is therefore **best-effort** (catches multi-token only when dims are correctly set).
  - **Root causes to fix (for review):** (1) `MultiHeadAttentionLayer.parse` weight extraction for
    real ViT / ImageStar models; (2) the empty-ImageStar produced upstream of LayerNorm (likely
    the MHA ImageStar reach, possibly compounded by the `single`-precision conversion warning);
    (3) the sound multi-token attention bound (item A above). Until then, **real-attention ViT
    verification should not be trusted** — use the FC-simulated path.

## Corrected attribution

There is **no "Ben Wooding" PR** in `sammsaski/n2v` (PR authors are Kiguli, HCWDavid, sammsaski;
commit authors David Wang / Hanchen David Wang / Samuel Sasaki). The relevant transformer PR is
**#12 by Kiguli — "Add ~50 transformer-era layer types + expand Sphinx docs"** (open, 2026-05-18).
Worth reviewing #12's layer set against NNV's for soundness parity. (If "Wooding" work exists, it's
in a different repo/fork — please point me at it.)

## References (the papers this is/should be based on)

- **Wei, Wu, Wu, Chen, Barrett, Farchi (2023)** — *Convex Bounds on the Softmax Function with
  Applications to Robustness Verification*, AISTATS. arXiv:2303.01713. (The sound/tight softmax
  bounds; exp-reciprocal + log-sum-exp decompositions.)
- **Shi, Zhang, Chang, Huang, Hsieh (2020)** — *Robustness Verification for Transformers*, ICLR.
  arXiv:2002.06622. (First closed-form linear self-attention bounds; cross-position dependency.)
- **Bonaert, Dimitrov, Baader, Vechev (2021)** — *Fast and Precise Certification of Transformers*
  (DeepT), PLDI.
- **Vertex-Softmax (2026)** — *Tight Transformer Verification via Exact Softmax Optimization*,
  arXiv:2605.10974. (Very recent; exact softmax.)
- α,β-CROWN / auto_LiRPA (linear relaxation framework NNV's bounds gesture at).

## How to re-verify locally (MCP)

```matlab
cd code/nnv; startup_nnv;
runtests('tests/soundness/test_transformer_soundness_regression.m')   % 7/7 expected
```

*Audit by Claude (Opus 4.8), autonomous session 2026-06-08. Changes verified on R2026a; soundness
claims are empirical (MC-containment), not formal proofs — treat the multi-token bound (item A) as
the priority for a rigorous follow-up.*
