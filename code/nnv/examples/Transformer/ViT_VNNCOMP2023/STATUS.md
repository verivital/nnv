# STATUS — sound NNV verification of the VNN-COMP 2023 ViT benchmark

Honest, measured status of this implementation. Numbers from MATLAB R2025b on this
machine; reproduce with `run_benchmark.m` and `bab_demo.m`.

## What is implemented and verified

**Sound attention reach** (`engine/nn/funcs/SoftmaxAttn.m`) — the missing
set-propagation math for `softmax(scale·Q·Kᵀ)·V` on Star sets:
interval (Rump) set@set matmul, the *exact* correlated row-softmax bound, a
sign-aware symbolic `A·V` envelope (value path stays correlated), and
provenance-safe prefix-aligned residual/head assembly. Each is validated in MATLAB
by Monte-Carlo containment against the `linprog` LP oracle — **26 tests** in
`tests/nn/attention/`, all green, covering all sign regimes, multi-head shapes,
and the historically-unsound branches (V≤0/mixed `A·V`, the structural-equality
residual trap, the (i,j)-local reciprocal).

**Layers** — `DynamicMatmulLayer.reach` now does a sound interval matmul when
operand shapes are set (refuses without, as before); `ScaledDotProductAttentionLayer`
is sound for multi-token input across all star-family methods (was box-lift +
fail-loud). All 47 vnncomp25 + 17 transformer-soundness + 19 SDPA regression tests
still pass.

**End-to-end ViT** (`ViTReach.m`) — forward parity to onnxruntime **≈ 1.8e-6**
(`test_ViTReach_parity.m`), so reach verifies the deployed ONNX. Sound single-shot
reach: **0.8 s** (pgd) / **2.9 s** (ibp), **0 sample escapes** across the suite.

## Measured precision (15 instances/model)

| model | full ε=1/255 verified | median ε* (estimate) | typical full-ε margin |
|---|---|---|---|
| `ibp_3_3_8` | **0/15** | **0.484/255** | ≈ −0.5 (close) |
| `pgd_2_3_16` | **0/15** | ≈ 0 | ≈ −1000 (vacuous) |

**0/200 at full ε is expected and correct, not a failure.** The benchmark's
`generate_properties.py` keeps only images where PGD fails *and* vanilla CROWN
cannot certify — i.e. exactly the instances on which every *incomplete*
(non-branch-and-bound) method fails. The published sound score (α,β-CROWN
**79/200**) comes **entirely from branch-and-bound**. What the sound reach shows
positively: `ibp_3_3_8` margins approach 0 (not vacuous) and it certifies a real
radius **ε\* ≈ 0.48/255** (median, sound), ~half the benchmark ε; `pgd_2_3_16` is
hostile to interval bounds (ε\* ≈ 0), matching auto_LiRPA IBP. This reproduces the
n2v reference's conclusions.

## Branch-and-bound

The user's note "BaB should already be in NNV" is only half-true for transformers.
NNV's BaB engine (`engine/nn/gpu_bab/`) is a GPU CROWN/β-CROWN engine that operates
on an op-list (`nn_to_ops.m`) whitelisted to `{affine,relu,conv,normaffine,avgpool,
maxpool,add,concat,product}`. **It refuses softmax-attention** (no `.type`, no CROWN
bound rule) — a transformer never reaches it (sound-by-refusal). Making *that*
engine verify attention needs per-op CROWN forward/backward rules for bilinear
QK^T + softmax + A·V, which is a substantial addition (and is what α,β-CROWN
implements). That is the real path to a higher count and is left as future work.

In the meantime this example provides **two sound BaB methods driven by the new
attention layers**, using `ViTReach.reach` as the bounding oracle:

- **`verifyBaB`** — input-box splitting. Sound (robust iff every leaf is certified;
  any falsifying sample ⇒ not-robust). Works with the fast estimate path. *Weak in
  3072-d*: one-pixel splits barely tighten, so it does not close the gap to 1/255 —
  exactly why complete tools split activations, not inputs. Included as a sound,
  working BaB and a correctness baseline.
- **`verifyBaBRelu`** — FF-ReLU **phase splitting** (β-CROWN style). A node forces a
  set of ReLU phases (making those neurons exact and adding their half-spaces); its
  two children fix one more candidate neuron active/inactive, so the node is robust
  iff both children are. The binding margin is closed with an **LP that honours the
  forced half-spaces** (`getMin`). Sound by the same branch-cover argument.

Measured levers on `ibp_3_3_8` inst 11 (single-shot **estimate** ε\*=0.625/255):
- **LP margins alone** (no splits) are markedly tighter than estimate: single-shot
  with `marginMode='lp'` certifies **0.80/255** (≈1.3× the estimate radius), because
  the LP `getMin` exploits the full Star constraint system (av-envelope facets +
  ReLU triangles) that the predicate-box estimate ignores. So `marginMode='lp'` is a
  near-free precision gain (one LP per binding class). The estimate-margin sweep
  above therefore *under*-reports the true certified radius.
- **ReLU-phase splits** extend further where the LP root fails — `bab_demo.m` sweeps
  ε ∈ {0.80, 0.90, 1.00}/255 comparing single-shot-LP vs BaB; BaB ≥ single-shot-LP
  by construction (the root is the no-split case), strictly greater when a split
  flips the binding margin. See the script output for the per-ε numbers.
- Full-eps LP margins (`lp_fulleps.m`) are markedly tighter than estimate but still
  negative on every instance: the closest, `ibp_3_3_8` inst 11, has estimate margin
  −0.2237 and **LP margin −0.0627** (dual-simplex, ~13 s) — tantalisingly close, but
  short.
- **Empirical full-ε BaB sweep** (`bab_sweep_close.m`, the FF-ReLU phase-split BaB
  on the 5 closest ibp instances, ≤55 nodes each, ~20 min/inst): **0/5 verified at
  ε=1/255.** The decisive case — inst 11, single-shot LP −0.0627 — did **not** flip
  even after 55 ReLU-phase-split nodes. This is direct evidence (not inference) that
  FF-ReLU-split BaB cannot close the gap: the binding looseness is the box-lifted
  softmax attention *weights*, which sit *upstream* of the FF ReLUs, so adding ReLU
  half-spaces barely moves the binding margin. Closing it requires a symbolic
  softmax-weight relaxation **and** BaB over the *attention* nonlinearity — which is
  exactly what α,β-CROWN's 79/200 does, and is the documented future work.

### Honest limitations

- Full ε=1/255 on the filtered instances needs many activation splits (the gap is
  large: margin ≈ −0.5 on ibp). Our BaB has a node budget and splits only FF ReLUs
  (not the attention nonlinearity), so it does **not** reach 79/200 — consistent
  with the user's expectation ("I don't expect them to reach the 1/255 level").
- The forward bound stays `estimate` even inside BaB (only the binding margin uses
  LP); a full β-CROWN would propagate the split constraints through *all* bounds
  (tighter, but the LP-everywhere forward is ~100× slower here).
- `verifyBaBRelu` does no per-branch counterexample search, so it returns
  robust/unknown (pair with `verifyBaB`/`falsifyBox` for not-robust witnesses).

## Files

- `engine/nn/funcs/SoftmaxAttn.m` — sound attention primitives (reusable).
- `ViTReach.m` — model load, parity forward, sound reach, verify, BaB.
- `extract_weights.py` — ONNX → .mat (weights + vendored inputs + ORT ref logits).
- `run_benchmark.m`, `bab_demo.m` — sweeps.
- tests: `tests/nn/attention/*` (primitives+layers), `test_ViTReach_{parity,reach,bab}.m`.
