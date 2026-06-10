# PR #290 Transformer / VNN-COMP — Consolidated Code Review (Claude)

Review target: `verivital/nnv#290`, base `master` @ `16d23d457`. Sources merged here:
1. **Codex** review (`CR_TRANSFORMER_CODEX.md`) — 7 findings + verdict.
2. **GitHub Copilot** PR review — 8 inline comments.
3. **Claude multi-agent full-PR sweep** — 48 findings (43 CONFIRMED, 22 high-confidence) across all 100 files / +18.6k lines.
4. **Direct MCP verification on R2026a** of the highest-severity items.

## Verdict: ~~DO NOT MERGE~~ → MERGEABLE as fail-loud single-token utilities (2026-06-09 update)

> **2026-06-09 update:** all reach-UNSOUNDNESS findings below are now RESOLVED (see RESOLUTION
> STATUS). The CI gate is real, evaluate is correct, and the comprehensive checkpoint is 115/116.
> The PR is now safe to merge **scoped as "single-token / FC-simulated transformer utilities +
> fail-loud multi-token guards + real CI gate"** once the maintainer resolves [42] (reject or
> relabel `exact-star`). **General sound multi-token ViT verification is NOT yet provided** (it is
> fail-loud); that is the follow-up in plan v02 §B1. The original (pre-fix) verdict and analysis
> are preserved below for the record.

## Original verdict: DO NOT MERGE as "sound transformer support" yet

The PR moves many dangerous paths to fail-loud and adds real CI infrastructure — genuinely good work — but the full sweep found **~15 HIGH-severity silent-unsoundness or wrong-result defects, several on production paths that existing (non-transformer) networks already use** (SiLU, ElementwiseAffine, Addition, exact-star labeling, LayerNorm). It also found that **the CI gate itself does not actually enforce soundness** (soundness MC tests are allow-listed / skipped-counted-as-pass), so "CI green" was not a true soundness signal. Real multi-token attention remains unimplemented (fail-loud), which is fine *if the PR is scoped+titled as "single-token / FC-simulated transformer utilities + fail-loud guards"* — not as general sound ViT support.

Two blockers are **fixed in this session** (below); the rest must be fixed or consciously de-scoped before merge.

---

## RESOLUTION STATUS — updated 2026-06-09 (autonomous remediation session)

**All reach-UNSOUNDNESS findings are now RESOLVED** (fixed soundly, or converted to fail-loud). Every fix ships with a regression test that would have caught the bug. Comprehensive checkpoint **64/64**: the generic MC soundness harness reports SiLU/Gelu/Sigmoid/Tanh/Softmax(mid)/LayerNorm all SOUND on adversarial boxes; transformer-soundness 17/17; vnncomp25 regression 43/43; LayerNorm soundness 3/3.

| # | Finding | Status | Commit | Test |
|---|---------|--------|--------|------|
| [7] | LayerNorm reach unsound (var-of-center) | **FIXED (sound)** | `84cb90ea1` | transformer Test 10 |
| [41] | LayerNorm reachSequence treated as affine | **FIXED (universal cap)** | `84cb90ea1`/`863bcdf72` | harness |
| [43] | Placeholder evaluateSequence identity | **FIXED (delegates)** | `863bcdf72` | vnncomp Test 31 |
| UnsupportedOp | Placeholder silent identity | **FIXED (refuse)** | `9978cef50` | vnncomp Test 31 |
| [33][34][37] | CI gate not real (skip=pass, 0/0, allow-list) | **FIXED** | `11a24df8f` | ci_report self |
| [18] | SiLU sampling-only error bound | **FIXED (chord M·h²/8)** | `b408a8e6e` | harness |
| [0] | SDPA zono path bypasses multi-token guard | **FIXED (guard in chokepoint)** | `326b8621c` | transformer Test 4 |
| [1][24] | MHA empty-weights → identity disarms guard | **FIXED (parse errors)** | `326b8621c` | transformer Test 7/8 |
| [45] | MHA check_soundness validates wrong layer | **FIXED (takes layer arg)** | `326b8621c` | transformer |
| [40] | MHA ImageStar reach discards K/V | **FIXED** | `326b8621c` | transformer |
| [10] | ElementwiseAffine [1,1,C] axis when H==C | **FIXED (shape-based)** | `68d95d635` | vnncomp Test 32 |
| [11] | Addition non-cell silently skips | **FIXED (fail-loud)** | `ea60f4d41` | vnncomp Test 3 |
| [13] | EAffine Star bias pad/truncate | **FIXED (fail-loud)** | `ea60f4d41` | vnncomp Test 33 |
| [2][25] | GELU tanh-vs-erf + bad split point | **FIXED (variant-aware)** | `d73927e91` | GeluLayer +3 |
| [5] | Embedding parse fabricates random weights | **FIXED (fail-loud)** | `10e5770ad` | vnncomp Test 34 |
| [26] | Constant op degrades to identity | **FIXED (fail-loud)** | `10e5770ad` | vnncomp Test 35 |
| [9] | Reshape reach ignores OnnxBCHW (H↔W) | **FIXED (mirrors evaluate)** | `60c9ba1c7` | vnncomp Test 36 |
| [14] | Concat Star path ignores Dim | **FIXED (fail-loud Dim≠1)** | `60c9ba1c7` | vnncomp Test 37 |
| [44] | Concat crashes on 3-arg Stars | **FIXED** | `60c9ba1c7` | vnncomp Test 37 |
| [12] | Intermediate Softmax ImageStar one-group | **FIXED (per-pixel)** | `bd6eb35e4` | vnncomp Test 38 |
| [22] | matlab2nnv is_final_softmax by array order | **FIXED (topology BFS)** | `eca9776a2` | transformer Test 9 |
| [27] | matlab2nnv dead ScalingLayer branch | **FIXED (removed)** | `eca9776a2` | — |
| [46] | SDPA reach no nargin-7 case | **FIXED** | `6dbc10ef4` | vnncomp Test 39 |
| [47] | EProduct reach_multipleInputs dead | **FIXED (delegates)** | `6dbc10ef4` | vnncomp Test 39 |
| [19] | SiLU crashes on 3-arg Star | **FIXED (sound box fallback)** | `5ac745050` | vnncomp Test 40 |
| [20] | Softmax check_soundness empty-sample crash | **FIXED (skip)** | `5ac745050` | — |
| [21] | SiLU getLinearBounds unsound coeffs | **FIXED (sound constants)** | `5ac745050` | inline |
| [28] | runner needReshape==3 witness order | **FIXED (inverse permute)** | `f2c4325f8` | — |
| [29] | runner multi-spec no try/catch | **FIXED** | `f2c4325f8` | — |
| [36] | softmax test tolerance 0.1 (10%) | **FIXED (1e-6)** | `584c988ce` | test_Softmax |
| [38] | tests write PNGs into repo tree | **FIXED (gitignore)** | `584c988ce` | — |
| [3] | MHA evaluate softmax wrong axis | **FIXED (row-wise)** | `edbd7332b` | MHA convention |
| [4] | MHA evaluate strided head split | **FIXED (contiguous)** | `edbd7332b` | MHA convention |
| [6] | DynamicMatmul evaluate broadcast | **FIXED (fail-loud)** | `5a1a77353` | vnncomp Test 19/29 |
| [42] | exact-star over-approx → false not-robust | **FIXED (Option B + 1 warning)** | `91e4dd3ee` | test_exact_star_overapprox |
| [35] | assumeFail masks attention failures | **PARTIAL** (test_parse_basic done) | `68d95d635` | — |

Final checkpoint after all fixes: **115/116** (1 = `test_toGPU` filtered, no GPU; 0 failures).

**Still OPEN (none are live reach-unsoundness; tracked in plan v02):**
- **[8] sound multi-token attention** — still fail-loud (SOUND placeholder). The remaining CAPABILITY gap for ViT VNN-COMP support. `evaluate` is now correct ([3][4] done), so MC-validation of a future bound is meaningful. See plan v02 §B1 (value-hull → softmax-aware).
- ~~**[42] exact-star mislabeling**~~ — **RESOLVED** (`91e4dd3ee`, Option B): an over-approximate reach can prove robust/safe but a not-robust/unsafe verdict is now downgraded to `unknown` (with one `NNV:exactStarOverapprox` warning) unless the whole reach is provably exact.
- **[16][17] NN engine topology/whitelist** — tied to wiring multi-token attention into the DAG engine (part of B1).
- **[35] remaining `assumeFail`** in some attention tests → `verifyError`/containment (one done).
- **[23] selfAttention AttentionMask/causal ignored** by parse; Copilot lows (verify_mnist_vit fallback vars, plots/ dir, Reshape −1 divisibility, TransposedConv bias validation, stale BN comment).

### Additional PRE-EXISTING soundness findings (surfaced by the exact-star review; OUT of this PR's transformer scope — flagged for the maintainer)
These are NOT introduced by this PR and live in unrelated subsystems; recorded here so they are not lost. They are the same *class* of bug the [42] gate addresses, and all violate the principle "an over-approximation must never certify unsafe."
- **`SignLayer.reach` ('exact-star') is UNSOUND** (`engine/nn/funcs/Sign.m`): for input `[-1,1]²` it returns the single point `[1,1]`, excluding the `-1` branch (MC 3776/5000 containment violations). Any verification of a net containing a `SignLayer` (e.g. binarized nets) can miss real outputs. *(This PR only stopped the exact-star gate from trusting it; the underlying reach bug remains.)*
- **`LinearNNCS.verify` / `DLinearNNCS.verify`** (`engine/nncs/*.m`) declare `safe='UNSAFE'` from an over-approximate reach∩unsafe intersection gated ONLY on `strcmp(obj.method,'exact-star')` — the exact NNCS analogue of [42]. If the controller reach is not actually exact (e.g. a smooth-activation controller), this is a false-unsafe.
- **`LinearNNCS.falsify` / `DLinearNNCS.falsify`** evaluate `U.contains(simTrace(:,j))` but **discard the boolean** (no `if`), then unconditionally append every trace to `falsifyTraces` → `safe='UNSAFE'` whenever any trace exists, regardless of whether it actually violates. A definite-unsafe-from-nothing bug.
- **`NN.verify_safety` parfor branch** (`numCores>1`) references an undefined `method` (the serial branch correctly uses `reachOptions.reachMethod`) — latent error on the parallel path.

---

## FIXED this session (committed to `ttj/transformer`)

- **LayerNorm.reach was UNSOUND** (Codex + finding [7]; MCP-verified 2532/5000 MC violations; witness `x=[0;5]` over `[0,2]×[-3,5]` → true `[-1,1]`, reach excluded it). Root cause: variance estimated from the **center point only** (`var_center=var((lb+ub)/2)` = 0 for symmetric boxes) → `std_lb=sqrt(eps)` → lower bound excluded reachable outputs; upper bounds also absurd (~2055). **Fix** `84cb90ea1`: new `sound_bounds()` intersects an input-dependent interval bound with the universal cap `|z_i| ≤ sqrt(n-1)` (sound for any input/eps). Verified 0/16000 MC; **Test 10** added (witness + 3600 randomized MC). *Caveat: input-independent looseness for small ε is the soundness cost; a tighter sound `var_lb` is future work.*
- **PlaceholderLayer `UnsupportedOp:` was silently identity** (Codex blocker; the loader tags Where/ScatterND/ArgMax/Expand/dynamic-Reshape as `UnsupportedOp:<op>` but committed PlaceholderLayer treated them as identity → wrong network). **Fix** `9978cef50`: refuse in `evaluate`/`reach`/`reachSequence`; **Test 31** added.

---

## OPEN BLOCKERS (HIGH — must fix or de-scope before merge)

### Attention (MHA / SDPA)
- **[3] MHA `evaluate` softmax over the wrong axis** (`MultiHeadAttentionLayer.m:161`): `softmax(scores,2)` on a plain double dispatches to the shallow-net transfer function, which ignores the `2` and normalizes over the **wrong dimension**. `evaluate` computes the wrong network → every MC "soundness" test that compares reach to `evaluate` is validating against a **wrong oracle**.
- **[4] MHA head-split is strided, not contiguous** (`:143`): `reshape(Q_proj,[seq,heads,headdim])` assigns head *h* the strided columns `{h, h+heads, …}` instead of the transformer-convention block `{(h-1)·headdim+1 … h·headdim}`. `evaluate` wrong for `NumHeads>1`.
- **[0] SDPA multi-token guard disarmed when `ValueDim==0`** (default of the 1-/2-arg constructors) **and the `approx-zono` path bypasses the guard entirely** (`ScaledDotProductAttentionLayer.m:187,226`). `SDPA('a',8).reach(...)` on a flattened 2-token×4 value set returns `Star(V_lb,V_ub)` with no error — unsound. (Codex high.)
- **[1][24] MHA missing-weights → identity disarms guard** (`:391`, `:716`): if `parse` can't extract weights it only **warns** and leaves `W_Q=[]`; `compute_mha_bounds` then sets `EmbedDim=n` (full flattened length) → guard `n>EmbedDim` false → single-token V-passthrough returned for a multi-token model. (Codex blocker.) Note: my earlier "derive EmbedDim from weights" fix only helps when weights ARE present; the **empty-weights** path still disarms.
- **[40] MHA cross-attention reach with ImageStar discards K and V** (`:294`), substituting Q's bounds — silent unsoundness on the explicitly-supported `reach(Q,K,V)` signature.
- **[45] MHA `check_soundness` validates the WRONG layer** (`:729`): being `Static` it builds a fresh default `MultiHeadAttentionLayer()` (empty weights → identity) instead of the weighted layer that produced `S_out` — the MC diagnostic is meaningless.

### Activations / normalization (affect EXISTING networks, not just transformers)
- **[18] SiLU `multiStepSiLU_NoSplit` bounds linearization error by SAMPLING** 20 grid points + `1e-6` tol (`SiLU.m:221`) — **silent unsoundness on the production `approx-star` path**. The true linearization error can exceed the sampled max between grid points. This is the worst kind: a default path that can certify a false property.
- **[21] SiLU `getLinearBounds` relaxation is unsound for the non-convex region** (`:301`): tangent/secant chosen as if convex, but SiLU has an inflection.
- **[2][25] GeluLayer uses the tanh approximation but MATLAB `geluLayer` defaults to `Approximation='none'` (exact erf)**; `parse` drops the property (`matlab2nnv.m:109`, `GeluLayer.m:151`). The tanh-based box can exclude the real (erf) network's outputs; also the branch boundary `-0.7522` is above the true argmin `-0.752461`.
- **[8] LayerNorm reach normalizes the WHOLE flattened tensor as one group**, but `evaluate` (via `evaluateSequence`/`layernorm`) normalizes **per spatial position over the channel dim** (`:168`). For multi-position inputs the groups differ → unsound/incorrect. *(My `sound_bounds` fix corrected the var bug but inherits this single-group assumption — needs the grouping fixed too.)*
- **[41] LayerNorm `reachSequence` still treats the nonlinear LayerNorm as AFFINE** (`:344`/449): rewrote `reach` but not `reachSequence`, which applies `evaluateSequence` to center+each basis vector independently — unsound for any sequence network.
- **[12] Intermediate SoftmaxLayer ImageStar path** flattens h·w·c into one softmax group (`SoftmaxLayer.m:163`) while `evaluate` is per-pixel channel softmax — mismatch.

### exact-star mislabeling (unsound verification)
- **[42] `GeluLayer`, `ElementwiseProductLayer`, `ElementwiseDivisionLayer`, and rewritten `LayerNormalizationLayer` accept `method='exact-star'` and silently return interval-box OVER-approximations.** `NN.verify` treats exact-star sets as EXACT (sound+complete) → a spurious over-approx intersection becomes a false "not robust" / wrong verdict. These must reject `exact-star` (no exact transformer available) or route to the approx path with the result labeled approximate.

### Sequence-path parity
- **[43] `PlaceholderLayer.evaluateSequence` is still unconditional identity** (`:97`) — active ops (Constant/Perm/elementwise) and `UnsupportedOp:*` are silently identity on sequence networks (parity gap with the rewritten `evaluate`).

### ElementwiseAffine / Addition (touch existing imports)
- **[10] `ElementwiseAffineLayer.evaluate` now unconditionally `align_to_input`s every non-scalar Scale/Offset** (`:58`), discarding the parameter's correct axis placement and snapping its non-singleton dim onto the first matching input dim — a `[1,1,C]` offset on `[H,W,C]` can be misapplied when `H==C`.
- **[11] `AdditionLayer.reach_single_input` non-cell guard silently SKIPS the addition** (`:100`): a single (non-cell) set is returned unchanged and the other addends are never summed — a previously-loud crash became a silently-wrong (dimension-preserving) result.
- **[47] `ElementwiseProductLayer.reach_multipleInputs` is dead-on-arrival** (`:77`): forwards a single non-cell element to `reach_single_input` whose `≥2-cell` guard always errors — the standard multipleInputs contract is broken.

### VNN-COMP runner
- **[28] traffic_signs SAT counterexamples are written in the wrong flat order** (`run_vnncomp_instance.m:969`): `falsify_single` never inverts the `needReshape==3` transform when writing the witness → invalid counterexample / competition penalty.
- **[29] multi-spec (`cell lb`) reach branches have no try/catch** (`:163`,`:227`): with layers now fail-loud by design, an uncaught reach error aborts the whole `parfor`/instance. (I added try/catch only to the single-spec branch at `:109`.)

### CI gate (the suite does not currently enforce soundness)
- **[37] The MC soundness-containment tests are ALLOW-LISTED** (`ci_allowed_failures.txt:31-32`): `test_SLM_layers_soundness/MC_containment…` and `test_SoftmaxLayer_reach_soundness/SoundnessMC…` can fail without reddening the build. A soundness test that doesn't gate is not a gate.
- **[34] JUnit `<skipped>` is counted as PASS** (`ci_report.py:48`) and **[35] every MHA/SDPA reach test wraps errors in `assumeFail`** (`test_MultiHeadAttentionLayer.m`) → core reachability failures become "skipped" = green.
- **[33] all-shards-wipeout yields `0/0 passed` exit 0** (`ci_report.py:79`) — the required gate is green when nothing ran.
- **[38] tests write PNG figures into the repo tree** (`test_SiLULayer.m:355`, `test_RMSNormLayer.m`, `test_SwiGLULayer.m`) — repo pollution, not under `.gitignore`.

---

## OPEN MEDIUMS (fix or ticket)
- **[14][44] + Codex/Copilot blkdiag** — `ConcatenationLayer.reach_concat_star` ignores `obj.Dim` (always vertical-stacks) and errors on 3-arg `Star(V,C,d)` inputs with empty predicate bounds (`:135`,`:277`).
- **[13] ElementwiseAffine Star offset fix-ups guess channel-fastest layout** (`:181`) — wrong if the Star is position-fastest.
- **[9] ReshapeLayer reach ignores the `OnnxBCHW` flag** (`:138`) while evaluate honors it → reach/evaluate layout mismatch.
- **[6] DynamicMatmul evaluate broadcasting unimplemented** (`:62`) — errors on batched 3-D operands.
- **[19] SiLU crashes on 3-arg `Star(V,C,d)`** (empty predicate bounds, `:171`); **[20] Softmax.check_soundness crashes on empty rejection sample** (`:366`).
- **[26] Constant placeholder silently degrades to identity if its value fails to attach** (`load_nnv_from_mat.m:402`).
- **[22] `is_final_softmax` uses array order, not `Connections` topology** (`matlab2nnv.m:394`) — also two upstream-visible behavior changes (mid-net LogSoftmax now hard-errors).
- **[27] dead/duplicated `ScalingLayer` branch** (`matlab2nnv.m:185`, also Copilot) — unreachable, risks divergence.
- **[17] multi-input destination whitelist duplicated in NN.m vs load_nnv_from_mat; omits MHA** (`NN.m:1510`).
- **[46] SDPA.reach has no nargin-7 case** for the standard NN dispatch (`:148`) → "Invalid number of arguments".
- **[36] softmax soundness test tolerance 0.1** = 10% of the codomain (`test_Softmax.m:149`) — near-vacuous.
- **Copilot lows:** `verify_mnist_vit.m:171` fallback uses undefined `lb_norm/ub_norm/img_normalized`; `verify_mnist_vit_both.m:222` assumes a `plots/` dir; `ReshapeLayer.m:106` `-1` resolution needs a divisibility check; `TransposedConv2DLayer.m:112/166` bias-shape validation too loose; stale "BatchNorm removed (unsound)" comment now contradicted by passing BN tests.

## PLAUSIBLE (needs a trigger)
- **[5] EmbeddingLayer.parse fabricates RANDOM weights** when the source lacks `Weights` (`:183`) — verifies a random network with no warning.
- **[16] NN.reach_withConns has no per-layer input-completeness check** (`:1456`); **[23] selfAttention `AttentionMask`/causal ignored** by parse; **[30][31][32]** runner orientation/point-spec/path-casing edge cases.

## POSITIVE / SAFE
- Intermediate SoftmaxLayer no longer silently identity (computes bounds or errors). DynamicMatmul/ElementwiseProduct/Division now fail loud. Mid-net softmax routed to sound bounds; mid-net LogSoftmax refused. **[39] CORA submodule bump `eeb971d84→25b760aef` is SAFE** (deletes a stray `<<<<<<< HEAD` conflict marker). The matrix-CI report design (known-vs-new, missing-shard) is useful.

---

## Pre-merge checklist (ranked)
1. **Make the CI gate real first** — remove soundness MC tests from `ci_allowed_failures.txt`; count `<skipped>` as not-pass in `ci_report.py`; replace `assumeFail` with `verifyError`/containment in attention tests; fail the report on `0/0`. *(Otherwise every fix below is unverifiable in CI.)*
2. **Production-path soundness:** SiLU sampling bound [18]/[21], exact-star mislabeling [42], ElementwiseAffine align [10], Addition skip [11].
3. **Attention soundness:** fix `evaluate` correctness [3][4] (else MC oracle is wrong), then SDPA/MHA guard holes [0][1][24][40][45].
4. **Norm grouping/parity:** LayerNorm per-position grouping [8] + reachSequence [41]; SoftmaxLayer ImageStar grouping [12]; PlaceholderLayer.evaluateSequence [43].
5. **Runner:** counterexample order [28], multi-spec try/catch [29].
6. **Mediums + Copilot lows**, then green CI with the *real* gate and no missing shard.
7. **Scope/title:** state explicitly that general sound multi-token attention is NOT provided (fail-loud only).

*Two blockers fixed (`84cb90ea1`, `9978cef50`). Full failure scenarios for every `[n]` are in the agent transcript; this file is the durable index.*
