# Merge-Readiness Plan — NNV `transformer` branch PR

*2026-06-09. Author: Claude (Opus 4.8). Provenance: a multi-agent workflow synthesis (its 5 parallel
investigators were transiently rate-limited, so the per-area deep-dives are from the synthesis + this
session's direct code reads + **MCP verification on R2026a**). Every root cause cites a verified
`file:line`. Companion docs: [`code/nnv/tests/soundness/TRANSFORMER_SOUNDNESS_AUDIT.md`](code/nnv/tests/soundness/TRANSFORMER_SOUNDNESS_AUDIT.md), [`CI_TEST_TRACKING.md`](CI_TEST_TRACKING.md).*

**Hard constraint:** the branch must not merge with any *silently unsound* code path. Soundness here =
the computed reachable set provably contains every actual output, verified empirically by Monte-Carlo
(MC) containment (these are over-approximations, not formal proofs).

**Verdict:** The branch is **mergeable under a fail-loud-gate strategy** after **one small required code
change** (the `EmbedDim` fix below, which closes B3+B4 together). Every known unsound path is then either
(a) fixed and MC-verified, or (b) guarded to *raise* instead of silently producing a wrong set. The one
genuinely broken end-to-end path (real-attention ViT) is fenced off by guards + documentation, not
shipped as a working-but-wrong feature.

---

## 1) MERGE BLOCKERS (fix-or-gate before merge)

Legend: **MC-containment** = N random samples of the concrete layer, assert every sample ∈ computed
`[lb,ub]`. **fail-loud** = the unsound branch now `error()`s, and a regression test asserts it errors.

### B1 — Intermediate (attention) softmax returned an unsound identity passthrough — **FIXED**
- **Root cause:** `SoftmaxLayer.reach` routed non-final softmax through `Softmax.reach_star_approx(IS, method='')`; the empty method errored and a `try/catch` silently returned the identity (input-domain) set. `engine/nn/layers/SoftmaxLayer.m`.
- **Status:** FIXED (`d5f3aa05e`) — now calls `Softmax.reach_star_approx_bounds(IS)` directly; unsupported set types raise. **MC-verified:** `final=false` ∈ [0,1], 0/3000 violations. Locked by Tests 1 & 3 in `tests/soundness/test_transformer_soundness_regression.m` (not in `ci_allowed_failures.txt` → gates CI).

### B2 — ONNX loader marked mid-network softmax as final — **FIXED**
- **Root cause:** `load_nnv_from_mat` built `SoftmaxLayer(name)` with default `IsFinalLayer=true`. `engine/utils/load_nnv_from_mat.m`.
- **Status:** FIXED (`d5f3aa05e`) — position-aware post-pass sets `IsFinalLayer=(i==n)`. Erring toward non-final is sound (only looser). All 28 vnncomp25 regression tests pass.

### B3 + B4 — Multi-token / real-attention reach (UNIFIED by one root cause) — **needs the EmbedDim fix**
**These two were diagnosed separately by the synthesis; MCP shows they are the *same* bug.**

- **MCP-VERIFIED root cause:** `matlab2nnv` on `mnist_vit_attention_model.mat` (a `DAGNetwork`) **extracts
  `W_Q` fine** ([512×512], non-empty, no warning) — weight extraction is **not** the problem. The bug is
  **`EmbedDim=0`, `HeadDim=0`**: `MultiHeadAttentionLayer.parse` (≈lines 623-627) reads `EmbedDim` from
  `NumChannels`/`NumKeyValueChannels`, but R2026a `selfAttentionLayer` exposes
  `InputSize`/`NumQueryChannels`/`NumKeyChannels`/`NumValueChannels`/`OutputSize` — so neither `isprop`
  matches and `EmbedDim` stays 0. Consequences:
  1. In `compute_mha_bounds`, the per-head loop uses `idx_start=(h-1)*HeadDim+1`, `idx_end=h*HeadDim`
     with `HeadDim=0` → range `1:0` = **empty** → empty head bounds → **empty ImageStar** → LayerNorm
     chokes → garbage bounds (class lb −134.9, max-other ub 155.9 @ ε=0.02) → 0/5 "verified"
     (`verify_mnist_vit_attention.m`).  *(This is "B4".)*
  2. The multi-token guard's `obj.EmbedDim > 0` test is false (EmbedDim=0), so it **doesn't fire** — a
     genuine multi-token ViT slips past as "single-token."  *(This is the "B3 guard hole.")*
- **Underlying unsoundness (B3):** even with dims correct, `compute_attention_bounds` returns single-token
  V-bounds (`ScaledDotProductAttentionLayer.m`, `MultiHeadAttentionLayer.m:compute_attention_bounds_single_head`),
  which is unsound for seq_len>1 (the output is a cross-token convex combination that can leave any single
  token's range). Guards added (`27311daf1`) raise `multiTokenUnsound` — but only when `EmbedDim>0`.

- **REQUIRED fix (small, ~3 lines, closes B3 + B4 together):** in `MultiHeadAttentionLayer.parse`, derive
  dims from the **extracted weights** (version-independent), after the weight block:
  ```matlab
  if isempty(L.EmbedDim) || L.EmbedDim == 0
      if ~isempty(L.W_Q), L.EmbedDim = size(L.W_Q, 1); end
  end
  if L.NumHeads > 0 && L.EmbedDim > 0, L.HeadDim = L.EmbedDim / L.NumHeads; end
  ```
  Keep the `NumQueryChannels`/`InputSize` property reads as a fallback. Apply the same `size(W_*,1)`
  derivation in `ScaledDotProductAttentionLayer` (`ValueDim`).
- **Effect:** `HeadDim>0` → no more empty sets; for a real (multi-token) ViT the input is `seq_len·EmbedDim`,
  so the guard **fires and fails loud** (sound refusal) instead of returning garbage. Single-token models
  compute bounds as before.
- **Resolution:** **fix the dim derivation (sound) + keep the multi-token guard (fail-loud).** The full
  *sound multi-token bound* is deferred to post-merge F1.
- **Effort:** S. **Risk:** low (single-token unit tests unaffected). **Soundness verified:** add a
  regression case — a real-ViT-shaped (multi-token) MHA must now `error('multiTokenUnsound')`, not return a
  verdict; the 19 SDPA + 24 MHA single-token tests must still pass.
- **Also mark** `verify_mnist_vit_attention.m` with an EXPERIMENTAL header ("known-broken pending F1/F2;
  use `verify_mnist_vit.m`"), so it isn't presented as a working demo.

### B5 — `matlab2nnv` softmax → PlaceholderLayer (identity), unsound for mid-network MATLAB softmax — **gate**
- **Root cause:** MATLAB import drops `SoftmaxLayer` to `PlaceholderLayer` (identity). `engine/utils/matlab2nnv.m:45-48`. Sound for a *final* argmax-on-logits spec; unsound for a *standalone mid-network* softmax (rare — MATLAB self-attention is self-contained in `selfAttentionLayer→MHA`).
- **Resolution:** **fail-loud-gate.** Mirror B2: when the dropped `SoftmaxLayer` is **not** the last layer, route it to the sound non-final `SoftmaxLayer` path or **raise** — never silently identity-pass a non-final softmax. **Effort:** S. **Verify:** assert a 2-softmax MATLAB net yields no non-final identity placeholder.

**After B3+B4 (EmbedDim fix) and B5 land, no silently-unsound path remains → hard constraint satisfied.**

---

## 2) VNN-COMP 2025 importer — finish vs defer

**Soundness-blocking (must be green to merge):**
- Tests 15/16 (intermediate softmax ∈ [0,1]; final softmax identity) — importer-side mirror of B1/B2; confirm they gate (not in `ci_allowed_failures.txt`). ✓ already.
- **Add MC-containment for the dynamic layers** `ElementwiseProductLayer`, `ElementwiseDivisionLayer`,
  `DynamicMatmulLayer` (existing Tests 17+ check shapes/values, not *reach* soundness). **Division is the
  trap:** interval division across zero is a classic unsoundness — assert the layer bounds correctly or
  **rejects** denominators straddling 0. This is the one importer soundness gap to close before merge.
- **Audit `PlaceholderLayer` routing for unknown ops** — an *active* (non-identity) op silently dropped to
  identity is the same unsoundness class as B5. Unknown ops must **error at import**, never load a partial net.

**Coverage-incomplete (defer; document as "partial VNN-COMP'25 support"):**
- The 3 categories blocked by MATLAB's `importNetworkFromONNX` (ONNX IR/opset): `cctsdb_yolo` (gather
  forward-prop), `lsnc_relu`, `traffic_signs_recognition`. The Python-importer path
  (`onnx2nnv.py → load_nnv_from_mat`) is the route to finish them, but it needs the model files + iteration
  — **post-merge (F4)**. `collins aerospace` produced "invalid SAT instances" — flag as a **correctness**
  concern to investigate, not just coverage.
- 26/29 categories load; `test_vnncomp25_regression` is 28/28. Benchmark breadth is a follow-up.

**Rule:** importer merges once the dynamic-layer MC tests are green. Breadth is post-merge.

---

## 3) Untracked files + CI gating

**Commit (part of the PR):** all `tests/soundness/*.m`, `tests/vnncomp25/*.m`, `tests/nn/**/*.m`; the CI
harness (`tests/run_shard.m`, `ProgressPlugin.m`, `ci_report.py`, `ci_allowed_failures.txt`,
`test_durations.csv`); the workflow (`.github/workflows/test-matrix.yml`); the docs
(`CI_TEST_TRACKING.md`, `TRANSFORMER_SOUNDNESS_AUDIT.md`, this plan).

**Gitignore / do-not-commit:** `code/nnv/test-results-ci/` (already ignored); the VNN-COMP
`examples/Submission/VNN_COMP2025/+<benchmark>/` data dirs and stray big binaries
(`ARCH-COMP2023/.../acc.mat`, `acc.pdf`); the multi-MB example model blobs
(`mnist_vit_*.mat/onnx`, `sentiment_*.mat/onnx`) — prefer a `download_models.m` + checksum (the
`SST2/download_sst2.m` pattern) over committing ~30 MB to an upstream PR. If kept, call them out in the PR
description.

**Drop:** the `TODO_*.md` scratch family (`TODO_SLM-pass01..15`, `-orig - Copy`, `TODO_TRANSFORMER-pass*`,
…) — local notes, not artifacts. And the untracked **demonstration** tests
`tests/soundness/{test_SLM_layers_soundness, test_SoftmaxLayer_reach_soundness}.m` — they *record*
unsoundness rather than pass/fail and break under `runtests` (cross-`%%` scoping). They're superseded for
CI by `test_transformer_soundness_regression.m`. Either convert to self-contained pass/fail or leave out;
until then keep them in `ci_allowed_failures.txt` (they are) but do **not** advertise them as the guarantee.

**Checks that become REQUIRED (branch protection on the matrix CI):**
- `test_transformer_soundness_regression` — REQUIRED, 7/7. The hard-constraint enforcer.
- `test_vnncomp25_regression` — REQUIRED (softmax Tests 15/16 + the new dynamic-layer MC tests).
- A shard **CRASH or MISSING shard** reds the build (`ci_report.py` already does this) — REQUIRED.
- **Documented review rule: no soundness test may ever be added to `ci_allowed_failures.txt`** — that
  would defeat the gate. The allow-list is only for env-flaky tutorials and explicitly-WIP SLM demos.
- `run_shard.m` auto-discovers `soundness/` + `vnncomp25/` via `TestSuite.fromFolder(...,'IncludingSubfolders',true)` — the soundness tests ride the matrix with no per-test wiring.

---

## 4) Post-merge follow-ups (issues; none block merge)

- **F1 (HIGH) — Sound multi-token attention bound.** Replace the B3 fail-loud refusal with a real bound.
  Cheapest sound version (valid because attention weights live in the simplex): per output coordinate,
  `out_(i,d) ∈ [min_j V_lb(j,d), max_j V_ub(j,d)]` (value-hull bounding box, ignores Q/K). Requires the
  reach API to preserve `[seq_len,dim]` instead of flattening. Verify by MC for seq_len∈{2,4,8}. Refs:
  **Shi 2020** (arXiv:2002.06622, closed-form linear attention bounds), **Wei 2023** (arXiv:2303.01713,
  softmax convex bounds), DeepT/Bonaert 2021, Vertex-Softmax 2026 (arXiv:2605.10974).
- **F2 (HIGH) — Real-attention ViT pipeline.** With F1 + the EmbedDim fix, make `verify_mnist_vit_attention.m`
  produce *sound* verdicts (then un-mark EXPERIMENTAL); revisit the `single`-precision conversion warning.
- **F3 (MED) — Tighter softmax bounds (Wei 2023 linear).** Current `Softmax.compute_softmax_bounds` is sound
  but loose (interval box → many UNKNOWN). Implement `compute_softmax_bounds_tight` (currently a sampling
  placeholder). Tightening is sound-preserving; verify it still MC-contains and strictly narrows the box.
- **F4 (LOW) — VNN-COMP'25 breadth** (the 3 blocked categories via the Python importer; per-benchmark tests).
- **F5 (LOW) — Convert/retire the SLM demonstration soundness tests.**

---

## 5) Sequenced roadmap (with verification commands)

One-time per MCP session: `cd code/nnv; startup_nnv; NNVVERSION()` (expect `NNV v3.0.0`).

0. **Baseline (confirm).** `runtests('tests/soundness/test_transformer_soundness_regression.m')` → 7/7;
   `runtests('tests/vnncomp25/test_vnncomp25_regression.m')` → green. If either fails, stop.
1. **The one REQUIRED code change — EmbedDim fix (closes B3+B4).** Edit `MultiHeadAttentionLayer.parse`
   (and `ScaledDotProductAttentionLayer`) to derive `EmbedDim`/`HeadDim` from `size(W_*,1)`. Re-run:
   `test_transformer_soundness_regression` (7/7 + new multi-token-errors case),
   `test_MultiHeadAttentionLayer` (24), `test_ScaledDotProductAttentionLayer` (19). Checkpoint:
   multi-token (incl. real-ViT shapes) errors; single-token unaffected.
2. **Fail-loud B5 + mark B4 experimental.** `matlab2nnv.m:45-48` don't placeholder-identity a non-final
   softmax; header-banner `verify_mnist_vit_attention.m`. `check_matlab_code` clean; assert no verdict from
   garbage bounds.
3. **Importer dynamic-layer MC tests (§2).** Add MC-containment for ElementwiseProduct/Division/DynamicMatmul;
   division bounds-or-rejects across zero. `test_vnncomp25_regression` green.
4. **Regression safety net.** `verify_mnist_vit.m` → 4/5 @ ε=0.5, 0 errors; `runtests('tests/soundness','IncludeSubfolders',true)` green.
5. **Repo hygiene (§3).** Drop scratch `.md`; gitignore benchmark data + (optionally) model blobs; verify
   `git diff --stat master...transformer` is a clean, reviewable surface; no soundness test in the allow-list.
6. **CI gate + branch protection.** Push; matrix green with the REQUIRED checks; enable protection → **MERGE.**

---

## What soundness holds at merge (honest statement)
- **Sound, MC-verified (0 violations):** intermediate/final softmax; BatchNorm; Gelu/LayerNorm/positional
  encodings (exact-affine); the FC-simulated-attention ViT pipeline (`verify_mnist_vit.m`, 4/5).
- **Sound by refusal (fail-loud, regression-locked):** multi-token SDPA/MHA reach, real-attention ViT,
  mid-network MATLAB softmax import — these **cannot produce a wrong verdict; they error**.
- **Not yet supported (F1/F2):** trustworthy multi-token / real-attention verification — the branch does not
  claim it.
- **Caveat:** soundness claims are **empirical (MC-containment), not formal proofs.** The value-hull
  multi-token bound (F1) is the priority for a rigorous follow-up.
