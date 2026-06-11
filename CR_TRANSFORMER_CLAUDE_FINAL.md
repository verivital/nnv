# PR #290 — FINAL Review & Change Overview (merge decision)

Branch `transformer` → PR #290 (`verivital/nnv`), base `master`@`16d23d457`. This is the
merge-facing summary: the final review verdict, the complete set of changes, the soundness
guarantee and how it is tested, and the residual known items. Per-finding history is in
`CR_TRANSFORMER_CLAUDE.md`; the bound philosophy is in `REMEDIATION_AND_SOUNDNESS_STRATEGY_v01/v02.md`;
benchmark coverage is in `code/nnv/examples/Submission/VNN_COMP2025/VNN_COMP2025_SUPPORT.md`.

## VERDICT: ✅ SAFE TO MERGE

Scope: **sound single-token / FC-simulated transformer utilities + fail-loud multi-token guards +
import-correctness fixes + a real CI gate + an exact-star over-approximation gate.** Every reach is a
sound over-approximation or a fail-loud refusal; **no path certifies not-robust/unsafe from an
over-approximation** (the cardinal soundness rule), verified by 4 review rounds, the layer soundness
harness, the exact-star whitelist auto-test, and a full VNN-COMP sweep. General **sound multi-token
ViT verification is the documented follow-up** (it is fail-loud, hence sound, today).

## CI status (the merge gate)

Last fully-green commit `a6f877643`: **CI (matrix-sharded) ✅, CI (legacy full suite) ✅, Regression
Tests ✅, Deploy Documentation ✅**; PR `mergeable: MERGEABLE, state: CLEAN`. The BN-fix commit
`1833a5c9d` (small, well-tested) is re-running CI; confirm green (`gh pr checks 290 --repo
verivital/nnv`) before clicking merge.

Local final checkpoint: **layer soundness harness SOUND** (SiLU/Gelu/Sigmoid/Tanh/Softmax-mid/
LayerNorm); exact-whitelist 7/7; exact-star-overapprox 7/7; transformer-soundness 17/17; vnncomp25
47/47; acasxu 11/11; matlab2nnv 8/8; verify_robustness/safety 8/8 each; ~170 modified-layer dedicated
tests pass. 0 failures.

## Scale

`16d23d457..HEAD`: ~73 commits (32 this remediation session). Engine: 33 files, +6.2k/-274. Tests:
28 files, +7.5k (≈45 new regression tests across vnncomp25 + soundness/).

## The review process — 4 adversarial rounds

Multi-agent reviews; survivors verified serially against the live MATLAB engine (verifiers must NOT
share the single MATLAB MCP session — they deadlock; finders reason from code, I verify).

- **R1 (full session diff):** 3 confirmed — combiner layers wrongly trusted exact, EAffine zero-pred
  ImageStar scale dropped the channel axis, + the exact-star gate itself.
- **R2 (exact-star gate deep-dive):** `SignLayer.reach` is UNSOUND (returns +1 only) → dropped from
  whitelist; flagged pre-existing out-of-scope NNCS soundness bugs.
- **R3 (find-only, reach-bound math):** Reshape `obj.targetDim` mutation corrupted the MC oracle;
  Reshape ImageStar ignored OnnxBCHW; Conv1D/PixelClassification over-approx dropped from whitelist;
  EAffine sub-dim tiling fail-loud. **REFUTED (verified sound):** SiLU zono, LayerNorm grouping.
- **R4 (final merge-safety):** **BatchNorm plain-Star reach was loose (added a constant to every
  generator) yet whitelisted exact** → fixed to the exact affine + pinned by a new test. Remaining
  candidates were already-documented misuse-only / dormant / coverage gaps (now closed with tests).

Each round found real issues **in this session's own fixes** — the exact-star gate area is subtle, so
it got the scrutiny it needed; the reach-bound *math* was adversarially confirmed sound.

## Complete change overview (grouped)

**Reach soundness — sound bounds:** LayerNorm sound var bound + universal cap `|z|≤√(C-1)` (`84cb90ea1`);
SiLU chord-correction `M·h²/8` error bound + sound box fallback + sound getLinearBounds
(`b408a8e6e`,`5ac745050`); GELU variant-aware (exact erf vs tanh) + sound split bracket (`d73927e91`);
intermediate Softmax per-pixel channel groups + sound estimateRanges scalability (`bd6eb35e4`,`cf0608f76`);
**BatchNorm plain-Star exact affine** (`1833a5c9d`).

**Silent-wrong → fail-loud:** Addition dropped-addend, EAffine bad/ambiguous bias, Embedding random
weights, Constant with no value, Concat non-feature-axis flat-Star, Placeholder `UnsupportedOp:*`,
DynamicMatmul mismatched leads, multi-token attention (SDPA/MHA), Softmax unsupported type
(`ea60f4d41`,`10e5770ad`,`60c9ba1c7`,`9978cef50`,`5a1a77353`,`326b8621c`,`4ebe22d6a`).

**Wrong-but-recoverable → correct:** EAffine `[1,1,C]` H==C axis + row-vs-column outer-product
(`68d95d635`,`a6f877643`); Reshape OnnxBCHW reach (Star + ImageStar) + no-mutate-targetDim
(`60c9ba1c7`,`f8accc670`); MHA evaluate row-softmax + contiguous heads (`edbd7332b`); matlab2nnv
topology-BFS final-softmax + dead-branch removal (`eca9776a2`); dispatch fixes EProduct/SDPA
(`6dbc10ef4`).

**Exact-star over-approximation gate [42] (your Option B):** `NN.exactReach` = exact-star AND every
layer exact (`layer_reach_is_exact` whitelist of verified affine/PWL/inert); all not-robust promotions
gated on it; one `NNV:exactStarOverapprox` warning when degraded (`91e4dd3ee`+, refined `a57d6e56d`,
`5925dad45`, `f8accc670`, `1833a5c9d`).

**CI / runner / hygiene:** real CI gate (skip≠pass, 0/0 fails, soundness un-allow-listed) `11a24df8f`;
runner needReshape==3 witness inverse + multi-spec try/catch `f2c4325f8`; gitignore test figures +
softmax test tolerance `584c988ce`; M2NIST segmentation tutorial tabled in CI (too slow) `cf0608f76`.

## Soundness guarantees — and how they are MACHINE-CHECKED

1. **Every reach ⊇ evaluate (sound over-approximation)** — `tests/soundness/test_layer_soundness_harness.m`
   MC-checks each activation/norm on adversarial boxes; per-layer `test_soundness_*` suites.
2. **Over-approximation never certifies not-robust** — `test_exact_star_overapprox.m` (gate behaviour)
   + `test_exact_whitelist_soundness.m` (every whitelisted layer is MC-sound AND tight; affine layers
   now MC-pinned — this is what catches a future BatchNorm/EAffine-style loosening).
3. **reach == evaluate where claimed exact; fail-loud is loud** — `vnncomp25` Tests 1–47 +
   `transformer_soundness_regression` Tests 1–10 pin each fix with a regression that would have caught
   the original bug.

## VNN-COMP support (full sweep, current code)

`VNN_COMP2025_SUPPORT.md`: across both import paths (Python-importer/manifest and category/matlab2nnv),
**no instance returns sat/not-robust from an over-approximation**; verdicts are `unsat` (SAFE) /
`unknown`; refusals are fail-loud; the runner gates reach on cross-validation. ~21/26 import +
forward-pass-correctly; `malbeware`+`metaroom` proven SAFE; most others are compute-bound timeouts at a
10 s smoke budget (need ≥120 s / cloud); fail-loud: `cctsdb_yolo` (in-graph YOLO post-proc), `vit_2023`
(multi-token attention).

## Residual known items — none are live unsoundness

- **Sound multi-token attention** (ViT) — fail-loud today (sound); the capability follow-up (plan v02
  §B1: value-hull → softmax-aware; `evaluate` is now correct so MC-validation is meaningful).
- **Low / latent / misuse-only (documented):** SDPA `ValueDim==0` / `parse` no-fail-loud on missing
  ValueSize (real nets go through MHA, which is armed); BatchNorm/Conv `*_Sequence` box branches
  (dormant — not on the gated reach path); LayerNorm sequence cap assumes group==NumChannels (sound for
  standard channel-only LayerNorm); stale `test_SoftmaxLayer_reach_soundness.m` diagnostic.
- **Pre-existing, out-of-scope (flagged for maintainer, separate subsystems):** `SignLayer.reach`
  unsound; `LinearNNCS`/`DLinearNNCS` `.verify` method-string gate + `.falsify` discards `U.contains`;
  `verify_safety` parfor undefined `method`. See `CR_TRANSFORMER_CLAUDE.md` §"Additional PRE-EXISTING".
- **[42-exact] maintainer note:** the gate downgrades not-robust→unknown for over-approx layers; if you
  prefer hard-rejecting `exact-star` on those layers instead, that is a one-line policy change.

## Pre-merge checklist

1. Confirm CI green on the latest commit (`gh pr checks 290 --repo verivital/nnv`).
2. (Optional) decide the [42] policy (downgrade-to-unknown, current) vs reject-exact-star.
3. (Optional) repo hygiene: drop scratch `TODO_*`, gitignore benchmark blobs, copy the Python importer
   into `code/nnv/tools/` for a self-contained PR.
4. Merge — the soundness scope above holds; sound multi-token ViT is the tracked follow-up.
