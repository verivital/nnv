# CI Test Tracking — failing & excluded tests

Single source of truth for tests that **fail** (kept visible in CI, to be resolved) or are
**excluded** from the standard matrix CI (run elsewhere). Update this as tests are fixed.

Machine-readable encodings (the actual gates):
- **Allowed (failing-but-known) tests:** [`code/nnv/tests/ci_allowed_failures.txt`](code/nnv/tests/ci_allowed_failures.txt) — failures here are visible but do **not** red the build.
- **Excluded tests:** the `gpuPatterns` / `anticipatedHeavyPatterns` / `confirmedCrashPatterns` lists in [`code/nnv/tests/run_shard.m`](code/nnv/tests/run_shard.m).
- The **report job** ([`ci_report.py`](code/nnv/tests/ci_report.py)) prints these every run and reds the build only on **NEW** failures or crashes.

_Last updated 2026-06-08 (post CI overhaul, commit 36a1483a)._

---

## A. FAILURES to RESOLVE (run in CI, classified "known/allowed", kept visible)

### A1. SLM / transformer soundness — **resolve as part of finishing transformer support**
**Why:** these tests share variables (`results`, `N_SAMPLES`) across `%%` sections, but
`runtests` runs each `%%` in a fresh workspace → they pass when run as a *script* but error
as unit tests (`Unrecognized function or variable 'results'`). Plus genuine in-progress
transformer verification. **Fix:** make each `%%` section self-contained (re-init its own
state) — or convert to a function-based test with a shared fixture — then remove from
`ci_allowed_failures.txt`. Tie to the transformer/ViT soundness work.

- [ ] `test_SLM_layers_soundness/test_____1_GeluLayer_____`
- [ ] `test_SLM_layers_soundness/test_____2_ScaledDotProductAttentionLayer_single_tokenVectorInterpretation______`
- [ ] `test_SLM_layers_soundness/test_____3_MultiHeadAttentionLayer_self_attention_Single_token______`
- [ ] `test_SLM_layers_soundness/test_____4_RotaryPositionalEmbedding_deterministicRotationPerPosition______`
- [ ] `test_SLM_layers_soundness/test_____5_SinusoidalPositionalEncoding_additive_Deterministic______`
- [ ] `test_SLM_layers_soundness/test_____6_LearnedPositionalEmbedding_additive______`
- [ ] `test_SLM_layers_soundness/test_____7_EmbeddingLayer_epsilonPerturbationAroundAToken_sEmbedding______`
- [ ] `test_SLM_layers_soundness/test_____SummaryTable_____`
- [ ] `test_SoftmaxLayer_reach_soundness/TestA_Mid_networkSoftmax_Unsound_ShowReach__InputPassthrough_`
- [ ] `test_SoftmaxLayer_reach_soundness/TestB_Last_layerSoftmax_Argmax_preservingCheck`
- [ ] `test_SoftmaxLayer_reach_soundness/SaveResults`

(Also listed for safety in the allowed file: `test_SLM_layers_soundness/Helpers`,
`.../MC_containment...`, `test_SoftmaxLayer_reach_soundness/SoundnessMC_containment...` —
these likely pass; remove if confirmed.)

### A2. Tutorial integration (`test_all_tutorial`) — Linux/CI environment fixes
**Why:** pass locally (Windows) but fail on the Linux runner. Concrete errors:
- `test11_NNCS_AEBS_reach_`, `test13_NNCS_InvertedPendulum`: *"Unable to change current folder
  to .../ACAS Xu/../../AEBS (nonexistent)"* — fragile **relative `cd`** between tutorial folders
  on the case-sensitive Linux FS.
- `test7_NN_MNIST_verify_`: *"Unrecognized field name 'net'"* — a loaded `.mat` lacks the
  expected `net` field on CI.
- `test15_Other_SetRepresentations`: errors (same family).

**Fix:** use absolute paths via `nnvroot()` instead of relative `cd`; verify the `.mat`
variable name. Not a code/R2026a bug.

- [ ] `test_all_tutorial/test11_NNCS_AEBS_reach_`
- [ ] `test_all_tutorial/test13_NNCS_InvertedPendulum`
- [ ] `test_all_tutorial/test15_Other_SetRepresentations`
- [ ] `test_all_tutorial/test7_NN_MNIST_verify_`

---

## B. EXCLUDED from standard CI (not failures — run elsewhere)

### B1. Heavy / memory-flaky → high-memory fallback
**Why:** empirically **flake on the 16 GB GitHub runner** — nondeterministic exit-255 crashes
*and* hangs under memory pressure. **Run via** `run_shard(k,N,'highmem')` on a big-RAM /
self-hosted runner (`NNV_HIGHMEM_RUNNER`), ACCRE SLURM, or locally (see `CLOUD_COMPUTE_OPTIONS.md`).
**Recover into standard later** by pinpointing the *few* truly-flaky sections (most run fine)
and narrowing `anticipatedHeavyPatterns` / using `confirmedCrashPatterns`.

62 sections across 13 files (pattern-matched in `run_shard.m`):
`test_VolumeStar`, `test_soundness_VolumeStar`, `test_Conv3DLayer`, `test_soundness_Conv3DLayer`,
`test_AveragePooling3DLayer`, `test_soundness_AveragePooling3DLayer`, `test_PixelClassificationLayer`,
`test_soundness_PixelClassificationLayer`, `test_regression_segmentation`,
`test_regression_verify_segmentation`, `test_soundness_comprehensive_network`, and the
heavy sections of `test_all_tutorial` (segmentation, etc.).

- [ ] (optional) Pinpoint the actual flaky sections (vs the ~45/60 that ran fine) to recover most into standard.

### B2. GPU tests → need a GPU runner
**Why:** no GPU on GitHub-hosted runners (error/incomplete). **Run on** a GPU-equipped machine.
19 sections: the `toGPU` section of ~18 layer tests + `test_gpu_VolumeStar`.

- [ ] (optional) add a GPU self-hosted runner lane if GPU verification needs CI coverage.

---

## D. CORA nonlinear reachability hangs on the Linux runner (class excluded)
**`test_all_tutorial/test13_NNCS_InvertedPendulum`** was alone in its shard, which hung →
definitive. The whole **CORA *nonlinear* family** hangs/runs-away on the Linux runner (vs
~minutes on Windows); excluding one at a time was whack-a-mole as bin-packing reshuffled, so
the class is excluded together (`confirmedCrashPatterns` in `run_shard.m`) and run via the
fallback/locally — **38 elements across:** `test_NonLinearODE_*` (evaluate / reach_zono /
stepReachStar), `test_hybridA`, `test_soundness_HybridA`, `test_soundness_NonlinearNNCS`,
`test_soundness_DNonlinearNNCS`, `test_regression_neural_ode`, and
`test_all_tutorial/{test10_NNCS_ACC, test13_NNCS_InvertedPendulum}`. Linear reachability
(`DLinearODE`, `LinearNNCS`, `DLinearNNCS`) is fine and stays.
- [ ] **ROOT FIX:** investigate why CORA `NonLinearODE` reachability hangs on Linux (solver
  tolerance/iterations vs Windows) + fix the tutorials' relative `cd` (use `nnvroot()`); then
  remove from `confirmedCrashPatterns`. Until then they run via the high-mem fallback / locally.
- Safety net: the `report` job reds the build on any **missing shard**, so a new hang can't
  slip through unseen (and `continue-on-error` lets a cleanly-crashed shard still upload the
  progress trace that names the offender).

## C. How to update
- **Fixed a failing test?** Run it under `runtests` (not as a script) to confirm it passes, then
  remove its line from `ci_allowed_failures.txt` and check it off here.
- **Recovered a heavy/GPU test for standard CI?** Narrow the pattern in `run_shard.m` (or move a
  specific Name out), confirm a green CI run, then update B1/B2 here.
- **A new failure appears in CI** (report shows it as "NEW"): triage — if it's known/expected,
  add it to `ci_allowed_failures.txt` and list it here; otherwise it's a real regression to fix.
- **Regenerate `test_durations.csv`** after adding many tests so sharding stays balanced.
