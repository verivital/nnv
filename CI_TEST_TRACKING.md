# CI Test Tracking — failing & excluded tests

Single source of truth for tests that **fail** (kept visible in CI, to be resolved) or are
**excluded** from the standard matrix CI (run elsewhere). Update this as tests are fixed.

Machine-readable encodings (the actual gates):
- **Allowed (failing-but-known) tests:** [`code/nnv/tests/ci_allowed_failures.txt`](code/nnv/tests/ci_allowed_failures.txt) — failures here are visible but do **not** red the build.
- **Excluded tests:** the `gpuPatterns` / `anticipatedHeavyPatterns` / `confirmedCrashPatterns` lists in [`code/nnv/tests/run_shard.m`](code/nnv/tests/run_shard.m).
- The **report job** ([`ci_report.py`](code/nnv/tests/ci_report.py)) prints these every run and reds the build only on **NEW** failures or crashes.

_Last updated 2026-06-08 (CI green: 947/949 pass; Linux-duration-balanced shards + infra-flake
retries on setup-matlab & run-command; raw-logs-API hang diagnosis — see §E)._

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

### B3. Exit-255 crashers on the 16 GB runner → high-memory fallback
**Why:** these crash MATLAB outright (exit 255, no clean error) on the standard runner — the
report's progress trace named them as the last `START` without a `DONE`. Likely memory
pressure / a solver crash on CI, not a code bug (they pass locally). Excluded via
`confirmedCrashPatterns`; run via `highmem`/locally.
- [ ] `test_conv_layer_perturbation/test_conv_layer_perturbation` (~154s; heavy conv perturbation)
- [ ] `test_all_tutorial/test8_NN_MNIST_verify_fc_` (crashes mid-tutorial on CI; passes locally)

---

## D. CORA nonlinear reachability hangs on the Linux runner (class excluded)
**`test_all_tutorial/test13_NNCS_InvertedPendulum`** was alone in its shard, which hung →
definitive. The whole **CORA *nonlinear* family** hangs/runs-away on the Linux runner (vs
~minutes on Windows); excluding one at a time was whack-a-mole as bin-packing reshuffled, so
the class is excluded together (`confirmedCrashPatterns` in `run_shard.m`) and run via the
fallback/locally — **≈40 elements across:** `test_NonLinearODE_*` (evaluate / reach_zono /
stepReachStar), `test_hybridA`, `test_soundness_HybridA`, `test_soundness_NonlinearNNCS`,
`test_soundness_DNonlinearNNCS`, `test_regression_neural_ode`, and **all three nonlinear-NNCS
tutorials** `test_all_tutorial/{test10_NNCS_ACC, test11_NNCS_AEBS, test13_NNCS_InvertedPendulum}`
(each was isolated alone in a shard that then hung — definitive). Linear reachability
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
  Best source is a green CI run's JUnit (actual Linux per-test seconds), not local timing:
  parse every `results-*.xml` `<testcase time>` and merge over the existing CSV (keep local
  for any test that didn't run on CI). Local Windows durations under-estimate Linux time and
  cause uneven shards.

---

## E. CI infrastructure flakiness & diagnosis (operational)

The standard matrix is GREEN for the suite (run 27176965670: 14/14 shards, 947/949 pass — only
the §A2 tutorial errors remain, classified known/allowed). The remaining red cause is **GitHub
free-runner infra flakiness**, NOT NNV:

- **What flakes:** `matlab-actions/setup-matlab@v2` (downloads MATLAB + 11 toolboxes, ~5 min)
  occasionally fails on a transient mirror/network error; more rarely MATLAB exits 255 at
  launch. With ~11 shards, P(some shard flakes) ≈ 50%/run. The failing shard **moves
  run-to-run** (observed: 7, then 5) — the tell-tale of infra vs a deterministic test.
- **Mitigation (in `test-matrix.yml`):** both the **Setup MATLAB** and **Run shard** steps
  retry once (`id` + `continue-on-error: true` on the first; a duplicate step gated on
  `steps.<id>.outcome == 'failure'`). This drops a shard's double-failure to ~0.5% → ~95% of
  runs green on the first attempt. Backstop for a persistently-bad runner: `gh run rerun --failed`.
- **Why `timeout-minutes: 15` + `continue-on-error` on run-command:** a hung shard's step is
  killed but the job still uploads its progress trace, so the report can name the hanger.

### Diagnosing a red run
1. **Report job** (auto, in the run summary): pass/fail counts, NEW vs known/allowed failures,
   **crashes** (named via the progress trace), and **missing shards**.
2. **A failed/missing shard → which step?**
   `gh run view <RID> --repo ttj/nnv --json jobs --jq '.jobs[]|select(.conclusion=="failure")|{name,steps:[.steps[]|{name,conclusion}]}'`
   - failed at *Setup MATLAB* (and its retry) = infra → `gh run rerun --failed`.
   - failed at *Run shard* = a test crashed/hung → name it (next).
3. **Name a crash/hang in a killed shard** (the artifact upload is skipped on a hard kill, so the
   report can't see it) via the **raw-logs API** — the per-job log endpoint returns the full
   run-command stdout even for cancelled jobs (the `gh run view --log` subcommand returns only
   teardown):
   ```bash
   jid=$(gh run view <RID> --repo ttj/nnv --json jobs \
        --jq '.jobs[]|select(.name=="tests <K>/<N>")|.databaseId')
   MSYS_NO_PATHCONV=1 gh api "/repos/ttj/nnv/actions/jobs/$jid/logs" \
        | sed 's/\r$//' | grep -aE "Running test_|code 255" | tail
   ```
   The last `Running test_…` before the kill = the offender. If deterministic, add it to
   `confirmedCrashPatterns` in `run_shard.m`; if it only flakes, leave it (infra).

### Caching the MATLAB install (the real fix for setup cost + flakiness)
`setup-matlab@v3` with `cache: true` (on all 3 setup steps) stores MATLAB + the 11 toolboxes
in the GitHub Actions cache. Measured on this repo (2026-06-08):
- **Cold run** (cache miss — e.g. first run, or after changing `release`/`products`): downloads
  from MathWorks **and** saves the cache. Wall ~1134s (~19 min), shard mean 787s.
- **Warm run** (cache hit): restores the install in **~55s** instead of the ~5 min download.
  Wall **~611s (~10 min)**, shard mean 382s — **~46% faster**, and with no MathWorks round-trip
  the download flake is gone on warm runs. Cache key:
  `matlab-cache-linux-x64-<MATLABversion>-<products-hash>`; it **fits under GitHub's 10 GB
  limit** (it saved successfully). The key auto-changes (→ a one-time cold rebuild) whenever you
  edit `release` or the `products` list.
- Shard count auto-scales to **~11 (cap 16)**: over-sharding to ~18 REGRESSED (more concurrent
  jobs drew GitHub's slow-runner tail — per-test time ~4× worse, wall 1349s vs warm-N11 611s),
  so ~11 is the empirical sweet spot (~10 min warm). `run_shard` packs by (per-test overhead
  ~3s + recorded JUnit time) since fixed fixture/JIT overhead dominates the sub-second tests.

### Rebuilding the cache on demand
The cache auto-rebuilds when `release`/`products` change. To FORCE a rebuild *without* a config
change (cache corruption, a MATLAB point-release that reuses the version string, periodic
refresh):
- **One-click (recommended):** run the workflow via `workflow_dispatch` with **`rebuild_cache:
  true`** — the `prepare` job deletes every `matlab-cache-*` entry, then the matrix run
  repopulates them:
  `gh workflow run test-matrix.yml --repo ttj/nnv --ref transformer -f rebuild_cache=true`
  (or the "Run workflow" button in the Actions UI, with "Force-rebuild…" checked).
- **Manual (no run):** `gh cache list --repo ttj/nnv` → `gh cache delete <id> --repo ttj/nnv`
  (needs `actions: write`); the next run then repopulates.
- The first run after any rebuild is a cold run; the setup + run-command retries cover a
  transient download flake during it.

### Alternatives considered (and why cache:true won)
- **`cache: true` (CHOSEN):** one-line, keeps the free public-repo auto-licensing, ~46% faster
  warm runs, fits the 10 GB cache. Best effort/payoff.
- **Prebuilt Docker image (MATLAB + toolboxes baked in):** fastest possible startup, but the
  free public-repo auto-license is a `setup-matlab` feature — a bare `mathworks/matlab` image
  needs its own batch-token/network license (the friction point from `CLOUD_COMPUTE_OPTIONS.md`)
  plus image build/maintenance. Reserve as the fallback **only if** the cache ever stops fitting
  the 10 GB limit (e.g. many more toolboxes).
- **Self-hosted runner with MATLAB preinstalled:** zero per-run setup, but leaves the free
  GitHub-hosted pool and needs campus-license/runner upkeep — that's the high-mem fallback lane,
  not the default.
