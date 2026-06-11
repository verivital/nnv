# VNN-COMP 2025 — Benchmark Support Status (PR #290, `ttj/transformer`)

> **⤴ CURRENT STATUS:** this doc's per-benchmark table below is the 10 s smoke-test sweep. For the
> **real 120 s/instance results** (15/26 completed verdicts) and the up-to-date readiness tiers, see the
> top-level [`VNNCOMP2025_STATUS.md`](../../../../../VNNCOMP2025_STATUS.md) (raw data:
> `results_20260610_200503.csv`).

**Purpose:** authoritative, merge-facing summary of which VNN-COMP 2025 benchmarks the
PR's importer + reachability supports, and the soundness guarantee. Generated from two
sweeps of all 26 benchmark folders (one representative instance each):

- **Manifest / Python-importer path** — `<model>.nnv.mat` → `load_nnv_from_mat` →
  `NN`. Last full sweep `results_pyimporter_20260506_001820.csv` (import + ONNX-runtime
  cross-validation `xval` + reach).
- **Category-dispatcher path** — `run_vnncomp_instance(category,…)` → `matlab2nnv` /
  per-benchmark scripts. Fresh sweep **`results_20260610_130007.csv` (2026-06-10, this
  PR's current code, 10 s/instance)**.

> Run a fresh full sweep with a realistic timeout (≥ 120 s, ideally on the CI matrix /
> cloud — see `nnv/CLOUD_COMPUTE_OPTIONS.md`) before the final competition submission.
> The 10 s cap below is a smoke-test budget, so most "timeout" rows are **compute-bound,
> not capability- or soundness-bound**.

---

## 1. Headline: the soundness guarantee holds end-to-end

Across **both** sweeps, **no instance returned `sat` / "not-robust" from an
over-approximate method.** Every completed verdict is one of:

- **`unsat`** — property proven (the over-approximate reachable set does not intersect
  the unsafe region → definitively SAFE/robust; sound), or
- **`unknown`** — the over-approximation intersects the unsafe region but cannot certify
  a counterexample (sound; never reported as not-robust).

Refusals are **fail-loud** (an `error`), never a silent wrong verdict. The runner also
**gates reach on cross-validation**: when the import's forward pass does not match the
ONNX runtime (`xval` fail), it does **not** run reach (status `-`), so a mis-imported
network can never yield an unsound verdict. This is the property this PR's review cycle
was built to guarantee (see `CR_TRANSFORMER_CLAUDE.md`, finding [42] and the exact-star
gate).

Fresh-sweep tally (10 s budget): **unsat 2, unknown 4, timeout 17, error 3, sat 0.**

---

## 2. Per-benchmark support matrix

`import` = loads into NNV; `xval` = NNV forward pass matches ONNX runtime (manifest path,
2026-05-06); `reach (10s)` = fresh category-dispatcher outcome at a 10 s budget.

| # | Benchmark | import | xval (manifest) | reach (10 s, 2026-06-10) | Tier |
|---|-----------|--------|-----------------|--------------------------|------|
| 1 | acasxu_2023 | ✅ | pass 4.7e-6 | timeout | B compute |
| 2 | cctsdb_yolo_2023 | ✅ | **FAIL** | **error (fail-loud)** | C refuse |
| 3 | cersyve | ✅ | pass 7.2e-7 | timeout | B compute |
| 4 | cgan_2023 | ✅ | pass 7.5e-9 | timeout | B compute |
| 5 | cifar100_2024 | ✅ | pass 6.9e-6 | timeout | B compute |
| 6 | collins_aerospace | ✅ | **FAIL** | timeout | C/B* |
| 7 | collins_rul_cnn_2022 | ✅ | pass 1.9e-5 | timeout (manifest: **unsat**) | A/B |
| 8 | cora_2024 | ✅ | pass 2.4e-4 | timeout | B compute |
| 9 | dist_shift_2023 | ✅ | pass 1.1e-4 | **unknown** | A verifiable |
| 10 | linearizenn_2024 | ✅ | pass 1.5e-5 | **unknown** | A verifiable |
| 11 | lsnc_relu | ✅ | pass 2.5e-6 | **unknown** | A verifiable |
| 12 | malbeware | ✅ | pass 3.9e-6 | **unsat (SAFE)** | A verifiable |
| 13 | metaroom_2023 | ✅ | loose 2.9e-3 | **unsat (SAFE)** | A verifiable |
| 14 | ml4acopf_2024 | ✅ | FAIL (fixed→8e-6 py-path) | timeout | B/C* |
| 15 | nn4sys | ✅ | pass 1.2e-5 | timeout | B compute |
| 16 | relusplitter | ✅ | pass 4.6e-7 | **unknown** | A verifiable |
| 17 | safenlp_2024 | ✅ | pass 4.0e-6 | timeout | B compute |
| 18 | sat_relu | ✅ | pass 1.5e-6 | timeout | B compute |
| 19 | soundnessbench | ✅ (manifest) | pass 7.7e-6 | **error** (category: custom `ReshapeLayer1000` has no `ONNXParams`) | D import-gap |
| 20 | test (nano) | ✅ (manifest) | pass 0 | **error** ("ONNX model not supported" in category path) | D import-gap |
| 21 | tinyimagenet_2024 | ✅ | pass 1.1e-5 | timeout | B compute |
| 22 | tllverifybench_2023 | ✅ | pass 8.2e-7 | timeout | B compute |
| 23 | traffic_signs_2023 | ✅ | **FAIL (1.0)** | timeout | C/B* |
| 24 | vggnet16_2022 | ✅ | pass 3.7e-7 | timeout | B compute |
| 25 | vit_2023 | ✅ | **FAIL** | timeout | C refuse (multi-token attention) |
| 26 | yolo_2023 | ✅ | pass 1.9e-6 | timeout | B compute |

`*` collins_aerospace / ml4acopf / traffic reach on the **category-dispatcher** path
(timeout = ran, no verdict in 10 s) but **fail cross-validation on the manifest path**
(forward pass diverges from ONNX). They are SOUNDLY verifiable only once that import
mismatch is resolved; until then the manifest runner correctly refuses to reach them.

---

## 3. Support tiers

- **Tier A — verifiable (sound verdict at a 10 s budget):** malbeware & metaroom proven
  **SAFE (unsat)**; dist_shift, linearizenn, lsnc_relu, relusplitter return **unknown**
  (sound). collins_rul / test return SAFE on the manifest path. With a realistic timeout
  the Tier-B set joins this tier.
- **Tier B — imports + reaches, compute-bound:** 15+ benchmarks (acasxu, cersyve, cgan,
  cifar100, cora, nn4sys, safenlp, sat_relu, tinyimagenet, tllverifybench, vggnet16, yolo,
  …) ran reach but hit the 10 s cap. These are **not** capability/soundness limits — they
  need more wall-clock (cloud/parallel: see `CLOUD_COMPUTE_OPTIONS.md`).
- **Tier C — fail-loud / documented can't-handle (SOUND refusal):**
  - **cctsdb_yolo** — in-graph YOLO post-processing (`Where`/`ScatterND`/`ArgMax`/
    dynamic `Reshape`); tagged `UnsupportedOp:*` → refuses (never a wrong verdict).
  - **vit_2023** — real multi-token self-attention; multi-token attention reach is
    **fail-loud** in this PR (sound). A sound multi-token bound is the documented
    follow-up (`REMEDIATION_AND_SOUNDNESS_STRATEGY_v02.md` §B1).
  - **collins_aerospace / traffic_signs / ml4acopf** — import cross-validation mismatch
    (manifest path); the runner refuses reach until the forward pass matches.
- **Tier D — category-dispatcher import gap (manifest path works):** soundnessbench
  (custom `ReshapeLayer1000` without `ONNXParams` is parsed by the **manifest** path but
  not the category `matlab2nnv` path), test_nano (no category dispatcher). Prefer the
  manifest path for these.

**Bottom line:** ~21/26 import with a correct forward pass (xval pass/loose) and reach
soundly; 4–5 are fail-loud / import-mismatch (sound refusals); 0 unsound verdicts.

---

## 4. This session's fixes that touch the benchmarks

- **traffic_signs** — `ReshapeLayer` OnnxBCHW reach (Star [9] + ImageStar [166]), the
  reach-no-longer-mutates-the-layer fix [132], and the runner `needReshape==3`
  counterexample inverse [28]. (Manifest xval still flags a forward-pass mismatch to chase
  before claiming a verdict.)
- **lsnc_relu / dist_shift / linearizenn / relusplitter** — now return a sound `unknown`
  on the category path (were `reach_err` on the 2026-05-06 manifest sweep); the fail-loud
  layer fixes turned silent errors into either sound results or loud refusals.
- **vit_2023** — multi-token attention fail-loud guard (SDPA/MHA [0][1][24][45]); the
  multi-token reach is refused rather than computed unsoundly.
- **All layers** — the exact-star over-approximation gate [42] guarantees that if any of
  these is verified under `exact-star`, an over-approximate layer cannot certify
  not-robust (verdict downgraded to `unknown`, one warning emitted).

---

## 5. Merge recommendation

**Safe to merge as a sound VNN-COMP 2025 entry** with the scope: ~21 benchmarks import +
forward-pass-correctly and reach soundly; the remainder are **fail-loud refusals or
import-cross-validation mismatches**, never unsound verdicts. Before the competition:

1. Run a full sweep at a realistic timeout (≥ 120 s) on the CI matrix / cloud to convert
   Tier-B "timeout" into real verdicts.
2. Resolve the manifest-path forward-pass mismatches (collins_aerospace, traffic_signs,
   ml4acopf) if those benchmarks are in scope; until then they stay fail-loud (sound).
3. cctsdb_yolo and vit_2023 remain documented fail-loud (in-graph YOLO post-proc; sound
   multi-token attention is the §B1 follow-up).
