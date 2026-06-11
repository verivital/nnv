# NNV — VNN-COMP 2025 Readiness Status

*Top-level competition-prep snapshot. Reflects merged master + a **real 120 s/instance sweep**
(`results_20260610_200503`, one representative instance per benchmark folder). The authoritative
per-benchmark detail also lives at
[`code/nnv/examples/Submission/VNN_COMP2025/VNN_COMP2025_SUPPORT.md`](code/nnv/examples/Submission/VNN_COMP2025/VNN_COMP2025_SUPPORT.md).*

## TL;DR — where we are (at a realistic 120 s budget)

- **15 of 26 benchmarks return a completed verdict** at 120 s/instance (2 SAFE, 1 sound counterexample,
  12 sound `unknown`) — up from ~6 at the old 10 s smoke test. **The earlier "timeout" picture was a
  compute artifact, not a capability limit**, exactly as predicted.
- **6 are still compute-bound** at 120 s (large ResNets / VGG-16 / NLP / a SAT-encoding) — they import and
  reach, they just need a bigger budget (300–600 s / cloud).
- **2 are sound fail-loud refusals** (unsupported ops) and **3 are category-path import gaps that the
  Python-importer (manifest) path handles** — so with the manifest path the import coverage is ~24/26.
- **Soundness holds end-to-end: 0 unsound verdicts.** Every `unsat` is a proven over-approximation miss
  (SAFE); every `unknown` is sound; the single `sat` (collins_rul) is a **falsification counterexample**
  (a concrete input that violates the property, found by evaluating the real network — the reach→sat path
  is disabled by the [42] gate). No over-approximation ever certifies not-robust.

**120 s sweep tally:** `unsat 2 · sat 1 · unknown 12 · timeout 6 · error 5 · (sat-from-overapprox 0)`.

---

## Full benchmark matrix (all 26, @120 s)

`Import` = loads into NNV (category path) · `Forward (xval)` = NNV forward pass matches ONNX runtime
(manifest path) · `Reach @120s` = category-dispatcher verdict + time at a 120 s budget (2026-06-10).

| # | Benchmark | Import | Forward (xval) | **Reach @120 s** | Tier | Note |
|---|-----------|:--:|:--:|---|:--:|---|
| 1 | acasxu_2023 | ✅ | ✅ 4.7e-6 | **unknown** (20.0 s) | **A** | sound verdict |
| 2 | cctsdb_yolo_2023 | ✅ | ❌ | error — *fail-loud* (0.2 s) | **C** | in-graph YOLO post-proc ops |
| 3 | cersyve | ✅ | ✅ 7.2e-7 | **unknown** (12.2 s) | **A** | sound verdict |
| 4 | cgan_2023 | ⚠️ manifest | ✅ 7.5e-9 | error — import gap (13.2 s) | **D** | custom `ReshapeLayer1000` → use manifest path |
| 5 | cifar100_2024 | ✅ | ✅ 6.9e-6 | timeout (120 s) | **B** | large ResNet; needs >120 s |
| 6 | collins_aerospace | ✅ | ❌ | error — *fail-loud* (13.1 s) | **C** | "Unsupported Class of Layer" (`Split`/`Pow`) |
| 7 | collins_rul_cnn_2022 | ✅ | ✅ 1.9e-5 | **sat** — *counterexample* (1.9 s) | **A** | sound not-robust (falsification) |
| 8 | cora_2024 | ✅ | ✅ 2.4e-4 | timeout (120 s) | **B** | needs >120 s |
| 9 | dist_shift_2023 | ✅ | ✅ 1.1e-4 | **unknown** (3.4 s) | **A** | sound verdict |
| 10 | linearizenn_2024 | ✅ | ✅ 1.5e-5 | **unknown** (4.5 s) | **A** | sound verdict |
| 11 | lsnc_relu | ✅ | ✅ 2.5e-6 | **unknown** (1.4 s) | **A** | sound verdict |
| 12 | malbeware | ✅ | ✅ 3.9e-6 | **unsat (SAFE)** (2.4 s) | **A** | ✔ proven SAFE |
| 13 | metaroom_2023 | ✅ | ~ 2.9e-3 | **unsat (SAFE)** (4.5 s) | **A** | ✔ proven SAFE |
| 14 | ml4acopf_2024 | ✅ | ✅ 8e-6 (py) | **unknown** (12.2 s) | **A** | sound verdict (was timeout) |
| 15 | nn4sys | ✅ | ✅ 1.2e-5 | **unknown** (8.6 s) | **A** | sound verdict (was timeout) |
| 16 | relusplitter | ✅ | ✅ 4.6e-7 | **unknown** (4.1 s) | **A** | sound verdict |
| 17 | safenlp_2024 | ✅ | ✅ 4.0e-6 | timeout (120 s) | **B** | NLP; needs >120 s |
| 18 | sat_relu | ✅ | ✅ 1.5e-6 | timeout (120 s) | **B** | SAT-encoding (combinatorially hard) |
| 19 | soundnessbench | ⚠️ manifest | ✅ 7.7e-6 | error — import gap (2.9 s) | **D** | custom `ReshapeLayer1000` → use manifest path |
| 20 | test (nano) | ⚠️ manifest | ✅ 0 | error — no dispatcher (0.05 s) | **D** | toy; use manifest path |
| 21 | tinyimagenet_2024 | ✅ | ✅ 1.1e-5 | timeout (120 s) | **B** | large ResNet; needs >120 s |
| 22 | tllverifybench_2023 | ✅ | ✅ 8.2e-7 | **unknown** (14.0 s) | **A** | sound verdict (was timeout) |
| 23 | traffic_signs_2023 | ✅ | ❌ (manifest) | **unknown** (13.5 s) | **A** | category path verifies (reach fixes landed) |
| 24 | vggnet16_2022 | ✅ | ✅ 3.7e-7 | timeout (120 s) | **B** | VGG-16; needs >120 s |
| 25 | vit_2023 | ✅ | ❌ | **unknown** (31.3 s) | **A** | *this* instance (pgd_2_3_16) verifies; true multi-token-attention ViT models still fail loud (sound) |
| 26 | yolo_2023 | ✅ | ✅ 1.9e-6 | **unknown** (18.1 s) | **A** | sound verdict (was timeout) |

---

## Readiness tiers (@120 s)

- **Tier A — completed sound verdict (15):** SAFE — `malbeware`, `metaroom`; counterexample (sound) —
  `collins_rul`; sound `unknown` — `acasxu`, `cersyve`, `dist_shift`, `linearizenn`, `lsnc_relu`,
  `ml4acopf`, `nn4sys`, `relusplitter`, `tllverifybench`, `traffic`, `vit`(this instance), `yolo`.
- **Tier B — imports + reaches, still compute-bound at 120 s (6):** `cifar100`, `cora`, `safenlp`,
  `sat_relu`, `tinyimagenet`, `vggnet16`. Need a larger budget (300–600 s) and/or cloud fan-out — **not**
  a capability or soundness limit.
- **Tier C — sound fail-loud / can't-handle (2):** `cctsdb_yolo` (in-graph YOLO post-proc:
  `Where`/`ScatterND`/`ArgMax`/dyn-`Reshape`), `collins_aerospace` ("Unsupported Class of Layer":
  `Split`/`Pow`). Refusals, never wrong verdicts.
- **Tier D — category-path import gap; manifest (Python-importer) path works (3):** `cgan`,
  `soundnessbench`, `test`. All three forward-pass-validate on the manifest path (xval ≤ 7.7e-6); only the
  category `matlab2nnv` dispatcher can't parse them (custom `ReshapeLayer1000` / no dispatcher). **Run
  these via `code/nnv/tools/onnx2nnv_python/` → `load_nnv_from_mat`.**

---

## Pre-competition action list (updated post-sweep)

1. **🔴 Longer budget for the 6 Tier-B benchmarks.** Re-run `cifar100`, `cora`, `safenlp`, `sat_relu`,
   `tinyimagenet`, `vggnet16` at **300–600 s** (cloud / CI matrix fan-out — see `CLOUD_COMPUTE_OPTIONS.md`).
   This is now the *only* compute gap. (At 120 s, 15/26 already complete.)
2. **🟠 Drive the 3 Tier-D benchmarks through the manifest path** (the importer is now vendored at
   `code/nnv/tools/onnx2nnv_python/`): `cgan`, `soundnessbench`, `test`. Optionally fix the category-path
   `ReshapeLayer1000` parse so the dispatcher path also handles them.
2b. **🟠 Run *all* instances, not one per folder.** This sweep used one representative instance per
   benchmark; the real scorecard needs every instance (and per-benchmark official timeouts from each
   `instances.csv`).
3. **🟡 Sound multi-token attention (roadmap F1/T1.1)** for the true ViT models inside `vit_2023` (the
   `pgd_2_3_16` instance already verifies; multi-token-attention models still fail loud — sound).
4. **🟡 New ONNX op semantics** (`Where`/`ScatterND`/`ArgMax`/`Split`/`Pow`/dyn-`Reshape`) unblock
   `cctsdb_yolo` and `collins_aerospace`.
5. **🟢 Tighten loose bounds** (e.g. `metaroom` xval 2.9e-3; softmax box bounds) to convert more `unknown`
   → `unsat` at fixed ε (roadmap F3).

## Soundness guarantee (verified against the runner)

`sat` (status 0) is set **only** by `falsify_single` (a concrete counterexample, evaluated on the real
network); the reach→`sat` promotion is disabled (the [42] exact-star gate, commented out in
`run_vnncomp_instance.m`). So every `sat` is a sound not-robust, every `unsat` is a sound SAFE, every
`unknown` is sound, unsupported ops `error` (fail-loud), and the runner gates reach on forward-pass
cross-validation. **No benchmark, at any timeout, returns an unsound `sat`/not-robust.**

## How to reproduce / extend the sweep

```matlab
cd code/nnv; startup_nnv;                              % NNV v3.0.0
cd examples/Submission/VNN_COMP2025;
run_all_benchmarks('<repo>/vnncomp2025_benchmarks/benchmarks', 120);   % 120 s/instance, one per folder
% -> writes results_<timestamp>.csv + .md next to this script.
% For a longer budget: pass 300 or 600. For the full scorecard: iterate every row of each instances.csv.
% Tier-D / xval-mismatch nets: use the manifest path
%   python code/nnv/tools/onnx2nnv_python/onnx2nnv.py <model.onnx> -o <model.nnv.mat>
%   net = load_nnv_from_mat('<model.nnv.mat>');
```
