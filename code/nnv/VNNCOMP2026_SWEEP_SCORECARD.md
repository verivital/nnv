# VNN-COMP 2026 Sweep Scorecard + Coverage Analysis

Source: full cloud sweep of the **2026 benchmark set** (run 27394691599, 2026-06-12),
`run_which=all`, per-instance timeout 120 s, dispatched from the #331 branch (this is
**before** the #334 VNN-LIB 2.0 parser merged, so the 2.0-only benchmarks error here).

## Headline

- **3373 instances** | **1242 solved** (856 sat / 386 unsat) | 1002 unknown | 510 timeout | **619 error**.
- **NNV-new solves 1242 vs the NNV-2025 submission's 1084** — net improvement from the year's work.
- **SOUNDNESS: clean.** Of 1213 definitive verdicts cross-checked against the field,
  **all 1213 agree with the majority; 0 false-sat, 0 false-unsat, 0 disagreement with
  alpha-beta-CROWN.** The sound-or-unknown posture (the −150 guard) is holding.

The score opportunity is almost entirely in the **619 errors** (whole benchmarks that
fail to *run*) and secondarily the 510 timeouts. Errors are 0 points, not −150, so this
is pure upside — recovering them never risks soundness.

## Errors by root cause (priority order)

### P1 + P2 — manifest-path benchmarks starved by TWO sweep-measurement bugs (FIXED): cgan2026, cgan_2023, soundnessbench, soundnessbench_2026, lsnc_relu, traffic_signs (~273 instances)
All six route to `load_manifest_net` (the xval-gated Python-importer path), which needs a
`<model>.nnv.mat` generated alongside the ONNX. Root-caused to two **sweep-harness** bugs
(NOT the verifier — the competition's own folder-name categories + decompressed ONNX dodge both):
1. **Manifest-gen glob mismatch.** The 2026 set ships ONNX **gzipped** (`<model>.onnx.gz`), but
   the sweep's manifest-gen step globbed `*.onnx` → with `nullglob` it matched nothing and
   generated **zero** manifests. lsnc_relu/traffic_signs then hit `NNV manifest not found`;
   cgan/soundnessbench's earlier matlab2nnv error was actually the *category* bug below routing
   them away from the manifest path entirely.
2. **Category map mismatch** (`run_all_benchmarks.m`): `soundnessbench`→`'soundness'` (the
   routing checks `contains(category,"soundnessbench")` → never matched), and the new 2026
   folders (cgan2026, soundnessbench_2026, …) weren't mapped → category `'?'` → no dispatch.
**FIXED** (this PR): manifest-gen now `gunzip`s the `.onnx.gz` first; the category map fixes
`soundnessbench` and defaults unmapped folders to the folder name (mirroring the competition).
**Validated locally**: the Python importer parses soundnessbench's fused `Gemm_To_ReshapeLayer`
cleanly into 17 layers, and `run_vnncomp_instance('soundnessbench', …)` now routes through the
manifest and reaches (import error gone) — a re-run will measure the recovered verdicts.
(NB: the *separate* `matlab2nnv` `ReshapeLayer.parse` `.ONNXParams` gap on fused Gemm+Reshape
layers is real but moot for these benchmarks, since the manifest path bypasses matlab2nnv.)

### P3 — "ONNX model not supported" (importer gaps): cgan2026 (overlaps P1), challenging_certified_training_2026 (30), relusplitter_2026 (120)
`matlab2nnv`/`onnx2nnv` rejects the model outright. Needs per-architecture importer work; assess
each model's unsupported op. relusplitter_2026 (120) is the largest single bucket.

### P4 — Explicit "not supported" / unsupported layer: cctsdb_yolo_2023 (39, "Working on supporting this one"), collins_aerospace_benchmark (6, "Unsupported Class of Layer" — yolov5 LReLU)
Object-detector nets (YOLO variants); known hard. Lower priority.

### 2.0-only benchmarks (1 instance each, error here): adaptive_cruise_control_non_linear_2026, isomorphic_acasxu_2026, monotonic_acasxu_2026, smart_turn_multimodal_2026
These errored because the sweep predates **#334**. With the 2.0 parser now on master, a re-run
would: parse the single-net 2.0 forms, and **gate** the multi-net (monotonic/isomorphic),
nonlinear (acc), and multimodal (smart_turn) ones to `unknown` (sound). Net: no errors, a few
`unknown`. Tiny instance counts, but confirms the 2.0 dispatch end-to-end.

## Regressions vs the NNV-2025 column (same-named benchmarks — investigate)

| benchmark | 2025 solved | new solved | new status | likely cause |
|---|--:|--:|---|---|
| cgan_2023 | 16 | 0 | 19 error | P1 ReshapeLayer |
| cifar100_2024 | 190 | 0 | 157 timeout (job failed) | slowdown / OOM; was all-unsat in 2025 |
| yolo_2023 | 71 | 0 | 72 unknown | reach no longer reaching a verdict |
| ml4acopf_2024 | 17 | 0 | 69 unknown | cp-star path |
| nn4sys | 17 | 2 | 105 unknown / 65 error | mixed |
| dist_shift_2023 | 54 | 5 | 67 unknown | reach tuning |

(2025 numbers are on the 2025 *benchmarks*; 2026 instances may differ, so treat as directional.)
cifar100 + tinyimagenet + malbeware were the 3 **job-level** failures (OOM/disk/timeout at the
runner level, not per-instance) — cifar100/tinyimagenet still emitted partial timeout rows;
malbeware produced no rows at all (failed before writing results — diagnose from the job log).

## Recommended next actions (ranked)

1. **P1 ReshapeLayer fix** — localize the `.ONNXParams`→`.Vars` access; validate by importing a
   cgan + a soundnessbench model locally. ~148 instances, sound, well-localized.
2. **P2 manifest-gen in the sweep** — ensure lsnc_relu/traffic_signs `.nnv.mat` are generated in
   the cloud before the run. 125 instances, infra-only.
3. **Timeout tuning** — acasxu (87 to), cora (68 to), cifar100/tinyimagenet: per-category reach
   method / time-budget tuning (allowed by the 2026 rules; key on `category` only).
4. **Re-run the sweep on master** (post-#334/#335) to confirm the 2.0 dispatch + the option-3
   leak fix, and to re-measure after P1/P2.
5. Investigate the yolo_2023 / ml4acopf / dist_shift unknown-regressions (reach reaching no verdict).

Soundness is the win to protect: every recovery above must stay sound-or-unknown.
