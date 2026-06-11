# NNV — VNN-COMP 2026 Readiness Plan

*Forward plan for NNV's VNN-COMP 2026 entry. Sourced from a deep `gh`/web sweep of the VNN-COMP 2026
governance repo, the 2026 benchmarks repo, and the VNN-LIB 2.0 spec (2026-06-11). Companion to the 2025
status in [`VNNCOMP2025_STATUS.md`](VNNCOMP2025_STATUS.md).*

## TL;DR — and the hard deadline

- **VNN-COMP 2026** is July 24–25, 2026 in Lisbon (SAIV @ FLoC). **🔴 Tool submission closes June 30, 2026
  AoE** — a hard deadline. NNV is among the voted-in tools.
- **The regular track (24 benchmarks) is ALL VNN-LIB 1.0** → reachable with NNV's **existing**
  `load_vnnlib.m` parser. **No VNN-LIB 2.0 parser is strictly required to score the regular track.**
- **VNN-LIB 2.0** is the big new capability: **4 extended-track benchmarks are 2.0-only**
  (`adaptive_cruise_control_non_linear_2026`, `isomorphic_acasxu_2026`, `monotonic_acasxu_2026`,
  `smart_turn_multimodal_2026`). 2.0 adds `declare-network` (multi-net), `declare-hidden` (intermediate-
  tensor constraints), multi-input/output, `X[i,j]` tensor indexing, a type system, `equal-to`/
  `isomorphic-to`, nested `and`/`or` (DNF query trees), a `(vnnlib-version 2.0)` header, and `.vnnlib.gz`
  compression. Spec: vnnlib.org + arXiv:2605.07451; Python `vnnlib` pip v1.0.2 parses 2.0 (incl. `.vnnlib.gz`).
- **Scoring** (issue #2 / rules.md): correct = +10, **incorrect = −150** (so an unsound verdict is
  catastrophic — NNV's sound-by-refusal posture is exactly right), timeout/error/unknown = 0. SAT needs a
  witness (line 2: `((VAR val) …)`) replayable through onnxruntime within 1e-3 rel / 1e-4 abs.

## Benchmark delta vs 2025

- **Regular track (24, VNN-LIB 1.0):** the 2025 set NNV already handles, plus refreshed/new ones —
  `acasxu_2023`, `cersyve`, `cgan2026`, **`challenging_certified_training_2026` (NEW)**, `cifar100_2024`,
  `collins_rul_cnn_2022`, `cora_2024`, `dist_shift_2023`, `linearizenn_2024`, `lsnc_relu`, `malbeware`,
  `metaroom_2023`, `ml4acopf_2024`, `nn4sys`, **`relusplitter_2026`**, `safenlp_2024`, `sat_relu`,
  **`soundnessbench_2026`**, `tinyimagenet_2024`, `tllverifybench_2023`, `traffic_signs_recognition_2023`,
  `vggnet16_2022`, `vit_2023`, `yolo_2023`. Bot re-pushed ONNX/csv June 4–10, so **re-validate all against
  the latest repo**.
- **Extended track (6):** `cctsdb_yolo_2023` + `collins_aerospace_benchmark` (1.0; NNV already fail-louds
  these soundly), and the **4 VNN-LIB-2.0-only** new ones above.
- **6 new for 2026:** `adaptive_cruise_control_non_linear_2026` (2.0), `cersyve` (control: cartpole/double-
  integrator/lane-keep/pendulum/pointmass), `challenging_certified_training_2026`, `isomorphic_acasxu_2026`
  (2.0, network-equivalence), `monotonic_acasxu_2026` (2.0), `smart_turn_multimodal_2026` (2.0, vision+control).
- Known repo issues to watch (vnncomp2026_benchmarks): #5 1.0 files in 2.0 folders, #4 wrong ONNX paths in
  `monotonic_acasxu` instances.csv, #2 `cgan2026` version ambiguity (repo folder may appear as `cgan_2026`).

## VNN-LIB 2.0 — what `load_vnnlib.m` needs (for the extended track)

NNV's `load_vnnlib.m` is an S-expression parser (keeps working for all 1.0). To add 2.0:
1. **Version dispatch:** read `(vnnlib-version 1.0|2.0)` header; branch parsing. Add **`.vnnlib.gz`
   decompression** (gunzip on load) — 2026 ships specs gzipped.
2. **`declare-input/output/hidden` + types + shapes:** parse named, typed (`float32`…), shaped tensors;
   replace the implicit single `X`/`Y` with named multi-input/output maps.
3. **Tensor indexing `X[i,j,k]`** (vs 1.0 `X_i_j_k`): map multi-dim subscripts to NNV's linear/predicate
   indices.
4. **`declare-hidden`:** constrain an intermediate ONNX **node** output by name — needs NNV to expose
   named intermediate reach sets (non-trivial; biggest lift, used by the attention/encoder 2.0 benchmarks).
5. **Nested `and`/`or` (DNF query trees):** build full disjunctive `HalfSpace` sets (1.0 had limited
   disjunction).
6. **`declare-network` / `equal-to` / `isomorphic-to`:** multi-network queries (teacher-student,
   ACAS-equivalence). For the isomorphic/monotonic ACAS benchmarks, parse + represent two networks.
- **Shortcut:** the Python `vnnlib` pip package (v1.0.2) already parses 2.0 — consider a thin Python→`.mat`
  bridge (like the existing `onnx2nnv_python` importer) that emits a 1.0-style bounds/spec `.mat` NNV can
  consume, rather than reimplementing the full 2.0 grammar in MATLAB by June 30.

## Tool submission (rules.md)

- Three bash scripts: **`install_tool.sh`** (once: deps/licenses), **`prepare_instance.sh`** (version,
  benchmark, onnx, vnnlib; ≤10 min; no analysis), **`run_instance.sh`** (…+ results path + timeout s;
  writes a single word: `sat`/`unsat`/`unknown`/`timeout`/`error`, + SAT witness). Submit via the TUM
  repeatability site; start unofficial dry-runs NOW to surface install/ONNX errors early.
- **Platform (pick one):** recommend **`m5.16xlarge` CPU** (64 vCPU, 256 GB) — NNV is CPU/LP/star-based and
  the large ViT/VGG/TinyImageNet nets need RAM the GPU box (61 GB) lacks.

## Prioritized NNV action list (deadline-driven)

1. **🔴 NOW — submission harness + dry run.** Write `install/prepare/run_instance.sh` around the existing
   `run_vnncomp_instance` path; verify the SAT **witness format** + the dual CE tolerance (1e-3/1e-4);
   start unofficial dry-runs on the TUM site. (Biggest schedule risk is install/format, not verification.)
2. **🔴 Regular-track refresh (24, VNN-LIB 1.0).** Re-validate the 13 already-wired benchmarks against the
   June ONNX/csv refresh, and finish loaders for `cgan2026`, `relusplitter_2026`,
   `challenging_certified_training_2026`, `cersyve`, `soundnessbench_2026`. Carry over the 2025 wins
   (fast-method/no-exact-star for the compute-bound nets; manifest path for import-gap nets).
2b. **🟠 Run the full cloud sweep** (`.github/workflows/vnncomp-sweep.yml`, now defaulting to the official
   repo) at a realistic timeout to get the real regular-track scorecard before June 30.
3. **🟡 VNN-LIB 2.0 support (extended track).** Add `.vnnlib.gz` + version dispatch + named typed
   multi-IO + tensor indexing to `load_vnnlib.m` — OR a Python `vnnlib`→`.mat` bridge (faster path). Defer
   `declare-hidden` (intermediate constraints) and multi-network (`isomorphic/monotonic acasxu`) unless
   time allows; these are the deepest lifts.
4. **🟢 −150 penalty discipline.** The catastrophic-incorrect rule makes NNV's "never certify
   not-robust/unsafe from an over-approximation" gate a competitive *advantage* — keep every reach
   sound-or-refuse, and make sure SAT comes only from validated falsification witnesses.

## What's unchanged / reusable from the 2025 work

The fast-method runner (no exact-star for compute-bound cats), the manifest (Python-importer) path for
import-gap nets, the soundness gate ([42]), and the cloud sweep workflow all carry directly into 2026 — the
regular track is largely the 2025 pipeline re-pointed at `vnncomp2026_benchmarks`.
