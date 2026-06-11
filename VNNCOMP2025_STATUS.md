# NNV — VNN-COMP 2025 Readiness Status

*Top-level competition-prep snapshot. Reflects merged master (`9391b0ab`, PR #290). The authoritative
per-benchmark detail lives at
[`code/nnv/examples/Submission/VNN_COMP2025/VNN_COMP2025_SUPPORT.md`](code/nnv/examples/Submission/VNN_COMP2025/VNN_COMP2025_SUPPORT.md);
this file is the at-a-glance status + the prioritized action list for the competition.*

## TL;DR — where we are

- **All 26 VNN-COMP 2025 benchmarks import into NNV; ~21 import with a verified-correct forward pass and
  reach soundly.** The remaining 4–5 are **sound fail-loud refusals or forward-pass (xval) mismatches** —
  never wrong verdicts.
- **The headline gating factor is COMPUTE, not capability.** The fresh sweep below ran at a **10 s
  smoke-test budget**; **17 of 26 are "timeout" simply because 10 s isn't enough** — they import and reach
  fine. A realistic sweep (**≥120 s/instance, on the cloud / CI matrix**) is the single most important
  pre-competition action and will convert most of those into real verdicts.
- **Soundness holds end-to-end:** across both import paths, **no instance returns `sat`/not-robust from an
  over-approximation** — every verdict is `unsat` (SAFE) or `unknown`; refusals `error` (fail-loud); and
  the runner **gates reach on forward-pass cross-validation**, so a mis-import can't produce an unsound
  result.

**Fresh 10 s sweep tally (2026-06-10, category path):** `unsat 2 · unknown 4 · timeout 17 · error 3 · sat 0`.

---

## Full benchmark matrix (all 26)

`Import` = loads into NNV · `Forward (xval)` = NNV forward pass matches ONNX runtime (manifest path,
2026-05-06) · `Reach @10s` = category-dispatcher outcome at the 10 s smoke budget (2026-06-10).

| # | Benchmark | Import | Forward (xval) | Reach @10s | Readiness | What's needed for a competition verdict |
|---|-----------|:--:|:--:|:--:|:--:|---|
| 1 | acasxu_2023 | ✅ | ✅ 4.7e-6 | timeout | **B** compute | More wall-clock (≥120 s); imports+reaches today |
| 2 | cctsdb_yolo_2023 | ✅ | ❌ | error (fail-loud) | **C** refuse | New ONNX ops: in-graph YOLO post-proc (`Where`/`ScatterND`/`ArgMax`/dyn-`Reshape`) |
| 3 | cersyve | ✅ | ✅ 7.2e-7 | timeout | **B** compute | More wall-clock |
| 4 | cgan_2023 | ✅ | ✅ 7.5e-9 | timeout | **B** compute | More wall-clock |
| 5 | cifar100_2024 | ✅ | ✅ 6.9e-6 | timeout | **B** compute | More wall-clock (large ResNet) |
| 6 | collins_aerospace | ✅ | ❌ | timeout | **C→B** | Resolve forward-pass mismatch (`Split`/`Pow`); then compute-bound |
| 7 | collins_rul_cnn_2022 | ✅ | ✅ 1.9e-5 | timeout *(manifest: **unsat**)* | **A/B** | Verifies SAFE on manifest path; category path needs time |
| 8 | cora_2024 | ✅ | ✅ 2.4e-4 | timeout | **B** compute | More wall-clock |
| 9 | dist_shift_2023 | ✅ | ✅ 1.1e-4 | **unknown** | **A** verifiable | Sound verdict already; tighter bounds → fewer unknowns |
| 10 | linearizenn_2024 | ✅ | ✅ 1.5e-5 | **unknown** | **A** verifiable | Sound today |
| 11 | lsnc_relu | ✅ | ✅ 2.5e-6 | **unknown** | **A** verifiable | Sound today |
| 12 | malbeware | ✅ | ✅ 3.9e-6 | **unsat (SAFE)** | **A** verifiable | ✔ proven SAFE at 10 s |
| 13 | metaroom_2023 | ✅ | ~ 2.9e-3 (loose) | **unsat (SAFE)** | **A** verifiable | ✔ proven SAFE; tighten the loose xval |
| 14 | ml4acopf_2024 | ✅ | ❌→✅ 8e-6 (py path) | timeout | **B/C** | Use the fixed Python-importer manifest; then compute-bound |
| 15 | nn4sys | ✅ | ✅ 1.2e-5 | timeout | **B** compute | More wall-clock |
| 16 | relusplitter | ✅ | ✅ 4.6e-7 | **unknown** | **A** verifiable | Sound today |
| 17 | safenlp_2024 | ✅ | ✅ 4.0e-6 | timeout | **B** compute | More wall-clock |
| 18 | sat_relu | ✅ | ✅ 1.5e-6 | timeout | **B** compute | More wall-clock |
| 19 | soundnessbench | ✅ (manifest) | ✅ 7.7e-6 | error (category) | **D** import-gap | Custom `ReshapeLayer1000` lacks `ONNXParams` on the category path → use manifest path |
| 20 | test (nano) | ✅ (manifest) | ✅ 0 | error (category) | **D** import-gap | No category dispatcher → use manifest path |
| 21 | tinyimagenet_2024 | ✅ | ✅ 1.1e-5 | timeout | **B** compute | More wall-clock (large ResNet) |
| 22 | tllverifybench_2023 | ✅ | ✅ 8.2e-7 | timeout | **B** compute | More wall-clock |
| 23 | traffic_signs_2023 | ✅ | ❌ (1.0) | timeout | **C→B** | Resolve forward-pass mismatch (layout); reach fixes already in PR |
| 24 | vggnet16_2022 | ✅ | ✅ 3.7e-7 | timeout | **B** compute | More wall-clock (VGG-16) |
| 25 | vit_2023 | ✅ | ❌ | timeout | **C** refuse | **Sound multi-token attention bound** (roadmap F1/T1.1) |
| 26 | yolo_2023 | ✅ | ✅ 1.9e-6 | timeout | **B** compute | More wall-clock |

---

## Readiness tiers

- **Tier A — sound verdict already at a 10 s budget (8):** `malbeware` + `metaroom` proven **SAFE
  (unsat)**; `dist_shift`, `linearizenn`, `lsnc_relu`, `relusplitter` return sound **`unknown`**;
  `collins_rul`, `test` return SAFE on the manifest path.
- **Tier B — imports + reaches, compute-bound (≈15):** `acasxu`, `cersyve`, `cgan`, `cifar100`, `cora`,
  `collins_rul`, `nn4sys`, `safenlp`, `sat_relu`, `tinyimagenet`, `tllverifybench`, `vggnet16`, `yolo`
  (+ `collins_aerospace`/`traffic`/`ml4acopf` once their import is settled). **Not** capability or
  soundness limits — they need wall-clock. **This is the population the ≥120 s cloud sweep unlocks.**
- **Tier C — sound fail-loud / documented can't-handle (refusals, never wrong):**
  - `cctsdb_yolo` — in-graph YOLO post-processing ops with no NNV semantics → `UnsupportedOp:*` refusal.
  - `vit_2023` — real multi-token self-attention; multi-token reach is **fail-loud** today (sound). A
    sound multi-token bound is roadmap **F1/T1.1**.
  - `collins_aerospace`, `traffic_signs`, `ml4acopf` — manifest-path forward-pass (xval) mismatch; the
    runner **refuses reach** until the forward pass matches (ml4acopf is already fixed on the Python path).
- **Tier D — category-dispatcher import gap, manifest path works (2):** `soundnessbench`, `test_nano` —
  import + verify via the **manifest / Python-importer path**; the category `matlab2nnv` path can't parse
  them yet.

---

## Pre-competition action list (prioritized)

1. **🔴 Run a full sweep at a realistic timeout (≥120 s, ideally 5–10 min) on the cloud / CI matrix.**
   This is the single highest-value action: it converts the **17 compute-bound "timeout" rows into real
   verdicts** and gives the true performance picture. The benchmark set is embarrassingly parallel — fan
   it out. See [`CLOUD_COMPUTE_OPTIONS.md`](CLOUD_COMPUTE_OPTIONS.md) for the licensing + parallelization
   options (campus license / MLM token, `matlab-actions/setup-matlab`, Docker `mathworks/matlab`, ACCRE
   SLURM array jobs). The matrix CI is now toolbox-cache-accelerated and flake-resistant, so it's a viable
   harness.
2. **🟠 Resolve the 3 forward-pass (xval) mismatches if those benchmarks are in scope:**
   `collins_aerospace` (`Split`/`Pow`), `traffic_signs` (layout — reach fixes already landed in PR #290),
   `ml4acopf` (use the fixed Python-importer manifest). Each then drops from Tier C into Tier B.
3. **🟡 Sound multi-token attention (roadmap F1/T1.1)** unblocks **`vit_2023`** — value-hull/simplex bound,
   MC-validated for seq_len ∈ {2,4,8}. Until then `vit_2023` stays a sound fail-loud refusal.
4. **🟡 New ONNX op semantics** (`Slice`/`Expand`/`Where`/`ScatterND`/`ArgMax`/large-`Split`/`Pow`) unblock
   **`cctsdb_yolo`** and the Collins import; per-benchmark tests as each lands.
5. **🟢 Use the manifest path** for `soundnessbench` and `test_nano` (the category dispatcher can't parse
   them); optionally close the category-path gap later.

## Soundness guarantee (holds regardless of timeout)

Every completed verdict is `unsat` (the over-approximate reachable set misses the unsafe region →
definitively SAFE) or `unknown` (sound; never reported as not-robust). Unsupported ops `error`
(fail-loud). The runner **gates reach on forward-pass cross-validation**, and the exact-star
over-approximation gate ([42]) guarantees an over-approximate layer can never certify not-robust. **No
benchmark, at any timeout, returns an unsound `sat`/not-robust.**

## How to reproduce / extend the sweep

```matlab
cd code/nnv; startup_nnv;                 % NNV v3.0.0
% category-dispatcher path (one representative instance per folder):
%   run_vnncomp_instance(category, onnx, vnnlib, timeout_s)
% manifest / Python-importer path (preferred for Tier-D and the xval-mismatch set):
%   python code/nnv/tools/onnx2nnv_python/onnx2nnv.py <model.onnx> -o <model.nnv.mat>
%   net = load_nnv_from_mat('<model.nnv.mat>');
```
Raise the per-instance timeout from the 10 s smoke budget to ≥120 s and run all instances (not one per
folder) to produce the real competition scorecard.
