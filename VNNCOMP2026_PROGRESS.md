# NNV — VNN-COMP 2026 participation: execution progress & direction

*Living tracker of how the [scoring strategy](VNNCOMP2026_STRATEGY.md) is being executed. Status as of
2026-06-11. Deadline: **tool submission June 30, 2026 AoE**; competition July 24–25, Lisbon. NNV is voted
in. The 2025 baseline we're improving on: **6th/7, 697.3 pts** (see
[VNNCOMP2025_RESULTS_ANALYSIS.md](VNNCOMP2025_RESULTS_ANALYSIS.md)).*

## Direction in one paragraph

The 2025 numbers say NNV's deficit is **executional, not theoretical** — CORA and nnenum use the same
reachability family and outscore us. The three levers, in priority order from the data: **(1) find more SAT**
(NNV got 354 vs ~1000 — weak random sampling), **(2) stop the −150 penalty leak** (19 invalid/missing
witnesses), **(3) go faster** (worst startup; timeouts on big nets). The soundness gate is a *competitive
advantage* under the −150 rule. So the plan is **fast-and-sound-or-`unknown`**: cheap gradient falsification
first, then a cheap-to-precise abstraction ladder with refinement — never the exponential exact-star on large
nets.

## Pillar status

| Pillar | Item | Status |
|---|---|---|
| **1 — Falsification (SAT)** | Gradient PGD/FGSM falsifier (feature inputs) | ✅ merged (#305) |
| | needReshape-aware PGD (image/permuted inputs) | 🟢 #306 (auto-merging) |
| | Format-from-layer fix + PGD time-budget (no budget starvation) | 🟢 #306 |
| | Per-category falsification config (nRand, restarts/steps, ε) | ⬜ next |
| | falsify ∥ reach race (parfeval), CW fallback for hard cases | ⬜ |
| **2 — Validity (no −150)** | `validate_witness` re-evaluates every SAT before emitting | ✅ merged (#305) |
| | needReshape-correct witness replay | 🟢 #306 |
| | >3-D (Image3D) permute robustness in flat↔net mapping (Copilot) | 🟢 #306 (1e9717b5e) |
| | **onnxruntime witness guard** (`validate_witness_onnx` + `onnx_replay_check.py`) | ✅ built + tested, awaiting runner wiring |
| **3 — Imprecise + refine + speed** | fast-method ladder (zono→abs-dom, no exact-star for big nets) | ✅ (earlier) |
| | adaptive `relaxFactor` per benchmark | ⬜ next |
| | CEGAR ReLU/input splitting refinement (`NN_reach_refinement`) | ⬜ |
| | bandit method selection + bounds caching | ⬜ |
| | GPU bound propagation / PGD-on-GPU | ⬜ (stretch) |
| **Submission** | `install/prepare/run_instance.sh` (R2026a) | ✅ merged (#305, VNN_COMP2026/) |
| | TUM repeatability dry-run | ⬜ (needs the eval site) |
| **VNN-LIB 2.0 (extended track)** | Python `vnnlib`-pip → `.mat` bridge | ⬜ (regular track needs none) |

## Measured results

- **Baseline (2025 official):** 1082 solved (354 SAT / 728 UNSAT), 6th/7, 697.3 pts.
- **2025-set sweep, pre-PGD (fast-method only):** 120 s/instance, one representative instance per folder →
  **1 sat**, 2 unsat, 14 unknown, rest error/timeout. See [VNNCOMP2025_STATUS.md](VNNCOMP2025_STATUS.md).
- **2025-set sweep, WITH the PGD falsifier (2026-06-11, `results_20260611_103123`):** same 120 s/instance,
  26 folders → **4 sat** (collins_rul, cora, safenlp, sat_relu), 2 unsat (malbeware, metaroom),
  12 unknown, 5 error, 3 timeout. **SAT 1→4 (4×)** from gradient falsification — the headline Pillar-1
  signal; unknowns 14→12. (Representative-instance sweep, not the full per-instance set — directional.)
- **Errors worth a follow-up fix** (from the same sweep): `cgan_2023` and `soundnessbench` both fail with
  `Unrecognized ... 'ONNXParams' for class '...ReshapeLayer1000'` — an importer gap in the ONNX Reshape
  path; `collins_aerospace` "Unsupported Class of Layer"; `cctsdb_yolo` known-WIP.

## Immediate next actions (deadline-ordered)

1. **Land #306** (needReshape PGD + time-budget) → gradient falsification covers image benchmarks.
2. **Per-category falsification + `relaxFactor` tuning** — raise `nRand`/PGD effort for the SAT-likely
   classes; set `relaxFactor` high for the large CNNs (speed). Small, safe, testable.
3. **Run the full cloud sweep** (`.github/workflows/vnncomp-sweep.yml`, official 2025 repo) at a realistic
   timeout to get the true regular-track scorecard with the new falsifier.
4. **Submission dry-run** on the TUM site to surface install/witness-format issues early.
5. **`relaxFactor`/bandit + (stretch) CEGAR refinement** for the UNSAT closer; **vnnlib2 bridge** for the
   4 extended-track 2.0 benchmarks if time allows.

## Key hardening (Pillar 2): onnxruntime witness guard — BUILT + TESTED ✅

`validate_witness` replays the candidate through **NNV's** forward (the same reshape/needReshape convention
the falsifier used). That catches numerical near-misses, but a *systematic* witness-encoding error (the
suspected source of several of the 19 −150 penalties — the runner's own "wrong counterexample writing?" note)
would pass NNV's self-consistent replay yet fail the competition's check. The competition validates the SAT
witness by replaying it through **onnxruntime on the original ONNX model**.

**Done (this session):** the definitive guard now exists and is verified end-to-end:
- `code/nnv/tools/onnx2nnv_python/onnx_replay_check.py` — replays a FLAT (ONNX-order) witness through the
  ORIGINAL .onnx via onnxruntime and reports `VIOLATED`/`OK`/`ERROR` for the unsafe region `G·y ≤ g`.
- `code/nnv/examples/Submission/VNN_COMP2025/validate_witness_onnx.m` — MATLAB wrapper that calls it.
  Verified on acasxu: violated→1, impossible-region→0, multi-HalfSpace any-violates→1.
- Two bugs found + fixed while bringing it up: (a) the Windows case-insensitive FS collapsed the temp
  `G.csv`/`g.csv` into one file (g clobbered G) — renamed to `Gmat.csv`/`gvec.csv`; (b) the Python loader
  now reshapes G against the model's *known* output dim instead of trusting the CSV's parsed rank.

**Remaining:** wire `validate_witness_onnx` into `run_vnncomp_instance.m`'s SAT path as the final gate before
emitting `sat` (emit `unknown` if the ONNX replay doesn't violate). Deferred until #306 merges — #306 also
edits `run_vnncomp_instance.m`, so wiring now would conflict. These three files sit on a separate branch off
master (not in #306).

## What NOT to spend time on

Breadth/coverage — NNV already participated in 15/16 benchmarks in 2025; coverage is fine. The score comes
from **SAT-finding, witness validity, and speed**, in that order.
