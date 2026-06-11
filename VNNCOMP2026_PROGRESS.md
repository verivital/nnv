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

*(2026-06-11 implementation push: PRs #306–#311 all merged/merging to `verivital/nnv:master`.)*

| Pillar | Item | Status |
|---|---|---|
| **1 — Falsification (SAT)** | Gradient PGD/FGSM falsifier (feature inputs) | ✅ merged (#305) |
| | needReshape-aware PGD (image/permuted inputs) + format-from-layer + time-budget | ✅ merged (#306) |
| | >3-D (Image3D) permute robustness in flat↔net mapping (Copilot) | ✅ merged (#306) |
| | **Per-category falsification config** (PGD restarts/steps/max_time + nRand table) | ✅ merged (#309) |
| | nRand tuned DOWN from measured ~15 ms/sample (acasxu falsify 21.7→8.6 s) | ✅ merged (#309) |
| | **PGD on the NN-manifest path** (numerical grad: finite-diff + SPSA) | 🟢 #311 (lsnc_relu unknown→**SAT**) |
| | falsify ∥ reach race (parfeval), CW fallback for hard cases | ⬜ (future) |
| **2 — Validity (no −150)** | `validate_witness` re-evaluates every SAT before emitting | ✅ merged (#305) |
| | needReshape-correct witness replay + >3-D robustness | ✅ merged (#306) |
| | **onnxruntime witness guard** (`validate_witness_onnx` + `onnx_replay_check.py`) | ✅ merged (#307); runner-wiring DEFERRED (see below) |
| **3 — Imprecise + refine + speed** | fast-method ladder (zono→abs-dom, no exact-star for big nets) | ✅ (earlier) |
| | **Fix 4 `reachOptionsList{1}`-overwrite bugs** (dist_shift, linearize, cora, malbeware) | ✅ merged (#309) |
| | **acasxu reach recovered** (ImageStar by input-layer type, not shape) | ✅ merged (#310) |
| | per-category `relaxFactor` (set in the #309 reach ladders) | ✅ merged (#309) |
| | CEGAR ReLU/input splitting refinement; bandit method selection; GPU | ⬜ (future) |
| **CI reliability** | Disk-reclaim before Setup MATLAB (fix "No space left on device" flake) | ✅ merged (#306 test-matrix, #308 ci+regression) |
| **Submission** | `install/prepare/run_instance.sh` (R2026a) | ✅ merged (#305) |
| | TUM repeatability dry-run | ⬜ (needs the eval site) |
| **VNN-LIB 2.0 (extended track)** | Python `vnnlib`-pip → `.mat` bridge | ⬜ (regular track needs none) |

### Deferred / next (judgment calls made this session)
- **onnxruntime guard runner-wiring**: built + merged as a standalone tool (#307), but NOT wired into the SAT path. Wiring it naively would emit `unknown` for a valid `sat` whenever the eval machine lacks onnxruntime (suppressing real SAT points). Needs a fallback that trusts `validate_witness` when onnxruntime is unavailable. Tracked.
- **importer ReshapeLayer `.Vars`** (cgan, soundnessbench ERROR): new-style `importONNXNetwork` reshape layers store the shape in `.Vars` (+ permutes), not `.ONNXParams`; `ReshapeLayer.parse` only handles the old style. Coverage work — strategy deprioritizes it.
- **acasxu exact-star scalability**: now that acasxu reach RUNS (#310), exact-star can blow up on its 6×50 ReLU net; competition timeout bounds it (easy instances solved, hard ones time out as before). A bounded/relax fallback could help.
- **Full validation sweep** with all of #306–311 to quantify the cumulative scorecard lift.

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
- **Post-implementation validation sweep (2026-06-11, all of #306–#311 on `master`):** focused
  correct-pairing sweep (via each benchmark's `instances.csv`) of the changed benchmarks →
  **4 sat** (lsnc_relu, cora, linearizenn, collins_rul), 4 unknown, **0 errors / 0 regressions**.
  lsnc_relu is the #311 win (unknown→SAT via numerical-gradient PGD); cora + linearizenn confirm the
  #309 reach-overwrite fixes did NOT regress (both still find SAT). (Naive onnx/vnnlib pairing — NOT using
  instances.csv — produces spurious dimension-mismatch errors; always pair via instances.csv.)

## Immediate next actions (deadline-ordered)

The 2026-06-11 push (**#306–#311, all merged**) landed the SAT/validity/speed core: gradient PGD
(feature + image + NN-manifest), per-category falsification config + measured `nRand` tuning, the four
`reachOptionsList{1}`-overwrite reach fixes, acasxu reach recovery, the onnxruntime witness guard, and CI
disk-reclaim across all workflows. Remaining, deadline-ordered:

1. **Run the full cloud sweep** (`.github/workflows/vnncomp-sweep.yml`, official 2025 repo, correct
   instances.csv pairing) at a realistic timeout to get the true regular-track scorecard with everything landed.
2. **Submission dry-run** on the TUM site to surface install/witness-format issues early.
3. **onnxruntime-guard runner-wiring with a fallback** (see below) — the definitive −150 closer.
4. **(stretch)** CEGAR refinement / bandit method selection for the UNSAT closer; **importer ReshapeLayer
   `.Vars`** for cgan/soundnessbench; **vnnlib2 bridge** for the 4 extended-track 2.0 benchmarks if time allows.

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

**Remaining:** the guard is **MERGED as a standalone tool (#307)** but NOT wired into the runner's SAT path.
Wiring it naively would emit `unknown` for a valid `sat` whenever the eval machine lacks onnxruntime
(suppressing real SAT points — a net loss). The wiring needs a 3-state result (violated / not-violated /
onnxruntime-unavailable) so it DOWNGRADES to `unknown` only on a confirmed non-violation and otherwise trusts
`validate_witness`. Tracked as the next Pillar-2 item.

## What NOT to spend time on

Breadth/coverage — NNV already participated in 15/16 benchmarks in 2025; coverage is fine. The score comes
from **SAT-finding, witness validity, and speed**, in that order.
