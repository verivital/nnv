# NNV — VNN-COMP 2026 Submission Readiness Review

*Generated 2026-06-11 by a read-only readiness review (no files modified, no MATLAB run). The authoritative rules were pulled live from the official `VNN-COMP/vnncomp2026` repo (`rules.md`). Sources cited inline by path/URL. Companion to `VNNCOMP2026_READINESS.md`, `VNNCOMP2026_STRATEGY.md`, `VNNCOMP2026_PERCATEGORY_TUNING.md`, `VNNCOMP2026_PROGRESS.md`.*

---

## 0. Framing: substantial prior work already exists

A large amount of the 2026 engineering is already done (gradient PGD falsifier, per-category tuning, witness validation, onnxruntime guard, the four reach-overwrite bug fixes — PRs #305–#319). This review focuses on what the **rules require** and what is **still stale or unfinished** in the `VNN_COMP2026/` scaffolding.

---

## (a) The 2026 rules that matter for us

Source: `https://github.com/VNN-COMP/vnncomp2026/blob/main/rules.md` + `README.md`.

### Scoring
- Correct unsat / correct sat: **+10** each.
- **Incorrect result: −150.** Timeout / Error / Unknown: **0**. No time bonus.
- Per-benchmark score = `100 × (your instance-score sum) / (max sum of any tool on that benchmark)`; overall = sum of per-benchmark percentages (benchmarks weighted equally). Ties → total runtime on solved instances.
- The **−150:+10 = 15:1 penalty ratio** makes NNV's sound-or-`unknown` posture a structural advantage — the single most important rule for our strategy.

### Tuning rule (CONFIRMED)
> "The benchmark name is provided, as **per benchmark tuning/settings are permitted (per instance settings are not, so do NOT use the onnx filename or vnnlib filename to customize the verification tool settings)**."

**Per-benchmark/per-category tuning is ALLOWED; per-instance is FORBIDDEN.** NNV's tuning table keys only on the `category` arg (never the onnx/vnnlib filename) → **compliant**. Caution: the `NNV_EPS_SHRINK` env hook must be a NO-OP (unset) in submission. Also: a group may submit multiple tools only if they "differ substantially … different parameters with the same tool are not permitted" → submit **one** configured NNV.

### Prepare-vs-run split
- `prepare_instance.sh`: capped at **10 minutes**, **"should not do any analysis"** (exceeding → that instance scored unknown). Called per instance; per-benchmark prep allowed; must not do verification or per-instance tuning. (Manifest/ONNX conversion is setup, not analysis → allowed here.)
- `run_instance.sh`: the **timed** phase (timeout = arg 6, seconds).
- Instances run one-by-one, alternating prepare→run. A whole category can be skipped by having `prepare_instance.sh` return nonzero.

### Tool submission interface
- `install_tool.sh v1` — once per image; deps/licenses/compile (manual license step allowed).
- `prepare_instance.sh v1 <category> <onnx> <vnnlib>` — ≤10 min, no analysis.
- `run_instance.sh v1 <category> <onnx> <vnnlib> <results_file> <timeout_seconds>`.

### Result-file format (⚠️ INTERNAL INCONSISTENCY IN THE RULES)
- "Input and Output Formats" section + the witness example: `{"sat","unsat","timeout","error","unknown"}`.
- The `run_instance.sh` script section says: `holds, violated, timeout, error, or unknown` (`holds`=`unsat`, `violated`=`sat`).
- **NNV writes `sat`/`unsat`/`unknown`** (`run_vnncomp_instance.m:361-377`), matching the Output-Formats section + the witness example, and this scored 697.3 pts in 2025 with the same runner. **Action: post a one-line clarifying question on the 2026 rules issue; keep `sat`/`unsat` unless told otherwise.**

### SAT witness format + tolerance
- On `sat`, line 2+ = `((X_0 val)(X_1 val)…(Y_0 val)…)`. No penalty if onnxruntime replay matches **< 1e-3 relative** on outputs **and** constraints met to **1e-4 absolute**. Missing witness on `sat` → penalty.
- NNV's `write_counterexample` (`run_vnncomp_instance.m:1073-1091`) emits exactly this shape at `%.16g` — **format matches**.

### Eval environment
- CPU `m5.16xlarge` (64 vCPU, 256 GB) — **recommended for NNV** (CPU/LP/star-based; large nets need the 256 GB the GPU boxes lack). GPU options: `p3.2xlarge` (V100, 61 GB), `g5.8xlarge` (A10G, 128 GB). **Pick ONE platform for all benchmarks.**
- Per-benchmark total runtime ≤ **6 hours**; per-instance timeout set by the benchmark proposer.

### Timeline
| Milestone | Date |
|---|---|
| Registration / rules / benchmark submit / voting | Mar 24 / May 1 / May 8 / May 29, 2026 (passed) |
| **Tool submission window** | **June 30, 2026** |
| Competition @ SAIV/FLoC, Lisbon | **July 24–25, 2026** |

Taylor Johnson is a General Chair → NNV can win categories/certificates but cannot receive cash prizes (organizer rule); honor-code expectations apply.

---

## (b) Gap analysis — `VNN_COMP2026/` ready vs. needs-updating

`VNN_COMP2026/` has the 3 shell scripts + README; the MATLAB logic is deliberately NOT forked — `run_instance.sh` calls `../VNN_COMP2025/execute.py` → `run_vnncomp_instance.m` (single source of truth). Sound design.

| Item | Status |
|---|---|
| MATLAB release R2024b→**R2026a** (`install_tool.sh`) | ✅ Updated |
| install_tool.sh: ONNX **+ PyTorch** converter | ✅ Updated (PyTorch needed by manifest importer) |
| install_tool.sh: warm NNV install (tbxmanager) | ✅ New — addresses 14.2 s startup |
| install_tool.sh / run_instance.sh: relative paths (vs hardcoded `/home/ubuntu/toolkit`) | ✅ Updated |
| prepare_instance.sh: kill stale procs, quieted, no analysis | ✅ Compliant |
| **`execute.py` hardcoded `/home/ubuntu/toolkit/code/nnv/` path (`execute.py:55`)** | 🔴 **STALE** — all 3 scripts depend on it; breaks if cloned elsewhere. **[FIXED 2026-06-11: now derives NNV root relatively, with `$NNV_ROOT` override.]** |
| Result string `sat`/`unsat`/`unknown` | ✅ (confirm vs `holds/violated` rules inconsistency) |
| onnxruntime witness guard | 🟡 Called in runner (`:121-131`); 3-state unavailable-case fallback to verify |
| `.vnnlib.gz` / VNN-LIB 2.0 | 🟡 Not implemented (only 4 extended-track 2.0-only benchmarks) |
| Sweep vs the **34-folder 2026 benchmark repo** | 🔴 Not yet swept (current sweeps used the 2025 set) |

**Confirmed 2026 benchmark set** (live, `VNN-COMP/vnncomp2026_benchmarks/benchmarks`, 34 folders): new 2026 entries incl. `adaptive_cruise_control_non_linear_2026`, `cgan2026`, `challenging_certified_training_2026`, `isomorphic_acasxu_2026`, `monotonic_acasxu_2026`, `relusplitter_2026`, `smart_turn_multimodal_2026`, `soundnessbench_2026`, plus duplicated old/new folders.

**Net:** shell scaffolding largely updated; the two items that can break a real run were (1) the hardcoded `execute.py` path (now fixed) and (2) no validation sweep against the actual 2026 repo yet.

---

## (c) Lessons from the 2025 retrospective worth carrying forward

1. **The gap is executional, not theoretical.** CORA (954.6) and nnenum (740.8) use the same reachability family and beat NNV (697.3). Levers: **SAT-finding > penalty leakage > speed > UNSAT precision.**
2. **#1 weakness — SAT-finding: 354 vs ~1000.** → gradient PGD/FGSM falsifier (done).
3. **19 incorrect/missing-CE leaks (−150 each)** — almost certainly invalid/missing witnesses → `validate_witness` re-eval + onnxruntime guard. **Highest-value carry-forward: one −150 erases 15 correct answers.**
4. Worst startup (14.2 s) → heavy setup once in install_tool.sh; minimal prepare.
5. **Never run `exact-star` on large/compute-bound nets** (safenlp/sat_relu timeouts) → `approx-zono`→`abs-dom` ladder.
6. Four `reachOptionsList{1}`-overwrite bugs (dist_shift/linearizenn/cora/malbeware) ran the wrong/exponential method — fixed (#309).
7. 4 compute-bound nets (cifar100, tinyimagenet, vggnet16; cora tractable with budget) are scalability limits — budget/cloud, not capability.
8. **Always pair onnx↔vnnlib via `instances.csv`**, never naive folder pairing.
9. **Soundness held end-to-end in 2025: 0 unsound verdicts.** Keep the gate (no not-robust from over-approximation).

---

## (d) Prioritized checklist (plan — execute deliberately)

**🔴 P0 — blocking correctness/interface**
1. ✅ **Fix the hardcoded path in `execute.py`** — DONE 2026-06-11 (relative NNV root + `$NNV_ROOT` override).
2. **Confirm result-string wording** (`sat`/`unsat` vs `holds`/`violated`) on the 2026 rules issue; keep `sat`/`unsat`.
3. **Verify the 3 scripts are committed executable** (`chmod +x`); `set -e` won't abort on benign `|| true`.

**🔴 P1 — score levers (mostly merged; verify on the eval path)**
4. **Run the full cloud sweep against the actual `vnncomp2026_benchmarks` repo** (34 folders), pairing via each `instances.csv`, realistic per-benchmark timeouts — re-point `benchmarks_repo` in `.github/workflows/vnncomp-sweep.yml`.
5. **Wire the onnxruntime witness guard with a 3-state fallback** (violated / not-violated / onnxruntime-unavailable) — downgrade to `unknown` ONLY on confirmed non-violation; keep `sat` when onnxruntime is unavailable.
6. **Confirm `install_tool.sh` provisions Python + onnxruntime** (the witness guard + the manifest importer need it on the eval VM).

**🟠 P2 — early risk surfacing**
7. **Unofficial TUM dry-run** (`https://vnn.repeatability.cps.cit.tum.de/`) ASAP — surfaces install/license/format errors (the top schedule risk, not verification).
8. **Re-validate per-category tuning** against renamed 2026 folders (`relusplitter_2026`, `cgan2026`, `soundnessbench_2026`); confirm `cctsdb_yolo` matches before `yolo`.
9. **ERROR-category triage** (cctsdb_yolo, cgan, collins_aerospace, soundnessbench): prefer the xval-gated **manifest path** (sound: worst case `unknown`, never −150) over a risky `matlab2nnv` importer fix — especially `soundnessbench`, purpose-built to catch unsound verifiers.

**🟡 P3 — stretch**
10. VNN-LIB 2.0 Python bridge for the 4 extended-track 2.0-only benchmarks.
11. Importer `ReshapeLayer .Vars` fix (coverage, not score) — or route via manifest.
12. CEGAR / bandit method ordering / GPU PGD — future.

---

## (e) Risks

1. **🔴 MATLAB R2026a licensing on the eval VM** — NNV needs commercial MATLAB + DL Toolbox (+ ONNX/PyTorch converters). Rules permit a manual license step, but a license file/server must actually be provisioned or **every instance errors (0 pts across the board)**. The single biggest existential risk — confirm the MathWorks competition license covers R2026a on `m5.16xlarge` before the dry-run.
2. **🔴 Hardcoded `execute.py` path** — FIXED.
3. **🟠 Python-importer (manifest) dependency** — `traffic_signs`, `lsnc_relu`, Tier-D (and any newly manifest-routed cgan/soundnessbench) need the vendored `onnx2nnv_python` env (Python + onnx + onnxruntime + numpy + scipy). Pin versions in `install_tool.sh`. (Manifest generation measured cheap: <3 s/model.)
4. **🟠 Result-string ambiguity** — low harm risk (2025 precedent) but confirm.
5. **🟠 Startup overhead** — `execute.py` starts a MATLAB engine per `run_instance`; the 14.2 s penalty may persist on short-timeout benchmarks. Consider a persistent/shared engine if timing data shows score bleed.
6. **🟢 GPU policy** — not a risk on `m5.16xlarge` CPU (recommended).
7. **🟢 −150 discipline** — every new importer fix is a fresh −150 hazard; keep `unknown` the default for any unvalidated verdict/witness path (esp. soundnessbench, collins_aerospace).

---

### Key file paths
- 2026 scripts: `code/nnv/examples/Submission/VNN_COMP2026/{install_tool.sh, prepare_instance.sh, run_instance.sh, README.md}`
- Maintained runner (single source of truth): `code/nnv/examples/Submission/VNN_COMP2025/run_vnncomp_instance.m` (result write `:361-377`, witness `:1073-1091`, onnx guard `:121-131`)
- Bridge (path fixed): `code/nnv/examples/Submission/VNN_COMP2025/execute.py`
- Planning docs: `VNNCOMP2026_{READINESS,STRATEGY,PERCATEGORY_TUNING,PROGRESS}.md`, `VNNCOMP2025_{RESULTS_ANALYSIS,STATUS}.md`
- Live rules: `https://github.com/VNN-COMP/vnncomp2026/blob/main/rules.md`; benchmarks: `https://github.com/VNN-COMP/vnncomp2026_benchmarks` (34 folders); eval site: `https://vnn.repeatability.cps.cit.tum.de/`

**Bottom line:** The submission is in good shape — most high-ROI engineering is merged, per-benchmark tuning is rules-compliant, and the soundness posture fits the −150 rule. The remaining work is **operational, not algorithmic**: the `execute.py` path (fixed), confirm R2026a licensing + the Python importer env on the eval VM, resolve `sat/unsat` wording, and run a TUM dry-run + a full **2026**-benchmark sweep before June 30.
