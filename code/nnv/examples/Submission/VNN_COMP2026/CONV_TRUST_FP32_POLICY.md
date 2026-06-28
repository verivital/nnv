# NNV_CONV_TRUST_FP32 — conv unsat-emit soundness policy

**Decision (Taylor, 2026-06-28): SHIP WITH IT ON for VNN-COMP 2026.** This doc records what the
policy is, why it is sound enough to ship, exactly what rides on it, how to toggle it, and the
post-competition exploration paths.

## What it does
For conv nets (cifar100, tinyimagenet, vggnet), NNV runs a **GPU-single FP32 batched-BaB "screen"**
(a raw CROWN lower bound on the spec margin, run entirely on the GPU). When the screen certifies
`robust` AND the falsify-first PGD attack found no counterexample, NNV **emits `unsat` directly from
the screen** and skips the ~120–194s FP64-CPU double-confirm. (`run_vnncomp_instance.m`, the
`i_gpu_bab_precheck` conv branch.)

## Why it is sound enough to ship (but is NOT a rigorous FP64 proof)
- **Two-mechanism model.** A `robust` FP32 lower bound + an independent adversarial attack (PGD)
  finding nothing is exactly the soundness model the GPU-winning tools use: **α,β-CROWN runs stock
  FP32 with PGD; NeuralSAT likewise.** We are not doing anything more aggressive than the field SOTA.
- **The FP32 rounding gap is ~1e-6**, far below the real robustness margins on these nets (the
  overnight P2 analysis measured intrinsic margins 1.4e-4 … 0.042; the thin-margin cifar recoveries
  have FP64 margins like +0.042).
- **Validated 0-false-robust** vs the α,β-CROWN gold set (23 cifar100/tinyimagenet incl gold-SAT).
- **Overnight defense-in-depth (2026-06-28):** the 5 thin-margin cifar resnet_large recoveries
  (idx 496/4385/5308/4757/8589) were independently **FP64 double-oracle confirmed** robust (unsat via
  3/7/15/9/175-node double). [Evidence lives in the team's SEPARATE internal status-tracking repo
  (not this nnv repo): `research/cifar_fp64_goldgate_2026-06-27/`.]
- **Residual risk (honest):** it is NOT an FP64 proof. If FP32 rounding ever flips a real margin on an
  instance PGD also misses → a false unsat → −150. The above evidence bounds, but does not eliminate,
  that risk. This is the standard FP32-verifier risk the whole field accepts.

## What rides on it
Essentially **all of cifar100's ~42% and tinyimagenet's ~40%** (≈ **80 breadth-points**) — their
unsats exist only via this path. Turning it OFF reverts them to `unknown`, because the sound FP64
reconfirm (120–194s) exceeds the 100s budget. So the choice is "ship the field-standard FP32 policy"
vs "≈0 conv unsats".

## How to toggle it (configurable, clean)
`NNV_CONV_TRUST_FP32` is now a **clean boolean** (via `i_envon()` in run_vnncomp_instance.m):
- `=1` / `true` / `on` / `yes` → **ON** (competition default; `vnncomp2026_env.sh`).
- `=0` / `false` / anything else / **unset** → **OFF** → conv goes straight to the sound FP64 double-confirm.
- (Historical footgun, fixed 2026-06-28: it used to read `~isempty(getenv(...))`, so `=0` did NOT
  disable — any non-empty value, including the string `'0'`, counted as ON.)

## Post-competition exploration paths (for tightening the soundness story)
1. **Broader FP64 verification (free, no risk):** run the FP64 double-oracle (set `NNV_CONV_TRUST_FP32=0`,
   generous timeout) over ALL conv unsats, not just the 5 — confirm 0 disagreements with the screen.
2. **Sound-FP32 emit (`NNV_SOUND_FP32_TIGHT`, default OFF):** the outward-rounded FP32 path — a provably
   sound lower bound (every CROWN bound widened by a rigorous roundoff bound) emitted from the GPU. The
   blocker is that the worst-case Higham γ_n widening is too loose on the wide-conv matmul; the
   **matmul measured-δ** (running-error analysis) is the lever to make it tight (documented in the
   team's separate internal status-tracking repo, memory `conv-sound-emit-diagnosis`; SUPERVISED-only,
   −150-sensitive).
3. **FP64-only (maximally conservative):** `NNV_CONV_TRUST_FP32=0`. Sound, but loses ~80 breadth-points
   at the 100s budget unless the FP64 confirm is sped up (the matmul-δ / batched-double work).
