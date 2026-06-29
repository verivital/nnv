#!/bin/bash
# VNN-COMP 2026 GLOBAL competition env -- THE SINGLE SOURCE for the always-on, category-INDEPENDENT
# soundness-policy verification flags. SOURCE this from every harness (run_instance.sh + every sweep
# orchestrator) instead of re-declaring the flags, so they can never drift apart. This file exists because
# of the 2026-06-24 cifar regression: sweep_lambda.sh had silently dropped these flags that run_instance.sh
# and robust_runner.m set, costing cifar100 its 20 unsats (23->2 decided).
#
#   usage:  source "<path>/vnncomp2026_env.sh"
#
# These four are deliberate COMPETITION OPT-INS (NOT code defaults), so the test suite / other callers of
# run_vnncomp_instance.m are unaffected unless they source this. Category-CONDITIONAL flags
# (safenlp NNV_FALSIFY_MAXTIME, cifar100 NNV_BAB_BETA_ITERS/NNV_CONV_BETA_FRONTIER/NNV_AMORT_ALPHA) are
# DELIBERATELY NOT here -- they are set per-category inside run_vnncomp_instance.m (PR #424), so they are
# harness-independent by construction and can't drift either. Per-instance budgets (NNV_REACH_BUDGET,
# NNV_ACAS_BAB_TIMECAP) are derived from the timeout by the caller, also not here.
# Source of truth for the full per-benchmark map: status-repo BENCHMARK_RUN_MATRIX.md.

# NNV_CONV_TRUST_FP32: emit the unsat verdict FROM the GPU-single batched-BaB screen (raw CROWN bound on the
#   GPU, ~5x faster than the FP64-CPU reconfirm), PGD falsify-first as the backstop -- the same two-mechanism
#   model alpha-beta-CROWN / NeuralSAT use. Validated 0 false-robust vs the alpha-beta-CROWN gold set
#   (23 cifar100/tinyimagenet incl gold-SAT). Forces the GPU screen ON. WITHOUT it -> 0 conv unsats.
#   SHIPPING POLICY (Taylor, 2026-06-28): KEEP ON for the competition. It governs ~all of cifar100's 42%
#   and tinyimagenet's 40% (~80 breadth-pts); turning it OFF reverts those to unknown (the sound FP64
#   reconfirm is 120-194s, over the 100s budget). SOUNDNESS basis + the post-competition exploration paths
#   (sound-FP32 emit NNV_SOUND_FP32_TIGHT, FP64-only, matmul measured-delta) are in CONV_TRUST_FP32_POLICY.md.
#   CONFIGURABLE: a CLEAN toggle now -- =1 ON, =0 (or unset) OFF -> straight to the sound FP64 double-confirm.
#   (Was a footgun: `=0` did NOT disable; it read non-empty=ON. Fixed via i_envon() in run_vnncomp_instance.m.)

# NNV_ORT_PYTHON: the dedicated python for execute.py's matlab.engine + the ORT/witness gate + onnx2nnv.
# An explicit value always wins (e.g. the Lambda dev box's taylor_venv). On the eval box, post_install.sh
# installs matlab.engine + the ONNX stack into a SYSTEM python via `sudo --break-system-packages` (system
# site-packages are readable by every user incl. the run user `ubuntu`, and persist). Since which system
# python carries the engine is decided at install time, we PROBE the system pythons at runtime (run as
# ubuntu) and pick the first whose matlab.engine imports. This sidesteps the smoke 269-276 failures: no venv
# (a venv symlinks to a per-user base python the run user cannot reach) and no /home/ubuntu artifact file
# (post_install runs as a user that cannot write /home/ubuntu).
if [ -z "${NNV_ORT_PYTHON:-}" ]; then
    for _c in /usr/bin/python3.13 /usr/bin/python3.12 /usr/bin/python3.11 /usr/bin/python3.10 /usr/bin/python3; do
        [ -x "$_c" ] && "$_c" -c "import matlab.engine" >/dev/null 2>&1 && { NNV_ORT_PYTHON="$_c"; break; }
    done
fi
export NNV_ORT_PYTHON="${NNV_ORT_PYTHON:-python3}"

export NNV_CONV_TRUST_FP32=1
# NNV_CONV_NO_STAR: a non-certifying conv precheck emits a FAST sound 'unknown' (Star never certifies these
#   resnets, it just burns the whole timeout). WITHOUT it -> every conv instance times out.
export NNV_CONV_NO_STAR=1
# NNV_CONV_FRONTIER: BaB nodes per GPU kernel.
export NNV_CONV_FRONTIER=512
# NNV_QUARANTINE_CPSTAR: strip the PROBABILISTIC (conformal) cp-star reach so a cp-star 'unsat' can never be
#   emitted as a sound verdict (cp-star is not an FP64 proof -> a latent -150 under the 0-wrong rule).
#   STRICTLY SOUND: only ever converts a cp-star unsat -> unknown; the FP64 / GPU-BaB / exact paths still
#   decide what they can.
export NNV_QUARANTINE_CPSTAR=1

# ---- PER-INSTANCE PARALLELISM PROFILE (2026-06-29) -----------------------------------------------------
# The COMPETITION runs ONE instance at a time on the whole g5.8xlarge (32 vCPU = 16 PHYSICAL cores + 1 A10G
# 24GB), so each instance should use the entire machine. Our LOCAL/LAMBDA sweeps run N=4 MATLAB sessions on
# one box (1-per-GPU), so the CPU knobs must be divided by N or they oversubscribe (4 sessions x all-cores
# BLAS = ~4x the cores). GPU/conv knobs are profile-INVARIANT: the sweep is already 1-instance-per-GPU
# (CUDA_VISIBLE_DEVICES=$gpu), so each sweep session already owns a full GPU's memory -- same budget as the
# competition. EVERY knob here is PERF-ONLY: serial==parallel verdict, larger batch only ever times-out/OOMs
# to a SOUND 'unknown' (no -150 risk). All use ':-' so an outer export (the sweep harness) wins.
# Switch profiles with NNV_PARALLEL_PROFILE=competition|sweep. Consumed by run_vnncomp_instance.m
# (NNV_BLAS_THREADS -> maxNumCompThreads; NNV_MAX_WORKERS -> the disjunct-level parfor pool size).
export NNV_PARALLEL_PROFILE="${NNV_PARALLEL_PROFILE:-competition}"
if [ "$NNV_PARALLEL_PROFILE" = "sweep" ]; then
    export NNV_MAX_WORKERS="${NNV_MAX_WORKERS:-1}"     # disjunct parfor -> serial (no per-session pool storm)
    export NNV_BLAS_THREADS="${NNV_BLAS_THREADS:-4}"   # ~cores/N; the sweep harness sets the exact value
else
    export NNV_MAX_WORKERS="${NNV_MAX_WORKERS:-auto}"  # auto = min(physical cores, #disjuncts)
    export NNV_BLAS_THREADS="${NNV_BLAS_THREADS:-auto}" # auto = MATLAB default (all physical cores)
fi
# GPU β-refine batch width (cifar100 conv TIER-2). Default 64 (validated on the 11GB sweep GPUs). On the
# competition A10G (24GB, all owned by the one instance) it can likely rise to ~96-128 to cut autodiff
# passes -- left configurable here; CONFIRM ON g5 (no OOM, within budget) before raising. See
# research/parallelism_g5_2026-06-29.md.
export NNV_CONV_BETA_FRONTIER_COMP="${NNV_CONV_BETA_FRONTIER_COMP:-96}"  # candidate g5 value (probe before use)
