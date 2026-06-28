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
# post_install.sh builds $HOME/.nnv_venv (matlab.engine + onnx stack + torch) because the eval box has
# anaconda first on PATH, so a bare `python3` is ambiguous and the MATLAB engine would not install into it
# (smoke 269: every result error_exit_code_1 from ModuleNotFoundError 'matlab'). An explicit NNV_ORT_PYTHON
# still wins (e.g. the Lambda dev box's taylor_venv). post_install records the VERIFIED working python path
# (the uv/3.12 venv, or a system-python fallback) in ~/.nnv_python_path -- read it so we use what actually works.
export NNV_ORT_PYTHON="${NNV_ORT_PYTHON:-$(cat "$HOME/.nnv_python_path" 2>/dev/null || echo "$HOME/.nnv_venv/bin/python")}"

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
