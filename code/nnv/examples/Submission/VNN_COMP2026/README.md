# NNV — VNN-COMP 2026 submission scripts

The three-script contract for the VNN-COMP 2026 measurement harness (TUM repeatability
site). Target image: Ubuntu + **MATLAB R2026a**. Recommended platform: **`m5.16xlarge`
CPU** (NNV is CPU/LP/star-based; the large ViT/VGG/TinyImageNet nets need the RAM the GPU
box lacks). See [`../../../../VNNCOMP2026_STRATEGY.md`](../../../../VNNCOMP2026_STRATEGY.md).

- **`install_tool.sh v1`** — once per image: installs the ONNX/PyTorch converters (mpm)
  and runs NNV `install` so the tbxmanager toolboxes are warm (keeps per-instance
  overhead low — NNV had the worst startup, 14.2 s, in 2025).
- **`prepare_instance.sh v1 <cat> <onnx> <vnnlib>`** — fast state cleanup only, no analysis.
- **`run_instance.sh v1 <cat> <onnx> <vnnlib> <results> <timeout>`** — runs one instance
  and writes a single-word result + a **validated** SAT witness.

## Self-contained 2026 submission

This folder owns the **full** MATLAB verification logic — the maintained runner
`run_vnncomp_instance.m` (gradient **PGD/FGSM falsifier** + **witness validation** — strategy
Pillars 1 & 2), its `execute.py` bridge, the falsifiers, and the sweep/first-pass tooling, all
local. `VNN_COMP2025` is frozen at its **original 2025 participation state** (the `vnncomp2025`
tag) for reproducibility; all 2026 development happens here.

## Pre-submission checklist (deadline June 30, 2026 AoE)

1. **Dry-run on the TUM site early** to surface install/format errors (the schedule risk is
   install + witness format, not verification).
2. **Witnesses:** every `sat` is re-validated (`validate_witness.m`) before it is written — this
   closes the 19 incorrect/missing-CE (−150) leaks from 2025.
3. **Regular track (24 benchmarks, all VNN-LIB 1.0)** runs on the existing parser — no 2.0 needed.
4. **Extended track (4 VNN-LIB-2.0-only):** add the `vnnlib`-pip → `.mat` bridge (readiness plan)
   if time allows; otherwise skip the 2.0-only benchmarks (regular track is where the score is).
5. Make the scripts executable: `chmod +x *.sh`.
