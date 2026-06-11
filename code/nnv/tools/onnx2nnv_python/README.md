# onnx2nnv (Python ONNX → NNV importer)

Standalone Python importer that converts an ONNX model into an NNV manifest
(`<model>.nnv.mat`) consumed by `code/nnv/engine/utils/load_nnv_from_mat.m`
(whose header references this path). Used by the VNN-COMP 2025 manifest path
to import + cross-validate networks the MATLAB `matlab2nnv` path cannot.

## Contents

- `onnx2nnv.py` — the importer: walks the ONNX graph, maps ops to NNV layers,
  and writes the `.nnv.mat` manifest (weights, layer specs, layout flags).
- `xvalidate.py` — a Python-side smoke check: runs the ONNX model under
  `onnxruntime` on N random inputs and prints each output's shape/range. It
  confirms the ONNX loads and runs; it does **not** itself compare against the
  NNV manifest. The authoritative ONNX-vs-NNV forward-pass cross-validation (the
  `xval` numbers in `VNN_COMP2025_SUPPORT.md`) is the **MATLAB-side** `xvalidate.m`,
  which runs `NN.evaluate`. The runner **gates reach on that xval** — a network
  whose forward pass diverges from ONNX is never reached, so a mis-import cannot
  produce an unsound verdict.
- `tests/test_onnx2nnv.py` — importer unit tests.

## Usage (sketch)

```bash
# produce the manifest (output path is an optional POSITIONAL arg; defaults to
# <model>.nnv.mat alongside the ONNX):
python onnx2nnv.py <model.onnx> [<model.nnv.mat>] [--vnnlib <spec.vnnlib>]
# Python ORT smoke check (prints output ranges; --mat defaults alongside the ONNX):
python xvalidate.py <model.onnx> [--mat <model.nnv.mat>] [--n 5] [--seed 0]
```

Then in MATLAB: `net = load_nnv_from_mat('<model.nnv.mat>')`.

## Notes

- Requires Python with `onnx`, `onnxruntime`, `numpy`, `scipy` (for `.mat` I/O).
- Generated manifests (`*.nnv.mat`) live next to the benchmark ONNX files and are
  **not** committed (regenerate with the commands above; `load_manifest_net`
  errors with regeneration instructions when one is missing).
- Provenance: copied into the repo from the workspace `tools/onnx2nnv_python/`
  to make the PR self-contained (per `PR_MERGE_READINESS_PLAN.md` §2).
