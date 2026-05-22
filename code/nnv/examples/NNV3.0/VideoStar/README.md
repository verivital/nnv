# VideoStar — Video Classification Verification (ZoomIn-4f)

Reachability-based robustness verification of a 3D-CNN video
classifier trained on the **MNIST-Video ZoomIn-4f** dataset (4-frame
"zoom in" sequences over MNIST digits). For each test sample and each
ε ∈ {1/255, 2/255, 3/255}, the verifier checks whether all
`L∞`-bounded perturbations of the input video preserve the predicted
class.

The technique extends NNV's Star and ImageStar set representations to
**VolumeStar** (a.k.a. VideoStar), supporting reachability through 3D
convolutional and pooling layers.

This folder is a curated subset of the broader video-stars benchmark
suite.

## References

- **VolumeStar / VideoStar (this work)**: Sasaki, S., Manzanas Lopez,
  D., Robinette, P.K., Johnson, T.T. *Robustness verification of video
  classification neural networks.* IEEE/ACM 13th International
  Conference on Formal Methods in Software Engineering (FormaliSE),
  2025.

## Layout

```
NNV3.0/VideoStar/
├── README.md
├── run_zoomin_4f.m              Canonical reference implementation (MATLAB)
├── run_videostar_zoomin4f.sh    Thin shell wrapper (calls the .m via -batch)
├── run_zoomin_4f.py             Optional Python entry point (matlab.engine)
├── requirements.txt             Python deps for run_zoomin_4f.py
├── models/
│   └── zoomin_4f.onnx           3D-CNN classifier (4-frame variant)
├── data/ZoomIn/
│   ├── mnistvideo_zoom_in_4f_test_data_seq.npy
│   └── mnistvideo_zoom_in_test_labels_seq.npy
├── npy-matlab/                  Vendored .npy reader (.m functions)
├── src/vvn/                     Vendored vvn helpers (verifyvideo.m, etc.)
└── results/                     Timestamped output (<yymmdd-HHMMSS>/)
```

## Entry points

Three ways to launch the same MATLAB pipeline. Use the one that fits
your workflow:

1. **MATLAB (canonical)** — `run_zoomin_4f.m` is the reference
   implementation; the other two entry points wrap it.
2. **Shell wrapper** — `run_videostar_zoomin4f.sh` calls
   `matlab -batch run_zoomin_4f`.
3. **Python (optional)** — `run_zoomin_4f.py` drives the same MATLAB
   code via `matlab.engine`. Adds CLI flags. Requires the MATLAB
   Engine for Python.

## Running

NNV must be installed (`code/nnv/install.m`) and on the MATLAB path
inside the Docker image.

### Default sweep (paper-aligned)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/VideoStar; run_zoomin_4f"
```

10 samples × 3 ε × `relax` algorithm × 30-min timeout per sample.

### Smoke (single sample, single ε)

```matlab
matlab -batch "cd code/nnv/examples/NNV3.0/VideoStar; \
    config.sampleIndices = 1:1; config.epsilon = [1/255]; \
    config.timeout = 300; run_zoomin_4f"
```

### Python entry point

```bash
cd code/nnv/examples/NNV3.0/VideoStar
pip install -r requirements.txt          # one-time
python run_zoomin_4f.py --algorithm relax --num-samples 10
```

## Configuration parameters

Edit the `CONFIGURATION` section near the top of
[`run_zoomin_4f.m`](run_zoomin_4f.m), or set fields in the `config`
struct before invoking the runner (the runner has a config-guard):

| Field                  | Default            | Notes |
|------------------------|--------------------|-------|
| `config.dsType`        | `'zoom_in'`        | `'zoom_in'` or `'zoom_out'` |
| `config.sampleLen`     | `4`                | Frames per video — only `4` ships in this folder |
| `config.verAlgorithm`  | `'relax'`          | `'relax'` or `'approx'` |
| `config.numClasses`    | `10`               | MNIST digit classes |
| `config.epsilon`       | `[1/255; 2/255; 3/255]` | ε values in normalized pixel units |
| `config.timeout`       | `1800`             | Per-sample timeout (s) |
| `config.sampleIndices` | `1:10`             | Which samples to verify |

## Outputs

A timestamped subfolder `results/<yymmdd-HHMMSS>/` is created per run.
Inside, one CSV per ε:

- `eps=1_255.csv`, `eps=2_255.csv`, `eps=3_255.csv`

Columns: `Sample Number, Result, Time, Method`.
Result codes: `1` = verified safe, `0` = violated, `2` = unknown,
`3` = timeout, `-1` = error.

This folder generates **no figures** — the paper's figures are
produced post-hoc from the CSVs.

## Expected runtime

Measured on an NVIDIA RTX 4060 host running the `nnv3.0` Docker image
(MATLAB R2024b), CPU verification:

- **Smoke** (1 sample × ε=1/255): **~30 s**
- **Default sweep** (10 samples × 3 ε × `relax`): **~12 minutes**
