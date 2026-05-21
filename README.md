# NNV

**Neural Network Verification Toolbox for MATLAB**

NNV is an open-source MATLAB toolbox for formal verification of deep neural networks and learning-enabled cyber-physical systems. It implements set-based reachability analysis using Star sets, ImageStars, VolumeStars, GraphStars, and more — supporting verification of feedforward, convolutional, recurrent, graph, and segmentation networks, as well as neural ODEs and neural network control systems.

Developed by the [VeriVITAL](https://www.verivital.com/) research group at Vanderbilt University.

## ATVA 2026 Artifact — NNV 3.0 Tool Paper

This repository includes the artifact-evaluation package for the **NNV 3.0
tool paper** at ATVA 2026. The reproducibility flow exercises six
experiments (FairNNV, ProbVer, GNNV, VideoStar, ModelStar, and a
NNV-vs-AIVL ToolComparison) and renders the paper's headline tables.

### Reviewer paths

| Path | When to choose |
|---|---|
| **Code Ocean** (recommended) | You want one-click reproduction with no local setup. AIVL and all toolboxes are pre-installed in the capsule at - ADD - |
| **Local Docker** | You have Docker installed and access to a MATLAB license. Build a self-contained R2025b image from this repo. |
| **Local MATLAB** | You already have MATLAB R2025b on the host and want to run experiments outside Docker. See [`code/nnv/examples/NNV3.0/README.md`](code/nnv/examples/NNV3.0/README.md). |

The walkthrough below covers the **Local Docker** path. Code Ocean
reviewers can skip to the capsule link above.

### AIVL Support Package

The MathWorks AI Verification Library is non-redistributable, so the
tarball is not in this repository. Three acquisition paths are supported:

1. **MATLAB Add-On Explorer** — any MATLAB user can install AIVL via
   **Home → Add-Ons → Get Add-Ons** ("AI Verification Library").
2. **ATVA 2026 AE reviewers** — a pre-built tarball is provided via a
   private Dropbox link in the HotCRP submission cover note.

Without AIVL, the NNV-side of the ToolComparison still runs end-to-end;
only the MathWorks rows are absent from the comparison table. NNV
claims in the paper are unaffected.

### Local Docker — step by step

> **Windows / PowerShell users.** The snippets below use bash syntax
> (`$PWD`, trailing `\` line continuations, single-quoted heredocs).
> To run them directly in native PowerShell, apply three substitutions:
>
> - Replace `$PWD` with `${PWD}` (curly braces) for clean expansion.
> - Quote each `-v` mount value as a single string, e.g.
>   `-v "${PWD}/results:/out"` instead of `-v "$PWD/results":/out`.
>   Without this the Windows drive-letter colon collides with the
>   `source:dest` separator and Docker fails with
>   `invalid reference format`.
> - Replace trailing `\` line continuations with backtick (`` ` ``),
>   or simply put each command on one line.
>
> Reviewers who prefer to keep the snippets verbatim can run them
> unchanged in **WSL2** or **Git Bash**.

**Step 1.** Clone the repository (or unzip the artifact ZIP):

```bash
git clone --recursive https://github.com/verivital/nnv.git
cd nnv
```

The `--recursive` flag pulls in the `npy-matlab` submodule used by
VideoStar. If you already cloned without it, run
`git submodule update --init --recursive`.

**Step 2.** (Optional, AIVL only) Stage the AIVL tarball if you want
the MathWorks-side rows of the ToolComparison experiment. Drop the
file at:

```
code/nnv/examples/NNV3.0/ToolComparison/utils/atva26-aivl.tar.gz
```

If you skip this step, the Dockerfile prints a warning during build and
the ToolComparison runs NNV-only.

**Step 3.** Build the Docker image (~15-25 min depending on network).
Pass your MATLAB license source via `--build-arg`:

```bash
docker build -t nnv3.0 --build-arg LICENSE_SERVER=<port>@<host> .
```

The build does not validate the license — it is consumed at first
`matlab` invocation inside the container. If you have a `network.lic`
file rather than a server, mount it at run time with
`-v /path/to/network.lic:/opt/matlab/R2025b/licenses/network.lic`
instead of `--build-arg`.

If you don't know your `port@host`, the helper scripts at
[`code/nnv/examples/NNV3.0/utils/`](code/nnv/examples/NNV3.0/utils/)
search the standard MATLAB license locations on your host and emit
any `port@host` values they find (`find-matlab-license.ps1` on Windows,
`find-matlab-license.sh` on Linux / macOS / WSL / Git Bash).

**Step 4.** (Optional) GPU sanity check. If you intend to run with
`--gpus all`, verify Docker can actually see your GPU before
committing to the ~30 min smoke run:

```bash
docker run --rm --gpus all nnv3.0 nvidia-smi
```

This should print your device name and driver version. If it errors
with `could not select device driver`, install/enable the NVIDIA
Container Toolkit on the host (or on Docker Desktop with WSL2,
enable GPU support in *Settings → Resources → WSL Integration*).

**Step 5.** Run the smoke test (~30 min on a 4-core / 16 GB host).
Mount a host directory at run time so the per-experiment outputs and
`summary.csv` persist after `--rm` removes the container. Pre-create
`results/` with permissive permissions so the in-container `matlab`
user can write to it:

```bash
mkdir -p results && chmod 777 results
docker run --rm --gpus all -v "$PWD/results":/out nnv3.0 \
    bash -c "bash run_all.sh && \
        cp -r /home/matlab/nnv/code/nnv/examples/NNV3.0/repeatability_logs /out/"
```

`run_all.sh` runs all six experiments sequentially, each in its own
MATLAB session so a failure in one doesn't lose the others. Drop
`--gpus all` if you don't have an NVIDIA GPU available to Docker:
ProbVer auto-skips on CPU-only hosts (it requires GPU + 48 GB RAM),
and the rest fall back to CPU. Two output artifacts land in
`results/repeatability_logs/`: `summary.csv` (per-experiment wall-clock
and status — expect 5 `ok` rows and `probver,skipped` on CPU, or 6 `ok`
rows with a GPU) and `run.log` (consolidated terminal output filtered
to status markers, per-instance verdicts, and final result tables; set
`NNV3_VERBOSE=1` to disable the filter for debugging).

**Step 6.** (Optional) Full reproduction (~5-7 h) renders the paper's
Tables 5, 6, and 7:

```bash
docker run --rm --gpus all -v "$PWD/results":/out \
    -e TOOLCOMPARISON_MODE=full nnv3.0 \
    bash -c "bash run_all.sh && \
        cp -r /home/matlab/nnv/code/nnv/examples/NNV3.0/repeatability_logs /out/"
```

### Full options

For per-experiment knobs (`NNV3_SKIP`, `NNV3_FORCE_GPU`,
`PROBVER_NUM_SAMPLES`, etc.), expected outputs cell-by-cell, host-side
MATLAB setup, troubleshooting, and reference timings, see:

- [`code/nnv/examples/NNV3.0/README.md`](code/nnv/examples/NNV3.0/README.md) — artifact-level run instructions
- [`code/nnv/examples/NNV3.0/EXPECTED_RESULTS.md`](code/nnv/examples/NNV3.0/EXPECTED_RESULTS.md) — per-table expected verdicts and timings
- [`code/nnv/examples/NNV3.0/ToolComparison/README.md`](code/nnv/examples/NNV3.0/ToolComparison/README.md) — NNV-vs-AIVL benchmark catalog

---

The rest of this README documents NNV as a general toolbox. For artifact
evaluation, the section above and the artifact-specific docs linked
above are the relevant entry points.

## Documentation

**Full documentation, tutorials, and API reference:**

**[verivital.github.io/nnv](https://verivital.github.io/nnv/)**

The documentation site includes:
- [Installation Guide](https://verivital.github.io/nnv/getting-started/installation.html) — setup on all platforms
- [User Guide](https://verivital.github.io/nnv/user-guide/index.html) — architectures, set types, verification methods, ONNX/VNNLIB, LP solvers
- [Examples & Tutorials](https://verivital.github.io/nnv/examples/index.html) — 10+ end-to-end walkthroughs across multiple domains
- [Developer Guide](https://verivital.github.io/nnv/developer/index.html) — extending NNV with new layers, set types, and methods
- [API Reference](https://verivital.github.io/nnv/api/index.html) — function signatures for all classes and utilities
- [Application Domains](https://verivital.github.io/nnv/application-domains.html) — aerospace, automotive, medical imaging, power systems, cybersecurity, fairness
- [Theoretical Foundations](https://verivital.github.io/nnv/theory/index.html) — Star reachability, ImageStar/VolumeStar, GraphStar, ModelStar, probabilistic verification, fairness

## Quick Start

```bash
git clone --recursive https://github.com/verivital/nnv.git
```

In MATLAB:

```matlab
cd('nnv/code/nnv')
install                  % Add NNV and dependencies to path
check_nnv_setup()        % Verify installation

% Your first verification
net = matlab2nnv(trained_network);
input_set = Star(lb, ub);
reachOptions.reachMethod = 'approx-star';
output_sets = net.reach(input_set, reachOptions);
result = net.verify_robustness(input_set, reachOptions, target);
% result: 1 = robust, 0 = not robust, 2 = unknown
```

**Requirements:** MATLAB 2023a+ with Computer Vision, Control Systems, Deep Learning, Image Processing, Optimization, Parallel Computing, Statistics and Machine Learning, Symbolic Math, and System Identification toolboxes.

See the [full installation guide](https://verivital.github.io/nnv/getting-started/installation.html) for detailed platform-specific instructions and optional packages.

## What's New in NNV3

NNV3 introduces major new verification capabilities:

- **VolumeStar**: Verification of video and 3D volumetric data (medical images)
- **GNNV (GraphStar)**: Formal verification of graph neural networks (GCN, GINE) for power systems
- **ModelStar**: Verification under weight perturbations (quantization, hardware errors)
- **FairNNV**: Formal fairness certification (counterfactual and individual fairness)
- **Probabilistic Verification**: Scalable conformal prediction-based analysis for intractable networks
- **Time-Dependent Networks**: Variable-length time series verification
- **Malware Detection Benchmark**: New cybersecurity verification domain

See [README_NNV3_CONTRIBUTIONS.md](README_NNV3_CONTRIBUTIONS.md) for technical details, or the [NNV3 Summary](https://verivital.github.io/nnv/nnv3.html) on the documentation site.

### NNV vs Other Tools

| Architecture | NNV3 | α,β-Crown | CORA | Marabou | nnenum |
|:------------|:----:|:---------:|:----:|:-------:|:------:|
| FFNN / CNN   | ✓ | ✓ | ✓ | ✓ | ✓ |
| RNN / LSTM   | ✓ | ✓ | — | — | — |
| Semantic Seg.| ✓ | ~ | — | — | — |
| GNN          | ✓ | — | ✓ | ~ | — |
| 3D / Video   | ✓ | — | — | — | — |
| Neural ODE   | ✓ | — | — | — | — |
| NNCS         | ✓ | ~ | ✓ | — | — |
| Weight Pert. | ✓ | — | — | — | — |
| Fairness     | ✓ | — | — | — | — |

See the [full 14-tool comparison](https://verivital.github.io/nnv/user-guide/architectures.html#nnv-vs-other-tools) on the documentation site.

## Execution Without Installation

NNV can be run online without installing MATLAB:

| Version | Platform | Link |
|---------|----------|------|
| **NNV3** (latest) | CodeOcean | https://codeocean.com/capsule/6810863/tree/v1 |
| NNV 2.0 (CAV 2023) | CodeOcean | https://doi.org/10.24433/CO.0803700.v1 |
| NNV 1.0 (CAV 2020) | CodeOcean | https://doi.org/10.24433/CO.0221760.v1 |
| Tutorials | MATLAB Online | [Try on MATLAB Online](https://matlab.mathworks.com/) |

## Related Tools

NNV integrates with [NNMT](https://github.com/verivital/nnmt) (model transformation), [HyST](https://github.com/verivital/hyst) (hybrid systems), and [CORA](https://github.com/TUMcps/CORA) (continuous reachability).

## Citation

If you use NNV in your research, please cite:

```bibtex
@inproceedings{nnv2_cav2023,
  author = {Lopez, Diego Manzanas and Choi, Sung Woo and Tran, Hoang-Dung and Johnson, Taylor T.},
  title = {NNV 2.0: The Neural Network Verification Tool},
  year = {2023},
  isbn = {978-3-031-37702-0},
  publisher = {Springer-Verlag},
  address = {Berlin, Heidelberg},
  url = {https://doi.org/10.1007/978-3-031-37703-7_19},
  doi = {10.1007/978-3-031-37703-7_19},
  abstract = {This manuscript presents the updated version of the Neural Network Verification (NNV) tool. NNV is a formal verification software tool for deep learning models and cyber-physical systems with neural network components. NNV was first introduced as a verification framework for feedforward and convolutional neural networks, as well as for neural network control systems. Since then, numerous works have made significant improvements in the verification of new deep learning models, as well as tackling some of the scalability issues that may arise when verifying complex models. In this new version of NNV, we introduce verification support for multiple deep learning models, including neural ordinary differential equations, semantic segmentation networks and recurrent neural networks, as well as a collection of reachability methods that aim to reduce the computation cost of reachability analysis of complex neural networks. We have also added direct support for standard input verification formats in the community such as VNNLIB (verification properties), and ONNX (neural networks) formats. We present a collection of experiments in which NNV verifies safety and robustness properties of feedforward, convolutional, semantic segmentation and recurrent neural networks, as well as neural ordinary differential equations and neural network control systems. Furthermore, we demonstrate the capabilities of NNV against a commercially available product in a collection of benchmarks from control systems, semantic segmentation, image classification, and time-series data.},
  booktitle = {Computer Aided Verification: 35th International Conference, CAV 2023, Paris, France, July 17--22, 2023, Proceedings, Part II},
  pages = {397--412},
  numpages = {16},
  keywords = {neural networks, cyber-physical systems, verification, tool},
  location = {Paris, France}
}
```

```bibtex
@inproceedings{nnv_cav2020,
  author = {Tran, Hoang-Dung and Yang, Xiaodong and Manzanas Lopez, Diego and Musau, Patrick and Nguyen, Luan Viet and Xiang, Weiming and Bak, Stanley and Johnson, Taylor T.},
  title = {NNV: The Neural Network Verification Tool for Deep Neural Networks and Learning-Enabled Cyber-Physical Systems},
  year = {2020},
  isbn = {978-3-030-53287-1},
  publisher = {Springer-Verlag},
  address = {Berlin, Heidelberg},
  url = {https://doi.org/10.1007/978-3-030-53288-8_1},
  doi = {10.1007/978-3-030-53288-8_1},
  abstract = {This paper presents the Neural Network Verification (NNV) software tool, a set-based verification framework for deep neural networks (DNNs) and learning-enabled cyber-physical systems (CPS). The crux of NNV is a collection of reachability algorithms that make use of a variety of set representations, such as polyhedra, star sets, zonotopes, and abstract-domain representations. NNV supports both exact (sound and complete) and over-approximate (sound) reachability algorithms for verifying safety and robustness properties of feed-forward neural networks (FFNNs) with various activation functions. For learning-enabled CPS, such as closed-loop control systems incorporating neural networks, NNV provides exact and over-approximate reachability analysis schemes for linear plant models and FFNN controllers with piecewise-linear activation functions, such as ReLUs. For similar neural network control systems (NNCS) that instead have nonlinear plant models, NNV supports over-approximate analysis by combining the star set analysis used for FFNN controllers with zonotope-based analysis for nonlinear plant dynamics building on CORA. We evaluate NNV using two real-world case studies: the first is safety verification of ACAS Xu networks, and the second deals with the safety verification of a deep learning-based adaptive cruise control system.},
  booktitle = {Computer Aided Verification: 32nd International Conference, CAV 2020, Los Angeles, CA, USA, July 21--24, 2020, Proceedings, Part I},
  pages = {3--17},
  numpages = {15},
  keywords = {Autonomy, Verification, Cyber-physical systems, Machine learning, Neural networks},
  location = {Los Angeles, CA, USA}
}
```

NNV3 paper is in preparation. For feature-specific citations (FairNNV, VolumeStar, GNNV, ModelStar, probabilistic verification), see the [How to Cite](https://verivital.github.io/nnv/publications/citing.html) page.

For the full list of 30+ publications using NNV, see the [Publications](https://verivital.github.io/nnv/publications/papers.html) page.

## Contributors

* [Diego Manzanas Lopez](https://mldiego.github.io/)
* [Hoang-Dung Tran](https://sites.google.com/view/v2a2/)
* [Taylor T. Johnson](http://www.taylortjohnson.com)
* [Neelanjana Pal](https://scholar.google.com/citations?user=3j_f-ewAAAAJ&hl=en)
* [Anne Tumlin](https://atumlin.github.io/)
* [Samuel Sasaki](https://sammsaski.github.io/)
* [Ben Wooding](https://woodingben.com)
* [Xiaodong Yang](https://scholar.google.com/citations?user=xe3Jr7EAAAAJ&hl=en)
* [Patrick Musau](https://pmusau17.github.io/)
* [Sung Woo Choi](https://scholar.google.com/citations?user=choi_sungwoo)
* [Stanley Bak](http://stanleybak.com/)
* [Weiming Xiang](https://xiangweiming.github.io/)
* [Luan Viet Nguyen](https://luanvietnguyen.github.io)

## Acknowledgements

This work is supported in part by AFOSR, DARPA, NSF.

## Contact

For any questions related to NNV, please add them to the [issues](https://github.com/verivital/nnv/issues) or contact [Taylor T. Johnson](mailto:taylor.johnson@vanderbilt.edu), [Samuel Sasaki](mailto:samuel.sasaki@vanderbilt.edu), [Anne Tumlin](mailto:anne.m.tumlin@vanderbilt.edu), or [Ben Wooding](mailto:ben.wooding@vanderbilt.edu).

Past contacts: [Diego Manzanas Lopez](mailto:diego.manzanas.lopez@vanderbilt.edu) and [Hoang-Dung Tran](mailto:trhoangdung@gmail.com).
