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
| **Local Docker — online licence** (recommended) | You have a MathWorks account but no network licence server. Sign in once via the browser; uses `Dockerfile.online`. |
| **Local Docker — network licence** | Your institution provides a `port@host` MATLAB licence server (or you have a 30-day eval network licence). Uses the standard `Dockerfile`. |
| **Local MATLAB** | You already have MATLAB R2025b on the host and want to run experiments outside Docker. See [`code/nnv/examples/NNV3.0/README.md`](code/nnv/examples/NNV3.0/README.md). |

Both Local Docker paths are described below in dedicated sections.

### AIVL Support Package

The MathWorks AI Verification Library is non-redistributable, so the
tarball is not in this repository. Two acquisition paths are supported:

1. **MATLAB Add-On Explorer** — any MATLAB user can install AIVL via
   **Home → Add-Ons → Get Add-Ons** ("AI Verification Library").
2. **ATVA 2026 AE reviewers** — a pre-built tarball is provided via a
   private Dropbox link in the HotCRP submission cover note.

Without AIVL, the NNV-side of the ToolComparison still runs end-to-end;
only the MathWorks rows are absent from the comparison table. NNV
claims in the paper are unaffected.

### Cross-cutting notes (apply to both Local Docker paths)

> **Windows / PowerShell users.** The bash snippets below use `$PWD`
> and trailing `\` line continuations. In native PowerShell:
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

> **Windows + Docker Desktop GPU passthrough.** Every `docker run`
> below uses `--gpus all` as the default (drop it on CPU-only hosts;
> ProbVer auto-skips and GNNV / VideoStar fall back to CPU). On
> Windows this needs the **WSL2 backend** in Docker Desktop *Settings
> → General* AND an **NVIDIA driver on the host** (not inside WSL2 —
> Windows-side `nvidia-smi` should already work). Sanity check from
> the host:
> `docker run --rm --gpus all nvidia/cuda:12.8.2-base-ubuntu22.04 nvidia-smi`.
> If that prints your device, `--gpus all` is wired up. From inside
> the MATLAB session, `gpuDeviceCount` returns `1` or more.

### Local Docker — Option A: Online MATLAB licence

For reviewers with a personal MathWorks account (`matlab.mathworks.com`
sign-in) but no network licence server. Builds on the official
`mathworks/matlab:r2025b` image and authenticates via a one-time
browser sign-in. The activation is cached in named Docker volumes so
subsequent runs reuse it without re-prompting.

Trade-off: experiments must run **inside the browser MATLAB session**
(the online sign-in is not visible to `matlab -batch` processes
spawned by `bash run_all.sh`). Three short `.m` entry points handle
this; see the full walkthrough in
[`code/nnv/examples/NNV3.0/ONLINE_LICENSE.md`](code/nnv/examples/NNV3.0/ONLINE_LICENSE.md).

**Step A1.** Clone (with submodules) and optionally stage the AIVL
tarball at `code/nnv/examples/NNV3.0/ToolComparison/utils/atva26-aivl.tar.gz`.

**Step A2.** Build the online image (~20 min, no licence needed):

```bash
docker build -t nnv3.0-online -f Dockerfile.online .
```

**Step A3.** Launch in browser mode + sign in (one-time, ~5 min):

```bash
docker run -it --rm --gpus all \
    -p 8888:8888 --shm-size=512M \
    -v nnv3-matlab-prefs:/home/matlab/.matlab \
    -v nnv3-matlab-mw:/home/matlab/.MathWorks \
    --name nnv3-setup nnv3.0-online -browser
```

Open <http://localhost:8888>, sign in with your MathWorks account.
In the browser MATLAB Command Window, run the three entry points in
turn:

```matlab
run('/home/matlab/nnv/code/nnv/examples/NNV3.0/setup_online.m')   % ~5 min
run('/home/matlab/nnv/code/nnv/examples/NNV3.0/run_smoke.m')      % ~20-25 min
run('/home/matlab/nnv/code/nnv/examples/NNV3.0/run_full.m')       % ~3-5 h, Tables 5,6,7
```

**Step A4.** Extract results to the host (from a separate shell):

```bash
docker cp nnv3-setup:/home/matlab/nnv/code/nnv/examples/NNV3.0/repeatability_logs .
```

Open `repeatability_logs/PAPER_COMPARISON.md` — that single file is
the entry point a repeatability committee opens to find every output
with explicit pointers to `results/ToolComparison/tables/table_main.{tex,txt}`
(paper Tables 5, 6, 7) and to each experiment's CSV / `.mat` outputs.

### Local Docker — Option B: Network MATLAB licence (port@host)

For reviewers whose institution provides a MATLAB licence server, or
who hold a 30-day MathWorks network-licence trial. Builds on
`mathworks/matlab-deps:R2025b`, installs MATLAB via `mpm` at build
time, and bakes the licence source into the image so `bash run_all.sh`
can run unattended end-to-end.

If you don't know your `port@host`, the helper scripts at
[`code/nnv/examples/NNV3.0/utils/`](code/nnv/examples/NNV3.0/utils/)
search the standard MATLAB licence locations on your host
(`find-matlab-license.ps1` on Windows,
`find-matlab-license.sh` on Linux / macOS / WSL / Git Bash).

**Step B1.** Clone (with submodules) and optionally stage the AIVL
tarball at `code/nnv/examples/NNV3.0/ToolComparison/utils/atva26-aivl.tar.gz`.

**Step B2.** Build the network-licence image (~15-25 min):

```bash
docker build -t nnv3.0 --build-arg LICENSE_SERVER=<port>@<host> .
```

The build does not validate the licence; it is consumed at first
`matlab` invocation. If you have a `network.lic` file rather than a
server, omit `--build-arg` and mount the file at run time:
`-v /path/to/network.lic:/opt/matlab/R2025b/licenses/network.lic`.

**Step B3.** (Optional) GPU sanity check before the long run:

```bash
docker run --rm --gpus all nnv3.0 nvidia-smi
```

**Step B4.** Run the smoke test (~30 min). Mount a host directory so
outputs persist after `--rm`:

```bash
mkdir -p results && chmod 777 results
docker run --rm --gpus all -v "$PWD/results":/out nnv3.0 \
    bash -c "bash run_all.sh && \
        cp -r /home/matlab/nnv/code/nnv/examples/NNV3.0/repeatability_logs /out/"
```

`run_all.sh` runs all six experiments sequentially, each in its own
MATLAB session so a failure in one doesn't lose the others. Two
artefacts land in `results/repeatability_logs/`: `summary.csv`
(per-experiment wall-clock and status — expect 5 `ok` rows and
`probver,skipped` on CPU, or 6 `ok` rows with a GPU) and `run.log`
(consolidated terminal output filtered to status markers, per-instance
verdicts, and final result tables; set `NNV3_VERBOSE=1` to disable the
filter for debugging).

**Step B5.** (Optional) Full reproduction (~5-7 h) renders the paper's
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

- [`code/nnv/examples/NNV3.0/README.md`](code/nnv/examples/NNV3.0/README.md) — artifact-level run instructions (both licence paths)
- [`code/nnv/examples/NNV3.0/ONLINE_LICENSE.md`](code/nnv/examples/NNV3.0/ONLINE_LICENSE.md) — online-licence walkthrough + troubleshooting
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

Earlier NNV versions can be run online without installing MATLAB:

| Version | Platform | Link |
|---------|----------|------|
| NNV 2.0 (CAV 2023) | CodeOcean | https://doi.org/10.24433/CO.0803700.v1 |
| NNV 1.0 (CAV 2020) | CodeOcean | https://doi.org/10.24433/CO.0221760.v1 |
| Tutorials | MATLAB Online | [Try on MATLAB Online](https://matlab.mathworks.com/) |

For NNV 3.0 (the ATVA 2026 artifact), use one of the Local Docker
paths described in the
[ATVA 2026 Artifact](#atva-2026-artifact--nnv-30-tool-paper) section
above.

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
