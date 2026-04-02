# NNV

**Neural Network Verification Toolbox for MATLAB**

NNV is an open-source MATLAB toolbox for formal verification of deep neural networks and learning-enabled cyber-physical systems. It implements set-based reachability analysis using Star sets, ImageStars, VolumeStars, GraphStars, and more — supporting verification of feedforward, convolutional, recurrent, graph, and segmentation networks, as well as neural ODEs and neural network control systems.

Developed by the [VeriVITAL](https://www.verivital.com/) research group at Vanderbilt University.

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
cd nnv
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
  booktitle = {Computer Aided Verification (CAV)},
  year = {2023},
  doi = {10.1007/978-3-031-37703-7_19}
}
```

```bibtex
@inproceedings{nnv_cav2020,
  author = {Tran, Hoang-Dung and Yang, Xiaodong and Manzanas Lopez, Diego and
            Musau, Patrick and Nguyen, Luan Viet and Xiang, Weiming and
            Bak, Stanley and Johnson, Taylor T.},
  title = {NNV: The Neural Network Verification Tool for Deep Neural Networks
           and Learning-Enabled Cyber-Physical Systems},
  booktitle = {Computer Aided Verification (CAV)},
  year = {2020},
  doi = {10.1007/978-3-030-53288-8_1}
}
```

NNV3 paper is in preparation. For feature-specific citations (FairNNV, VolumeStar, GNNV, ModelStar, probabilistic verification), see the [How to Cite](https://verivital.github.io/nnv/publications/citing.html) page.

For the full list of 30+ publications using NNV, see the [Publications](https://verivital.github.io/nnv/publications/papers.html) page.

## Contributors

* [Hoang-Dung Tran](https://sites.google.com/view/v2a2/)
* [Diego Manzanas Lopez](https://mldiego.github.io/)
* [Anne Tumlin](https://atumlin.github.io/)
* [Samuel Sasaki](https://sammsaski.github.io/)
* [Ben Wooding](https://woodingben.com)
* [Neelanjana Pal](https://scholar.google.com/citations?user=3j_f-ewAAAAJ&hl=en)
* [Weiming Xiang](https://xiangweiming.github.io/)
* [Stanley Bak](http://stanleybak.com/)
* [Patrick Musau](https://pmusau17.github.io/)
* [Xiaodong Yang](https://scholar.google.com/citations?user=xe3Jr7EAAAAJ&hl=en)
* [Luan Viet Nguyen](https://luanvietnguyen.github.io)
* [Taylor T. Johnson](http://www.taylortjohnson.com)

## Acknowledgements

This work is supported in part by AFOSR, DARPA, NSF.

## Contact

For questions, please [open an issue](https://github.com/verivital/nnv/issues) or contact [Diego Manzanas Lopez](mailto:diego.manzanas.lopez@vanderbilt.edu).
