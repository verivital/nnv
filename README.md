# NNV
Matlab Toolbox for Neural Network Verification

This toolbox implements reachability methods for analyzing neural networks and control systems with neural network controllers in the area of autonomous cyber-physical systems (CPS).

## Related tools and software

This toolbox makes use of the neural network model transformation tool ([nnmt](https://github.com/verivital/nnmt)) and for closed-loop systems analysis, the hybrid systems model transformation and translation tool ([HyST](https://github.com/verivital/hyst)), and the COntinuous Reachability Analyzer ([CORA](https://github.com/TUMcps/CORA)).

## Execution without installation:
NNV can be executed online without installing Matlab or other dependencies through [CodeOcean](https://www.codeocean.com) via the following CodeOcean capsules:

_Latest_
* CAV 2023 Tool Paper: https://doi.org/10.24433/CO.0803700.v1
  
_Previous_
* CAV 2020 ImageStar paper version: https://doi.org/10.24433/CO.3351375.v1
* CAV 2020 Tool paper version: https://doi.org/10.24433/CO.0221760.v1
* Earliest version: https://doi.org/10.24433/CO.1314285.v1

# Installation:

1. Clone or download the NNV toolbox from (https://github.com/verivital/nnv)
    
    Note: to operate correctly, nnv depends on other tools (CORA, NNMT, HyST, onnx2nnv), which are included as git submodules. As such, you must clone recursively, e.g., with the following:
    ```
    git clone --recursive https://github.com/verivital/nnv.git
    ```

2. If running in Ubuntu, install MATLAB and proceed to run the provided installation script (then, skip to step 6). 

    ```
    chmod +x install_ubuntu.sh
    ./install_ubuntu.sh
    ``` 
    
3. For MacOS and Windows, please install MATLAB (2023a or newer) with at least the following toolboxes:
   * Computer Vision
   * Control Systems
   * Deep Learning
   * Image Processing
   * Optimization
   * Parallel Computing
   * Statistics and Machine Learning
   * Symbolic Math
   * System Identification
   
4. Install the following support package
       [Deep Learning Toolbox Converter for ONNX Model Format](https://www.mathworks.com/matlabcentral/fileexchange/67296-deep-learning-toolbox-converter-for-onnx-model-format)
       
       Note: Support packages can be installed in MATLAB's HOME tab > Add-Ons > Get Add-ons, search for the support package using the Add-on Explorer and click on the Install button.


5. Open MATLAB, then go to the directory where NNV exists on your machine, then run the `install.m` script located at /nnv/
    
    Note: if you restart Matlab, rerun either install.m or startup_nnv.m, which will add the necessary dependencies to the path; you alternatively can run `savepath` after installation to avoid this step after restarting Matlab, but this may require administrative privileges

6. Optional installation packages

   a. To run verification for convolutional neural networks (CNNs) on VGG16/VGG19, additional support packages must be installed:

   * [VGG16](https://www.mathworks.com/matlabcentral/fileexchange/61733-deep-learning-toolbox-model-for-vgg-16-network)

   * [VGG19](https://www.mathworks.com/help/deeplearning/ref/vgg19.html)
       
    b) To run MATLAB's neural network verification comparison, an additional support package is needed (used in CAV'2023 submission):
        
    * [Deep Learning Toolbox Verification Library](https://www.mathworks.com/matlabcentral/fileexchange/118735-deep-learning-toolbox-verification-library)

    c) To load models from other deep learning frameworks, please install the additional support packages:

    * TensorFlow and Keras: [Deep Learning Toolbox Converter for TensorFlow Models](https://www.mathworks.com/matlabcentral/fileexchange/64649-deep-learning-toolbox-converter-for-tensorflow-models)
    * PyTorch: [Deep Learning Toolbox Converter for PyTorch Models](https://www.mathworks.com/matlabcentral/fileexchange/111925-deep-learning-toolbox-converter-for-pytorch-models)
        
## Uninstallation:

Open MATLAB, then go to the `/code/nnv/` folder and execute the `uninstall.m` script.

# Getting started with NNV 

### [Tutorial](code/nnv/examples/Tutorial)

To get started with NNV, let's take a look at a tutorial containing examples demonstrating:

__NN__

* Robustness verification on the MNIST dataset.
    * Includes model training and several verification examples.
* Robustness verification on the GTSRB dataset.
    * Includes model training and robustness verification.
* Comparisons of exact (sound and complete) and approximate (sound and incomplete) methods using Star sets
    * Visualize the size difference on the output sets and the computation times for each method.  
* Robustness analysis of a malware classifier (BODMAS Dataset).

  
__NNCS__
  
* Reachability analysis of an inverted pendulum.
* Safety verification example of an Adaptive Cruise Control (ACC) system.
* Safety verification of an Automated Emergency Braking System

And more! Please go to the [tutorial description](code/nnv/examples/Tutorial/readme.md) for more details!

### Examples

In addition to the examples from the tutorial, there are more examples in the 'code/nnv/examples/' folder, including:

**Semantic Segmentation**

* [Robustness analysis of semantic segmentation NNs](code/nnv/examples/NN/SemanticSegmentation/M2NIST)

**Recurrent Neural Networks**

* [Robustness analysis of RNNs](code/nnv/examples/NN/RNN)

**Neural Ordinary Differential Equations**

* [Reachability analysis of neural ODEs](code/nnv/examples/NN/NeuralODEs)

And more other [NN](code/nnv/examples/NN) and [NNCS](code/nnv/examples/NNCS) examples.

### Tests

To run all the tests, one can run the following command from 'code/nnv/tests/' folder:

`runtests(pwd, 'IncludeSubfolders', true);`

### Paper Publications and Competitions

All the code for the publication using NNV, including competitions that NNV has participated in (e.g. [VNNCOMP](https://sites.google.com/view/vnn2023) and [ARCH-COMP](https://cps-vo.org/group/ARCH/FriendlyCompetition)) can be found in the [examples/Submissions](code/nnv/examples/Submission) folder. For a selected subset of publications, [tags](https://github.com/verivital/nnv/tags) were created. To reproduce those results, please download NNV using the corresponding tag or use the corresponding CodeOcean capsules.


## Contributors

* [Hoang-Dung Tran](https://sites.google.com/view/v2a2/)
* [Diego Manzanas Lopez](https://mldiego.github.io/)
* [Neelanjana Pal](https://scholar.google.com/citations?user=3j_f-ewAAAAJ&hl=en)
* [Weiming Xiang](https://xiangweiming.github.io/)
* [Stanley Bak](http://stanleybak.com/)
* [Patrick Musau](https://pmusau17.github.io/)
* [Xiaodong Yang](https://scholar.google.com/citations?user=xe3Jr7EAAAAJ&hl=en) 
* [Luan Viet Nguyen](https://luanvietnguyen.github.io)
* [Taylor T. Johnson](http://www.taylortjohnson.com)

## References

The methods implemented in NNV are based upon or used in the following papers:

* Diego Manzanas Lopez, Sung Woo Choi, Hoang-Dung Tran, Taylor T. Johnson, "NNV 2.0: The Neural Network Verification Tool". In: Enea, C., Lal, A. (eds) Computer Aided Verification. CAV 2023. Lecture Notes in Computer Science, vol 13965. Springer, Cham. [https://doi.org/10.1007/978-3-031-37703-7_19]

* Anne M Tumlin, Diego Manzanas Lopez, Preston Robinette, Yuying Zhao, Tyler Derr, and Taylor T Johnson. "FairNNV: The Neural Network Verification Tool For Certifying Fairness. In Proceedings of the 5th ACM International Conference on AI in Finance (ICAIF '24)". Association for Computing Machinery, New York, NY, USA, 36–44. [https://doi.org/10.1145/3677052.3698677]

* Taylor T. Johnson, Diego Manzanas Lopez and Hoang-Dung. Tran, "Tutorial: Safe, Secure, and Trustworthy Artificial Intelligence (AI) via Formal Verification of Neural Networks and Autonomous Cyber-Physical Systems (CPS) with NNV," 2024 54th Annual IEEE/IFIP International Conference on Dependable Systems and Networks - Supplemental Volume (DSN-S), Brisbane, Australia, 2024, pp. 65-66, [https://doi.org/10.1109/DSN-S60304.2024.00027]

* Preston K. Robinette, Diego Manzanas Lopez, Serena Serbinowska, Kevin Leach, and Taylor T Johnson. "Case Study: Neural Network Malware Detection Verification for Feature and Image Datasets". In Proceedings of the 2024 IEEE/ACM 12th International Conference on Formal Methods in Software Engineering (FormaliSE) (FormaliSE '24). Association for Computing Machinery, New York, NY, USA, 127–137. [https://doi.org/10.1145/3644033.3644372]

* Hoang-Dung Tran, Diego Manzanas Lopez, and Taylor Johnson. "Tutorial: Neural Network and Autonomous Cyber-Physical Systems Formal Verification for Trustworthy AI and Safe Autonomy". In Proceedings of the International Conference on Embedded Software (EMSOFT '23). Association for Computing Machinery, New York, NY, USA, 1–2. [https://doi.org/10.1145/3607890.3608454]

* Neelanjana Pal, Diego Manzanas Lopez, Taylor T Johnson, "Robustness Verification of Deep Neural Networks using Star-Based Reachability Analysis with Variable-Length Time Series Input", to be presented at FMICS 2023. [https://arxiv.org/pdf/2307.13907.pdf]

* Mykhailo Ivashchenko, Sung Woo Choi, Luan Viet Nguyen, Hoang-Dung Tran, "Verifying Binary Neural Networks on Continuous Input Space using Star Reachability," 2023 IEEE/ACM 11th International Conference on Formal Methods in Software Engineering (FormaliSE), Melbourne, Australia, 2023, pp. 7-17, [https://doi.org/10.1109/FormaliSE58978.2023.00009]

* Hoang Dung Tran, Sung Woo Choi, Xiaodong Yang, Tomoya Yamaguchi, Bardh Hoxha, and Danil Prokhorov. "Verification of Recurrent Neural Networks with Star Reachability". In Proceedings of the 26th ACM International Conference on Hybrid Systems: Computation and Control (HSCC '23). Association for Computing Machinery, New York, NY, USA, Article 6, 1–13. [https://doi.org/10.1145/3575870.3587128]

* Diego Manzanas Lopez, Taylor T. Johnson, Stanley Bak, Hoang-Dung Tran, Kerianne Hobbs, "Evaluation of Neural Network Verification Methods for Air to Air Collision Avoidance", In AIAA Journal of Air Transportation (JAT), 2022 [http://www.taylortjohnson.com/research/lopez2022jat.pdf]

* Diego Manzanas Lopez, Patrick Musau, Nathaniel Hamilton, Taylor T. Johnson, "Reachability Analysis of a General Class of Neural Ordinary Differential Equations", In 20th International Conference on Formal Modeling and Analysis of Timed Systems (FORMATS), 2022 [http://www.taylortjohnson.com/research/lopez2022formats.pdf]

* Hoang-Dung Tran, Neelanjana Pal, Patrick Musau, Xiaodong Yang, Nathaniel P. Hamilton, Diego Manzanas Lopez, Stanley Bak, Taylor T. Johnson, "Robustness Verification of Semantic Segmentation Neural Networks using Relaxed Reachability", In 33rd International Conference on Computer-Aided Verification (CAV), Springer, 2021. [http://www.taylortjohnson.com/research/tran2021cav.pdf]

* Hoang-Dung Tran, Patrick Musau, Diego Manzanas Lopez, Xiaodong Yang, Luan Viet Nguyen, Weiming Xiang, Taylor T.Johnson, "NNV: A Tool for Verification of Deep Neural Networks and Learning-Enabled Autonomous Cyber-Physical Systems", 32nd International Conference on Computer-Aided Verification (CAV), 2020. [http://taylortjohnson.com/research/tran2020cav_tool.pdf]

* Hoang-Dung Tran, Stanley Bak, Weiming Xiang, Taylor T. Johnson, "Towards Verification of Large Convolutional Neural Networks Using ImageStars", 32nd International Conference on Computer-Aided Verification (CAV), 2020. [http://taylortjohnson.com/research/tran2020cav.pdf]

* Stanley Bak, Hoang-Dung Tran, Kerianne Hobbs, Taylor T. Johnson, "Improved Geometric Path Enumeration for Verifying ReLU Neural Networks", In 32nd International Conference on Computer-Aided Verification (CAV), 2020. [http://www.taylortjohnson.com/research/bak2020cav.pdf]

* Hoang-Dung Tran, Weiming Xiang, Taylor T. Johnson, "Verification Approaches for Learning-Enabled Autonomous Cyber-Physical Systems", The IEEE Design & Test 2020. [http://www.taylortjohnson.com/research/tran2020dandt.pdf]

* Hoang-Dung Tran, Patrick Musau, Diego Manzanas Lopez, Xiaodong Yang, Luan Viet Nguyen, Weiming Xiang, Taylor T.Johnson, "Star-Based Reachability Analysis for Deep Neural Networks", The 23rd International Symposium on Formal Methods (FM), Porto, Portugal, 2019, Acceptance Rate 30%. . [http://taylortjohnson.com/research/tran2019fm.pdf]

* Hoang-Dung Tran, Feiyang Cei, Diego Manzanas Lopez, Taylor T.Johnson, Xenofon Koutsoukos, "Safety Verification of Cyber-Physical Systems with Reinforcement Learning Control",  The International Conference on Embedded Software (EMSOFT), New York, October 2019. Acceptance Rate 25%. [http://taylortjohnson.com/research/tran2019emsoft.pdf]

* Hoang-Dung Tran, Patrick Musau, Diego Manzanas Lopez, Xiaodong Yang, Luan Viet Nguyen, Weiming Xiang, Taylor T.Johnson, "Parallelizable Reachability Analysis Algorithms for FeedForward Neural Networks", In 7th International Conference on Formal Methods in Software Engineering (FormaLISE), 27, May 2019 in Montreal, Canada, Acceptance Rate 28%. [http://taylortjohnson.com/research/tran2019formalise.pdf]

* Diego Manzanas Lopez, Patrick Musau, Hoang-Dung Tran, Taylor T.Johnson, "Verification of Closed-loop Systems with Neural Network Controllers (Benchmark Proposal)", The 6th International Workshop on Applied Verification of Continuous and Hybrid Systems (ARCH2019). Montreal, Canada, 2019. [http://taylortjohnson.com/research/lopez2019arch.pdf]

* Weiming Xiang, Hoang-Dung Tran, Taylor T. Johnson, "Output Reachable Set Estimation and Verification for Multi-Layer Neural Networks", In IEEE Transactions on Neural Networks and Learning Systems (TNNLS), 2018, March. [http://taylortjohnson.com/research/xiang2018tnnls.pdf]

* Weiming Xiang, Hoang-Dung Tran, Taylor T. Johnson, "Reachable Set Computation and Safety Verification for Neural Networks with ReLU Activations", In In Submission, IEEE, 2018, September. [http://www.taylortjohnson.com/research/xiang2018tcyb.pdf]

* Weiming Xiang, Diego Manzanas Lopez, Patrick Musau, Taylor T. Johnson, "Reachable Set Estimation and Verification for Neural Network Models of Nonlinear Dynamic Systems", In Unmanned System Technologies: Safe, Autonomous and Intelligent Vehicles, Springer, 2018, September. [http://www.taylortjohnson.com/research/xiang2018ust.pdf]

* Reachability Analysis and Safety Verification for Neural Network Control Systems, Weiming Xiang, Taylor T. Johnson [https://arxiv.org/abs/1805.09944]

* Weiming Xiang, Patrick Musau, Ayana A. Wild, Diego Manzanas Lopez, Nathaniel Hamilton, Xiaodong Yang, Joel Rosenfeld, Taylor T. Johnson, "Verification for Machine Learning, Autonomy, and Neural Networks Survey," October 2018, [https://arxiv.org/abs/1810.01989]

* Specification-Guided Safety Verification for Feedforward Neural Networks, Weiming Xiang, Hoang-Dung Tran, Taylor T. Johnson [https://arxiv.org/abs/1812.06161]

#### Cite

```
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
booktitle = {Computer Aided Verification: 35th International Conference, CAV 2023, Paris, France, July 17–22, 2023, Proceedings, Part II},
pages = {397–412},
numpages = {16},
keywords = {neural networks, cyber-physical systems, verification, tool},
location = {Paris, France}
}
```

```
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
booktitle = {Computer Aided Verification: 32nd International Conference, CAV 2020, Los Angeles, CA, USA, July 21–24, 2020, Proceedings, Part I},
pages = {3–17},
numpages = {15},
keywords = {Autonomy, Verification, Cyber-physical systems, Machine learning, Neural networks},
location = {Los Angeles, CA, USA}
}
```

### Acknowledgements

This work is supported in part by AFOSR, DARPA, NSF.

### Contact

For any questions related to NNV, please add them to the issues or contact [Diego Manzanas Lopez](mailto:diego.manzanas.lopez@vanderbilt.edu) or [Hoang Dung Tran](mailto:trhoangdung@gmail.com).

