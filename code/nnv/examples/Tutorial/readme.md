# Neural Network Verification for Medical Imaging Analysis - NNV Tutorial at SPIE 2025

Previous tutorials at

 - DSN 2024
 - 2023 IEEE IAVVC
 - EMSOFT'23 (Embedded Systems Week 2023)


## Getting Started
We recommend using MATLAB Online (**prior license not required!**) to participate in the interactive portion of this tutorial, but one can also use a local installation of MATLAB and NNV to follow along as well.

Licenses available here, you just need to have a Mathworks account and access Matlab online through this link (thanks to Mathworks!): https://www.mathworks.com/licensecenter/classroom/4715800/

###### Note: this license is not required if participants already have a valid license for all the toolboxes listed in [installation instructions](/README.md#installation).

It is important that participants complete this section before arriving at the tutorial, if you want to participate in the interactive portions. The tutorial begins with a lecture that does not require any involvement.

#### MATLAB Online Instructions

Create a copy of NNV into your MathWorks account (personal MATLAB Drive):
[NNV Online](https://drive.mathworks.com/sharing/980b6270-7a14-451c-924c-ae07745fa279)
  - Click on `Open in MATLAB Online` -> `Copy Folder`
  - This will prompt you to log into your account (or register if you don’t have one)
###### Note: copying NNV online may take anywhere from 15 minutes to a couple of hours.

#### Standalone Installation

To run NNV locally, please follow these [instructions](/README.md#installation).


## Interactive Tutorial

During this tutorial, one will learn how to represent convex sets as Star sets, load neural networks into NNV, compute reachability analysis of classification neural networks, compute the certified robustness of neural networks, and more!

_Setup_

Open [MATLAB Online](https://www.mathworks.com/licensecenter/classroom/4715800/) (or MATLAB), then go to the directory where NNV exists on your machine (/code/nnv/), then run the `startup_nnv.m` script ([startup_nnv.m](/code/nnv/startup_nnv.m)).
    
###### Note: if you restart MATLAB, rerun `startup_nnv.m`, which will add the necessary dependencies to the path; you alternatively can run `savepath` after installation to avoid this step after restarting Matlab, but this may require administrative privileges.

## Medical Imaging Tutorial (SPIE)

* Robustness verification of 2D images
    * Provide examples for verifying dataset (Certified robust accuracy) and single samples [SPIE/Classification2D](SPIE/Classification2D)
 
* Robustness verification of 3D patient data
    * Provide examples for verifying dataset (Certified robust accuracy) and single samples [SPIE/Classification3D](SPIE/Classification3D)
 
* Robustness verification of lesion semantic segmentation
    * Provide examples for verifying single 2D slices from MRI data [SPIE/Segmentation](SPIE/Segmentation)


### Additional examples from other domains

#### Neural Networks (NN)

* Robustness verification on the MNIST dataset.
    * Robustness verification example using a model with fully-connected and ReLU layers (MLP). [[verification](NN/MNIST/verify_fc.m)][[training](NN/MNIST/training_fc.m)]
    * Robustness verification example using a model with Convolutional, Pooling, Batch Normalization, ReLU, and fully-connected layers (CNN). [[verification](NN/MNIST/verify.m)][[training](NN/MNIST/training.m)]
    * Certified robustness of a neural network classifier over the test dataset (1000 samples). [[MLP](NN/MNIST/verify_fc_allTest.m)][[CNN](NN/MNIST/verify_allTest.m)]

* Robustness verification on the GTSRB dataset [[examples](NN/GTSRB)]
    * Includes training and verification scripts (robust, and not robust examples).
* Comparison of exact (sound and complete) and approximate (sound and incomplete) methods using Star sets. [[exact vs approx](NN/compareReachability/reach_exact_vs_approx.m)]
* Robustness analysis of a malware classifier (BODMAS Dataset) [[malware examples](NN/malware)]
    * Includes input set constructions depending on the type of adversarial attack. The perturbation can be applied over continuous features, discrete, or on all input features. 

#### Neural Network Control Systems (NNCS)

* Reachability analysis of an inverted pendulum. [[reachability code](NNCS/InvertedPendulum/reach_invP.m)]
* Safety verification example of an Adaptive Cruise Control (ACC) system.
    * [Training](NNCS/ACC/Training%20and%20testing). This requires installing Simulink.
    * [Safety Verification](NNCS/ACC/Verification/verify.m).
* Safety verification of an Automated Emergency Braking System ([AEBS](NNCS/AEBS))
    * This system contains several neural networks, so a custom control loop is included in the verification script.

In addition, we have also prepared a set of exercises for the participants:
* Verification of ACAS Xu neural network with ONNX and VNNLIB files: [ONNX & VNNLIB exercise](NN/ACAS%20Xu/exercise_vnnlib_onnx.m)
* Robustness verification of a neural network: [NN exercise](NN/GTSRB/exercise_verify_robustness.m)
* Safety verification of a linearized plant model of an adaptive cruise control system: [NNCS exercise](NNCS/ACC/Exercise/exercise_reachability_nncs.m)
