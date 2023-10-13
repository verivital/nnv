# Verification and Validation of Neural Networks in Automated Vehicles (NNV Tutorial) - E2023 IEEE IAVVC

Previous tutorial at EMSOFT'23 (Embedded Systems Week 2023)

## Getting Started
We recommend using MATLAB Online (**prior license not required!**) to participate in this tutorial, but one can also use a local installation of MATLAB and NNV to follow along as well.

It is important that participants complete this section before arriving at the tutorial.

#### MATLAB Online Instructions

Registered participants should have received a link with a MATLAB license to use in this tutorial. If a participant has not received it, 
please contact the organizers. 

###### Note: this license is not required if participants already have a valid license for all the toolboxes listed in [installation instructions](/README.md#installation).

Create a copy of NNV into your MathWorks account (personal MATLAB Drive):
[NNV Online](https://drive.matlab.com/sharing/6861a9c9-f818-45ab-b7cc-380bbfb2be0e)
  - Click on `Open in MATLAB Online` -> `Copy Folder`
  - This will prompt you to log into your account (or register if you don’t have one)
###### Note: copying NNV online may take anywhere from 15 minutes to a couple of hours.

#### Standalone Installation

To run NNV locally, please follow these [instructions](/README.md#installation).


## Interactive Tutorial

During this tutorial, one will learn how to represent convex sets as Star sets, load neural networks into NNV, 
compute reachability analysis of classification neural networks, compute the certified robustness of neural networks and more!

_Setup_

Open [MATLAB Online](https://workshop-matlab.mathworks.com/) (or MATLAB), then go to the directory where NNV exists on your machine (/code/nnv/), then run the `startup_nnv.m` script ([startup_nnv.m](/code/nnv/startup_nnv.m)).
    
###### Note: if you restart MATLAB, rerun `startup_nnv.m`, which will add the necessary dependencies to the path; you alternatively can run `savepath` after installation to avoid this step after restarting Matlab, but this may require administrative privileges.


#### Neural Networks (NN)

* Robustness verification on the MNIST dataset.
    * [Robustness verification of a single image](code/nnv/examples/Tutorial/NN/MNIST/verify_fc.m) using a [model with fully-connected and ReLU layers](code/nnv/examples/Tutorial/NN/MNIST/training_fc.m).
    * [Robustness verification of a single image](code/nnv/examples/Tutorial/NN/MNIST/verify.m) using a [model with Convolutional, Pooling, Batch Normalization, ReLU, and fully-connected layers](code/nnv/examples/Tutorial/NN/MNIST/training.m)
    * [Certified robustness](code/nnv/examples/Tutorial/NN/MNIST/verify_fc_allTest.m) of a neural network classifier using a [model with fully-connected and ReLU layers](code/nnv/examples/Tutorial/NN/MNIST/training_fc.m).

* Robustness verification on the GTSRB dataset.
    * Includes [training](code/nnv/examples/Tutorial/NN/GTSRB/train.m) and [verification](code/nnv/examples/Tutorial/NN/GTSRB/verify_robust_27.m) scripts as well.
* Comparison of exact (sound and complete) and approximate (sound and incomplete) methods using Star sets [exact vs approx](code/nnv/examples/Tutorial/NN/compareReachability/reach_exact_vs_approx.m)

#### Neural Network Control Systems (NNCS)

* Reachability analysis of an [inverted pendulum](code/nnv/examples/Tutorial/NNCS/InvertedPendulum/reach_invP.m).
* Safety verification example of an Adaptive Cruise Control (ACC) system.
    * [Training](code/nnv/examples/Tutorial/NNCS/ACC/Training%20and%20testing). This requires installing Simulink.
    * [Safety Verification](code/nnv/examples/Tutorial/NNCS/ACC/Verification/verify.m).
* Safety verification of an Automated Emergency Braking System ([AEBS](code/nnv/examples/Tutorial/NNCS/AEBS))
    * This system contains several neural networks, so a custom control loop is included in the verification script.

In addition, we have also prepared a set of exercises for the participants:
* Verification of ACAS Xu neural network with ONNX and VNNLIB files: [ONNX & VNNLIB exercise](code/nnv/examples/Tutorial/NN/ACAS%20Xu/exercise_vnnlib_onnx.m)
* Robustness verification of a neural network: [NN exercise](code/nnv/examples/Tutorial/NN/GTSRB/exercise_verify_robustness.m))
* Safety verification of a linearized plant model of an adaptive cruise control system: [NNCS exercise](code/nnv/examples/Tutorial/NNCS/ACC/Exercise/exercise_reachability_nncs.m)
