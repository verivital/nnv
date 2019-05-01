
# Closed-Loop Systems with Neural Network Controllers
Contains all the benchmarks presented at ARCH 2019. These benchmarks consists on a continuous-time dynamical model controlled
by a feedforward neural network with Rectified Linear Units (ReLU) and linear activation functions.

There are six main benchmarks (ACC, Acrobot, CartPole, Inverted Pendulum, MPC Quadrotor, and Mountain Car) introduced here in 
addition to the benchmarks previously presented by Souradeep Duta and tested using Sherlock in <i>
Reachability Analysis for Neural Feedback Systems using Regressive Polynomial Rule Inference, </i> by Souradeep Dutta, Xin Chen and Sriram Sankaranarayanan HSCC 2019, Montreal, Canada 
https://github.com/souradeep-111/Neural-Network-Controller-Verification-Benchmarks-HSCC-2019. 

### Requirements

MATLAB R2018b 

* Simulink

* Stateflow

* Deep Learning Toolbox

### Simulation files

In the main folder for each benchmark we find a simulink model "xxxxxx".slx that contains the neural network 
controller (Simulink block), and the system's plant (Stateflow block). This is the file we will use to simulate the system by simply
clicking on the <b> PLAY </b> button. All the simulations has an initial state point chosen from the initial state sets described in the paper.

### Plant files

In the main folder for each benchmark there is a dynamics.m file that contains the differential equations that define this particular system, as well as the SpaceEx (cfg and xml files) formats for the plants. 

### Controller files

There are two types of files that contain the neural network information, "controller".mat and "controller".onnx. The first one is the 
format used by NNV https://github.com/verivital/nnv, and the latter one correspond to the Open Neural Network Exchange format (ONNX) 
https://github.com/onnx/onnx. 


#### Author

1. Diego Manzanas Lopez
2. Patrick Musau

#### Contact

diego.manzanas.lopez@vanderbilt.edu

patrick.musau@vanderbilt.edu
