# nnv
Matlab Toolbox For Neural Network Verification

This toolbox implements reachability methods for analyzing neural networks, particularly with a focus on closed-loop controllers in autonomous cyber-physical systems (CPS).

# installation:
    1) install matlab mpt toobox and submodules
       Look at: www.tbxmanager.com.

       command: 1) tbxmanager install mpt mptdoc
                2) tbxmanager install lcp hysdel cddmex clpmex glpkmex fourier sedumi yalmip

    2) clone or download nnv toolbox from https://github.com/trhoangdung/nnv

    3) add nnv toolbox into matlab using setpath in maltab

# running tests and examples

    go into tests/examples folders to run the scripts for testing/analyzing examples.

# contributors

* Hoang-Dung Tran
* Weiming Xiang
* Patrick Musau
* Diego Manzanas Lopez
* Xiaodong Yang
* [Taylor T. Johnson](http://www.taylortjohnson.com)

# references

The methods implemented in this toolbox are based upon the following papers.

1. Weiming Xiang, Hoang-Dung Tran, Taylor T. Johnson, "Output Reachable Set Estimation and Verification for Multi-Layer Neural Networks", In IEEE Transactions on Neural Networks and Learning Systems (TNNLS), 2018, March. [http://taylortjohnson.com/research/xiang2018tnnls.pdf]

2. Weiming Xiang, Hoang-Dung Tran, Taylor T. Johnson, "Reachable Set Computation and Safety Verification for Neural Networks with ReLU Activations", In In Submission, IEEE, 2018, September. [http://www.taylortjohnson.com/research/xiang2018tcyb.pdf]

3. Weiming Xiang, Diego Manzanas Lopez, Patrick Musau, Taylor T. Johnson, "Reachable Set Estimation and Verification for Neural Network Models of Nonlinear Dynamic Systems", In Unmanned System Technologies: Safe, Autonomous and Intelligent Vehicles, Springer, 2018, September. [http://www.taylortjohnson.com/research/xiang2018ust.pdf]

4. Reachability Analysis and Safety Verification for Neural Network Control Systems, Weiming Xiang, Taylor T. Johnson [https://arxiv.org/abs/1805.09944]

5. Specification-Guided Safety Verification for Feedforward Neural Networks, Weiming Xiang, Hoang-Dung Tran, Taylor T. Johnson [https://arxiv.org/abs/1812.06161]
