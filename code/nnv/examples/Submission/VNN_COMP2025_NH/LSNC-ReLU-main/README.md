# VNN Benchmark: LSNC-ReLU

This is a benchmark for [VNN-COMP 2025](https://sites.google.com/view/vnn2025), a simplied version of [LSNC 2024](https://github.com/shizhouxing/LSNC_VNNCOMP2024). We simplify the setting to let more tools participate in this benchmark. Compared to last year, **only ReLU** nonlinearity is presented in the network (no trig functions, leakyrelu, etc.). The input dimension is **low** (dim=6), and only **small shallow networks** are used (275 parameters in total). The vnnlib specification is also **simple and similar to existing benchmarks**.

## Generate Specifications
```bash
python generate_properties.py SEED
``` 

## Background

The NN controller and Lyapunov function are trained using the technique presented in
[CT-BaB](https://arxiv.org/pdf/2411.18235). Compared to a general nonlinear dynamics as in last year, we utilized a neural-network dynamics this year. More specifically, the dynamics of this model is given by a shallow neural network with ReLU activations that are trained to approximate the 2D quadrotor dynamics. The verification goal is to certify the Lyapunov stability of the NN controller in a nonlinear dynamical system within a certain Lyapunov sublevel set ($V < c$). Specifications for the benchmark are randomly generated and consist of random rings of the form $\{x: c_1 < V(x) < c_2\}$. 

## References 

More details of the problem can be found in our paper, and please kindly cite the paper if you find this benchmark useful. 

Zhouxing Shi, Cho-Jui Hsieh, Huan Zhang "[Certified Training with Branch-and-Bound:
A Case Study on Lyapunov-stable Neural Control](https://arxiv.org/pdf/2411.18235)".

Details on the neural network dynamics can be found in the paper

Hongkai Dai, Benoit Landry, Lujie Yang, Marco Pavone, Russ Tedrake "[Lyapunov-stable neural-network control](https://arxiv.org/pdf/2109.14152)".
