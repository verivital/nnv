
# NAV Benchmark

## Property:
The control goal is to navigate a robot to a goal region while avoiding an obstacle.
Time horizon: `t = 6s`. Control period: `0.2s`.

Initial states:

    x1 = [2.9, 3.1]
    x2 = [2.9, 3.1]
    x3 = [0, 0]
    x4 = [0, 0]

Dynamic system: [dynamics.m](./dynamics.m)

Goal region ( t=6 ):

    x1 = [-0.5, 0.5]
    x2 = [-0.5, 0.5]
    x3 = [-Inf, Inf]
    x4 = [-Inf, Inf]

Obstacle ( always ):

    x1 = [1, 2]
    x2 = [1, 2]
    x3 = [-Inf, Inf]
    x4 = [-Inf, Inf]

## Networks:

We provide two networks:
- The first network is trained with standard (point-based) reinforcement learning: `nn-nav-point.onnx`
- The second network is trained set-based to improve its verifiable robustness by integrating reachability analysis into the training process: `nn-nav-set.onnx`

Reference set-based training: https://arxiv.org/abs/2401.14961


