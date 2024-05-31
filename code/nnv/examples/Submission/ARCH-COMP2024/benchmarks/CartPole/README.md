# Proposed ARCH Benchmark - CartPole

This benchmark is proposed for the ARCH Friendly Competition 2024.

## Benchmark

We consider a pendulum (pole) mounted on a movable cart (= CartPole). The cart can be moved by a controller. The carts postition x1, its velocity x2, the angle of the pole x3 (with 0, 2*pi being the upright postion) and its angular velocity x4 define the state vector of the 4-dimensional system. The system starts in a middle postition of the cart, with the pendulum in the upright position. The controllers goal is to counteract slight deviations in the starting values and balance the pendulum in the middle of the track.

## Specifications and Dynamics

The system dynamics can be found in ```dynamics.m```. The safe states are defined as a stable upright position, which has to be reached after 8 seconds and has to be hold for at least 2 seconds. The controllers step size is 0.02 seconds.
