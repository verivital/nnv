function [dx] = dynamics_sp(t,x,a)
% Ex_single_pendulum
% constants as per the documentation l=0.5, m=0.5; g= 1; c=0;
% description of the pararameters in these equations can be found here
% https://github.com/amaleki2/benchmark_closedloop_verification/blob/master/AINNC_benchmark.pdf

l=0.5; m=0.5; g= 1; c= 0;

dx(1,1) = x(2);
dx(2,1) = g/l * sin(x(1)) + (a - c*x(2))/(m*l^2);
dx(3,1) = 20;