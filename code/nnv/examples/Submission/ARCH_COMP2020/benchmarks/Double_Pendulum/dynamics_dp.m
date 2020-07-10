function [dx] = dynamics_dp(t,x,T)
% Ex_double pendulum 
% constants as per the documentation m=0.5, L=0.5, c= 0, g=1.0
% description of the pararameters in these equations can be found here
% https://github.com/amaleki2/benchmark_closedloop_verification/blob/master/AINNC_benchmark.pdf
% Equations with different parameters of m, L, c,g, can be found in
% double_pendulum_solver



th1 = x(1);
th2 = x(2);
u1 = x(3);
u2 = x(4);
T1 = T(1);
T2 = T(2);

dx(1,1) = x(3);
dx(2,1) = x(4);
dx(3,1) = 4*T1 + 2*sin(th1) - (u2^2*sin(th1 - th2))/2 + (cos(th1 - th2)*(sin(th1 - th2)*u1^2 + 8*T2 + 2*sin(th2) - cos(th1 - th2)*(- (sin(th1 - th2)*u2^2)/2 + 4*T1 + 2*sin(th1))))/(2*(cos(th1 - th2)^2/2 - 1));
dx(4,1) = -(sin(th1 - th2)*u1^2 + 8*T2 + 2*sin(th2) - cos(th1 - th2)*(- (sin(th1 - th2)*u2^2)/2 + 4*T1 + 2*sin(th1)))/(cos(th1 - th2)^2/2 - 1);
end