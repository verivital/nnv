function [dx]=car_dynamics(t,x, u)

mu=0.1; % friction parameter

dx(1,1)=x(2);
dx(2,1)= -0.5 * x(2)- mu * x(2)^2;
dx(3,1) = x(4);
dx(4,1) = -0.5 * x(4) - mu * x(4)^2  + u;
