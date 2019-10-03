function [dx]=dynamics(t,x,u)

dx(1,1)= x(2);
dx(2,1) = -9.8 * x(3) + 1.6 * x(3)^3 + x(1) * x(4)^2;
dx(3,1) = x(4);
dx(4) = u;

