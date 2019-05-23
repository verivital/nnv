function [dx]=dynamics(t,x,u)

dx(1,1)=x(2) + 0.5 * x(3)^2;
dx(2,1) = x(3);
dx(3,1) = u;

