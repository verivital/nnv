function [dx]=dynamics(t,x,u)

dx(1,1)= -x(1)^3 + x(2);
dx(2,1) = x(2)^3 + x(3);
dx(3,1) = u;

