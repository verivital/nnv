function [dx]=dynamics(t,x,u)

dx(1,1)= -x(1) + x(2) - x(3);
dx(2,1) = -x(1) * (x(3) + 1) - x(2);
dx(3,1) = -x(1) + u;

