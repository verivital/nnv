function [dx]=dynamics(t,x,u)

dx(1,1)=x(2);
dx(2,1) = -x(1) + 0.1 * sin(x(3));
dx(3,1) = x(4);
dx(4,1) = u;

