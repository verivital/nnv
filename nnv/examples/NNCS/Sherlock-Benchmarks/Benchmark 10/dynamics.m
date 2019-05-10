function [dx]=dynamics(t,x,u)

dx(1,1)= x(4) * cos(x(3));
dx(2,1) = x(4) * sin(x(3));
dx(3,1) = u(2);
dx(4,1) = u(1);

