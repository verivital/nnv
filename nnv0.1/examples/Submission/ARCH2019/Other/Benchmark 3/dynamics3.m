function [dx]=dynamics(t,x,u)

dx(1,1)= -x(1) * (0.1 + (x(1) + x(2))^2);
dx(2,1) = (u + x(1)) * (0.1 + (x(1) + x(2))^2);

