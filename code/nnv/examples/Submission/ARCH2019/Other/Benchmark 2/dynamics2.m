function [dx]=dynamics(t,x,u)

dx(1,1)=x(2);
dx(2,1) = u*x(2)^2 - x(1);

