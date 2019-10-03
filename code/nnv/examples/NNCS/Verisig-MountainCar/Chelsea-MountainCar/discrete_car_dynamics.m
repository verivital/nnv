function [dx]=discrete_car_dynamics(t,x,u, T)
T = [];
% always needs T like this? it is weird to run CORA
dx(1,1)=x(2)*0.1 + x(1);
dx(2,1)= (-0.025*cos(3*x(1)) + 0.015 * u)*0.1 + x(2);


