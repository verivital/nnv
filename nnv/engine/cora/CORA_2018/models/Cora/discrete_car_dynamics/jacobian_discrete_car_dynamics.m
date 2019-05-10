function [A,B]=jacobian_discrete_car_dynamics(x,u,T)

A=[1,1/10;...
(3*sin(3*x(1)))/400,1];

B=[0;...
3/2000];

