function [dx]=discrete_car_dynamics(t,x,u,T)
% Note that t and T is required for reachability analysis
T = [];
dx(1,1)=x(1) + x(2);
dx(2,1)= -0.0025*cos(3*x(1)) + 0.0015 * u + x(2);
end