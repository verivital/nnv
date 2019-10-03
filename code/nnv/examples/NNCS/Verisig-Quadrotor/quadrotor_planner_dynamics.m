function [dx]=quadrotor_planner_dynamics(t,x,u)

dx(1,1)=x(2);
dx(2,1)= -0.0025*cos(3*x(1)) + 0.0015 * u;


