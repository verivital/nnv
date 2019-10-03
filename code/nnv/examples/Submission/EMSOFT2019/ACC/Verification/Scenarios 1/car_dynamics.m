function [dx]=car_dynamics(t,x,a_ego)

mu=0.0001; % friction parameter

% x1 = lead_car position
% x2 = lead_car velocity
% x3 = lead_car internal state

% x4 = ego_car position
% x5 = ego_car velocity
% x6 = ego_car internal state

% lead car dynamics
a_lead = -2; 
dx(1,1)=x(2);
dx(2,1) = x(3);
dx(3,1) = -2 * x(3) + 2 * a_lead - mu*x(2)^2;
% ego car dyanmics
dx(4,1)= x(5); 
dx(5,1) = x(6);
dx(6,1) = -2 * x(6) + 2 * a_ego - mu*x(5)^2;
