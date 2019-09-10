% reachability analysis for nonlinear car model for adaptive cruise control system 
% This function is called every step T = 0.1s to compute the reachable set
% of two cars in the system. 


% the nonlinear dynamics of the car is: 
% x_lead' = v_lead
% v_lead' = -0.5 * v_lead - mu * v_lead^2 
% x_ego' = v_ego
% v_ego' = -0.5 * v_ego - mu * v_ego ^2 + u

% mu: is the friction parameter
% u : is the control input from the neural network controller, u is the
% acceleration set point

% Dung Tran : 11/19/2018
% 

car = NonLinearODE(4, 1, @car_dynamics);

lb = [49; 25; 9; 20];
ub = [51; 25.2; 11; 20.2];

B1 = Box(lb, ub);

init_set = B1.toZono; % initial set of state

lb = 0.9;
ub = 1.0;
B1 = Box(lb, ub);

input_set = B1.toZono; % input set

timeStep = 0.01;
tFinal = 0.1;


[R, reachTime] = car.reach_zono(init_set, input_set, timeStep, tFinal);