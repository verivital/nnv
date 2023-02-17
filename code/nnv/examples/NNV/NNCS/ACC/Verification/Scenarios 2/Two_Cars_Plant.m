% Dung Tran: 11/13/2018
% Scenarios: lead car acceleration = 0
% expectation: the ego car reduces the velocity to keep a safe distance
% between two cars.

% combine with position: dot{x} = v; 
% the model of two cars is 

% model of lead car
% x_lead' = v_lead; this is lead car
% v_lead' = -0.5*v_lead + a_lead = -0.5*v_lead 
 

% model of ego car
% x_ego' = v_ego; this is ego car
% v_ego' = -0.5*v_ego + a_ego, a_ego is come from the controller

% distance between two cars
% dis = x_lead - x_ego 

% safe distance between two cars, see here 
% https://www.mathworks.com/help/mpc/examples/design-an-adaptive-cruise-control-system-using-model-predictive-control.html

% dis_safe = D_default + t_gap * v_ego;
% the safety specification is: dis >= dis_safe

% the final state space model of two cars has 4 state variables

% x = [x_lead v_lead x_ego v_ego]'


A = [0 1 0 0; 0 -0.5 0 0; 0 0 0 1; 0 0 0 -0.5];
B = [0; 0; 0; 1];
C = [1 0 -1 0; 0 1 0 -1; 0 0 0 1]; % feedback relative distance, relative velocity, longitudinal velocity

D = [0; 0; 0]; 

sys = LinearODE(A, B, C, D);

Ts = 0.1;

sysd = sys.c2d(Ts);







