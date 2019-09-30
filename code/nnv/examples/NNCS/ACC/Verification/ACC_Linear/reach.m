% Reachability analysis for Linear ACC model
% Dung Tran: 9/30/2019



%% System model
% x1 = lead_car position
% x2 = lead_car velocity
% x3 = lead_car internal state

% x4 = ego_car position
% x5 = ego_car velocity
% x6 = ego_car internal state

% lead_car dynamics
%dx(1,1)=x(2);
%dx(2,1) = x(3);
%dx(3,1) = -2 * x(3) + 2 * a_lead - mu*x(2)^2;

% ego car dynamics
%dx(4,1)= x(5); 
%dx(5,1) = x(6);
%dx(6,1) = -2 * x(6) + 2 * a_ego - mu*x(5)^2;


A = [0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 -2 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 -2];
B = transpose([0 0 2 0 0 0; 0 0 0 0 0 2]);
C = [1 0 0 -1 0 0; 0 1 0 0 -1 0; 0 0 0 0 1 0]; % feedback relative distance, relative velocity, longitudinal velocity

D = [0 0; 0 0; 0 0]; 

plant = LinearODE(A, B, C, D);


%% Controller
load controller_3_20.mat;

n = length(weights);
Layers = [];
for i=1:n - 1
    L = LayerS(weights{1, i}, bias{i, 1}, 'poslin');
    Layers = [Layers L];
end
L = LayerS(weights{1, n}, bias{n, 1}, 'purelin');
Layers = [Layers L];
Controller = FFNNS(Layers); % feedforward neural network controller


%% NNCS

feedbackMap = [0]; % feedback map, y[k] 
ncs = NNCS(Controller, plant, feedbackMap); % the neural network control system
