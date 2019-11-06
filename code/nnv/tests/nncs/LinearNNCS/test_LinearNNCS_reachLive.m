% Reachability analysis for Linear ACC model
% Dung Tran: 10/4/2019




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


% let x7 = -2*x(3) + 2 * a_lead -> x7(0) = -2*x(3)(0) + 2*alead
% -> dx7 = -2dx3


A = [0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 0 0 0 1; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 -2 0; 0 0 -2 0 0 0 0];
B = [0; 0; 0; 0; 0; 2; 0];
C = [1 0 0 -1 0 0 0; 0 1 0 0 -1 0 0; 0 0 0 0 1 0 0];  % feedback relative distance, relative velocity, longitudinal velocity
D = [0; 0; 0]; 

plant = LinearODE(A, B, C, D); % continuous plant model

% the neural network provides a_ego control input to the plant
% a_lead = -2 


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

ncs = LinearNNCS(Controller, plant); % a discrete linear neural network control system

%% Initial Set of states and reference inputs

% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

ref_input = [30; 1.4];

% initial condition of the plant

% initial position of lead car x_lead
x_lead = [90 100];
% initial condition of v_lead
v_lead = [32 35];
% initial condition of x_internal_lead
internal_acc_lead = [0 0];
% initial condition of x_ego
x_ego = [10 11]; 
% initial condition of v_ego
v_ego = [30 30.2];
% initial condition of x_internal_ego
internal_acc_ego = [0 0];
% initial condition for new introduced variable 
x7_0 = [-4 -4]; % x7 = -2*x(3) + 2 * a_lead -> x7(0) = -2*x(3)(0) + 2*alead = -2*0 + 2*-2 = -4


x1 = x_lead;
lb = [x1(1); v_lead(1); internal_acc_lead(1); x_ego(1); v_ego(1); internal_acc_ego(1); x7_0(1)];
ub = [x1(2); v_lead(2); internal_acc_lead(2); x_ego(2); v_ego(2); internal_acc_ego(2); x7_0(2)];
init_set = Star(lb, ub);


%% Live Reachability Analysis 1

numSteps = 10; 
method = 'approx-star';
%method = 'exact-star';
numCores = 1;
plantNumOfSimSteps = 20;
% plot on-the-fly the distance between two cars x1-x4 and the safe distance d_safe = D_default + t_gap * v_ego = 10 + 1.4 * x(5)
output_mat = [1 0 0 -1 0 0 0]; % plot on-the-fly the distance between two cars x1 - x4 
output_vec = [];
boundary_mat = [0 0 0 0 1.4 0 0]; % plot on-the-fly the safe distance 
boundary_vec = [10];
figure;
ncs.reachLive(init_set, ref_input, numSteps, method, numCores, output_mat, output_vec, 'blue', boundary_mat, boundary_vec, plantNumOfSimSteps);

%% Live Reachability Analysis 2

numSteps = 10; 
method = 'approx-star';
%method = 'exact-star';
numCores = 1;
% plot on-the-fly the distance between two cars x1-x4 and the safe distance d_safe = D_default + t_gap * v_ego = 10 + 1.4 * x(5)
output_mat = [1 0 0 -1 0 0 0; 0 0 0 0 1 0 0]; % plot on-the-fly the distance between two cars x1 - x4  versus ego car speed
output_vec = [];
boundary_mat = [0 0 0 0 1.4 0 0; 0 0 0 0 1 0 0]; % plot on-the-fly the safe distance vs. ego car speed
boundary_vec = [10;0];
figure;
ncs.reachLive(init_set, ref_input, numSteps, method, numCores, output_mat, output_vec, 'blue', boundary_mat, boundary_vec);

%% END