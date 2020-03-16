% Plot reachable sets for Discrete Linear ACC model
% Dung Tran: 10/22/2019

path_out = [path_results(), filesep, 'ACC', filesep];	

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

plant = LinearODE(A, B, C, D, 0.1, 20); % continuous plant model

plantd = plant.c2d(0.1); % discrete plant model

% the neural network provides a_ego control input to the plant
% a_lead = -5 


%% Controller
load controller_5_20.mat;

weights = network.weights;
bias = network.bias;
n = length(weights);
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

ncs = DLinearNNCS(Controller, plantd); % a discrete linear neural network control system


%% Initial Set of states and reference inputs

% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

ref_input = [30; 1.4];

% initial condition of the plant

% initial position of lead car x_lead
x_lead = [90 92];
% initial condition of v_lead
v_lead = cell(10, 1);
v_lead{1, 1} = [29 30];
v_lead{2, 1} = [28 29];
v_lead{3, 1} = [27 28];
v_lead{4, 1} = [26 27];
v_lead{5, 1} = [25 26];
v_lead{6, 1} = [24 25];
v_lead{7, 1} = [23 24];
v_lead{8, 1} = [22 23];
v_lead{9, 1} = [21 22];
v_lead{10, 1} = [20 21];

% initial condition of x_internal_lead
internal_acc_lead = [0 0];
% initial condition of x_ego
x_ego = [30 31]; 
% initial condition of v_ego
v_ego = [30 30.5];
% initial condition of x_internal_ego
internal_acc_ego = [0 0];

a_lead = -5;
% initial condition for new introduced variable 
x7_0 = [2*a_lead 2*a_lead]; % x7 = -2*x(3) + 2 * a_lead -> x7(0) = -2*x(3)(0) + 2*alead = -2*0 + 2*-2 = -4

N = length(v_lead);
for i=1:N
    lb = [x_lead(1); v_lead{i}(1); internal_acc_lead(1); x_ego(1); v_ego(1); internal_acc_ego(1); x7_0(1)];
    ub = [x_lead(2); v_lead{i}(2); internal_acc_lead(2); x_ego(2); v_ego(2); internal_acc_ego(2); x7_0(2)];
    init_set(i) = Star(lb, ub);
end



%% Compute reachable set
reachPRM.init_set = init_set(1);
reachPRM.ref_input = ref_input;
reachPRM.numSteps = 50;
reachPRM.reachMethod = 'approx-star';

if ~exist('numCores')
    reachPRM.numCores = 4;
else
    reachPRM.numCores = numCores;
end

ncs.reach(reachPRM);

save([path_out, 'linear_ACC_ncs_1.mat'], 'ncs');

reachPRM.init_set = init_set(6);
ncs.reach(reachPRM);

save([path_out, 'linear_ACC_ncs_6.mat'], 'ncs');	

%% Plot output reach sets: actual distance vs. safe distance

% plot reachable set of the distance between two cars d = x1 - x4 vs. ego
% car's velocity

figure;
h1 = subplot(2,1,1);
load([path_out, 'linear_ACC_ncs_1.mat']);	
map_mat = [1 0 0 -1 0 0 0];
map_vec = [];
ncs.plotOutputReachSets('blue', map_mat, map_vec);
hold on;
% plot safe distance betwenen two cars: d_safe = D_default + t_gap * v_ego;
% D_default = 10; t_gap = 1.4 
% d_safe = 10 + 1.4 * x5; 

map_mat = [0 0 0 0 1.4 0 0];
map_vec = [10];
ncs.plotOutputReachSets('red', map_mat, map_vec);
title(h1,'Actual Distance (blue) vs. Safe Distance (red)');
xlabel(h1, 'Control Time Steps');
ylabel(h1, 'Distance');
xticks(h1, [0:5:50])

h2 = subplot(2,1,2);
load([path_out, 'linear_ACC_ncs_6.mat']);	
map_mat = [1 0 0 -1 0 0 0];
map_vec = [];
ncs.plotOutputReachSets('blue', map_mat, map_vec);
hold on;
% plot safe distance between two cars: d_safe = D_default + t_gap * v_ego;
% D_default = 10; t_gap = 1.4 
% d_safe = 10 + 1.4 * x5; 

map_mat = [0 0 0 0 1.4 0 0];
map_vec = [10];
ncs.plotOutputReachSets('red', map_mat, map_vec);
title(h2,'Actual Distance (blue) vs. Safe Distance (red)');
xlabel(h2, 'Control Time Steps');
ylabel(h2, 'Distance');
xticks(h2, [0:5:50])

saveas(gcf, [path_out, 'figure4_plot_linear_ACC_reachSets.png']);

%% END OF SCRIPT