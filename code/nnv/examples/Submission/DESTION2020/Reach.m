%% Script to generate the graphs for DESTION 2020
% We are going to evaluate the adaptive cruise control system with
% different controllers and different scenarios. We will also evaluate the
% discrete ACC and the nonlinear dynamics.
clc;clear;close all;

%% Load all components
% Load controller (NN)
load controller_3_20.mat;
weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = [];
for i=1:n - 1
    L = Layer(weights{1, i}, bias{i, 1}, 'ReLU');
    Layers = [Layers L];
end
L = Layer(weights{1, n}, bias{n, 1}, 'Linear');
Layers = [Layers L];

% Load plant (will use CORA for nonlinear dynamics)
Controller = FFNN(Layers); % feedforward neural network controller
Plant = NonLinearODE(6, 1, @car_dynamics);
Plant.set_timeStep(0.01); % time step for reachability analysis of the plant
Plant.set_tFinal(0.1); % Ts = 0.1, sampling time for control signal from neural network controller
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

feedbackMap = [0]; % feedback map, y[k] 
% Create neural network control object
ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

%% Initialize component and set conditions and variables
% initial condition of the Plant

% x = [x_lead v_lead x_internal_lead x_ego v_ego x_internal_ego]'

lb = [94; 32; 0; 10; 30; 0];
ub = [95; 32.2; 0; 11; 30.2; 0];

init_set = Star(lb, ub);

% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

lb1 = [30; 1.4];
ub1 = [30; 1.4];

input_ref = Star(lb1, ub1);

N = 50;  % number of control steps

n_cores = 4; % number of cores 

t_gap = 1.4;
D_default = 10;

% safety specification: distance > alp * safe_distance

map = [1 0 0 -1 0 0; 0 0 0 0 1 0]; % get distance between two cars and velocity of ego car

% plot safe_distance vs. velocity
alp = 1;
map1 = [0 0 0 0 alp*t_gap 0; 0 0 0 0 1 0]; % safe distance and velocity of ego car

%% Reachability analysis of NNCS

[P1, reachTime1] = ncs.reach('approx-star', init_set, input_ref, n_cores, N);

% plot intermediate reach set (all reachable set of the plant) vs. time
reachSet = ncs.plant.intermediate_reachSet;
reachSet = [init_set reachSet]; % add init_set into the reachable set
times = 0:0.01:N*0.1;

dis = [];
safe_dis = [];
ego_vel = [];
lead_vel = [];

for i=1:length(reachSet)
    dis = [dis reachSet(i).affineMap([1 0 0 -1 0 0], [])];
    safe_dis = [safe_dis reachSet(i).affineMap([0 0 0 0 alp*t_gap 0], alp*D_default)];
    ego_vel = [ego_vel reachSet(i).affineMap([0 0 0 0 1 0], [])];
    lead_vel = [lead_vel reachSet(i).affineMap([0 1 0 0 0 0], [])];
end

%% Plot results

figure;
Star.plotRanges_2D(dis, 1, times, 'blue'); % plot distance between two cars versus time 
hold on;
Star.plotRanges_2D(safe_dis, 1, times, 'red'); % plot safe distance versus time
hold on;
xlabel('Time (s)');
ylabel('Distance (m)');
title('Actual distance (blue) vs. safe distance (red)');
grid;
set(gca,'FontSize',13)
saveas(gcf, 'reachSet.png');
