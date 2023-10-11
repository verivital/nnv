%% Safety verification using reachability analysis of ACC

%% Show scenario
imshow("ACC_scenario.png");
% Example from: https://www.mathworks.com/help/mpc/ug/adaptive-cruise-control-using-model-predictive-controller.html

%% Begin verification
% Load components
net = load_NN_from_mat('controller_5_20.mat');
% define step-size for reachability of the plant
reachStep = 0.05;
% define control step size for nncs
controlPeriod = 0.1;
% define output matrix (C) of the plant
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
% Create plant
plant = NonLinearODE(6,1,@dynamicsACC, reachStep, controlPeriod, output_mat);
% Create NNCS
nncs = NonlinearNNCS(net,plant);

% Final time (would be part of specifications)
tF = 5; % seconds
% define time vector (used for plotting later on)
time = 0:controlPeriod:5;
% number of control steps
steps = length(time);

% Define reference inputs to controller (fixed duing verification)
% 30 -> target speed
% 1.4- > time gap between vehicles
input_ref = [30;1.4];

% Initial set
% x = [x_ego, v_ego, a_ego, x_lead, v_lead, a_lead)
lb = [90; 32; 0; 10; 30; 0];
ub = [110; 32.2; 0; 11; 30.2; 0];
init_set = Star(lb,ub);

% Input set (no action to start with)
lb = 0;
ub = 0;
input_set = Star(lb,ub);

% Define reachabilty analysis parameters 
reachPRM.ref_input = input_ref; % reference input to controller
reachPRM.numSteps = steps-1; % number of control steps to compute 
reachPRM.init_set = init_set; % initial conditions
reachPRM.numCores = 1; % number of cores for reachability computation (controller)
reachPRM.reachMethod = 'approx-star'; % reach method for controller

% Compute reachability analysis
[R,rT] = nncs.reach(reachPRM);
disp("Time to compute ACC reach sets: " +string(rT));

%% Visualize results

% Transform reach set into actual distance vs safe distance
t_gap = 1.4;
D_default = 10;
outAll = [];
safe_dis = [];
% Transfrom intermediate reachsets from cora to NNV
nncs.plant.get_interval_sets();
% Get intermediate reach sets and transform them
allsets = nncs.plant.intermediate_reachSet;
for i=1:length(allsets)
    outAll = [outAll allsets(i).affineMap(output_mat,[])]; % [relative distance, relative velocity, ego-car velocity]
    safe_dis = [safe_dis allsets(i).affineMap([0 0 0 0 t_gap 0], D_default)]; % safety distance at every reach step
end
times = reachStep:reachStep:tF;

% Plot results
f = figure;
Star.plotRanges_2D(outAll,2,times,'b');
hold on;
Star.plotRanges_2D(safe_dis,1,times,'r');
xlabel('Time (s)');
ylabel('Distance (m)');
