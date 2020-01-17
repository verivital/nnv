%% Script to generate the graphs for DESTION 2020
% We are going to evaluate the adaptive cruise control system with
% different controllers and different scenarios. We will also evaluate the
% discrete ACC and the nonlinear dynamics.
clc;clear;close all;

%% Load all components and create NNCS
load controller_5_20.mat;

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

Controller = FFNN(Layers); % feedforward neural network controller
Plant = NonLinearODE(6, 1, @car_dynamics);
Plant.set_timeStep(0.01); % time step for reachability analysis of the plant
Plant.set_tFinal(0.1); % Ts = 0.1, sampling time for control signal from neural network controller
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system


%% Initialize components and set conditions

% initial position of lead car x_lead
x_lead = [90 92];
% initial condition of v_lead
% v_lead = cell(6, 1);
% v_lead{1, 1} = [29 30];
% v_lead{2, 1} = [28 29];
% v_lead{3, 1} = [27 28];
% v_lead{4, 1} = [26 27];
% v_lead{5, 1} = [25 26];
% v_lead{6, 1} = [24 25];
v_lead = [25.5 26.5];

% initial condition of x_internal_lead
internal_acc_lead = [0 0];
% initial condition of x_ego
x_ego = [10 11]; 
% initial condition of v_ego
v_ego = [30 30.2];
% initial condition of x_internal_ego
internal_acc_ego = [0 0];

% n = length(v_lead);

% init_set(n) = Star();


lb = [x_lead(1); v_lead(1); internal_acc_lead(1); x_ego(1); v_ego(1); internal_acc_ego(1)];
ub = [x_lead(2); v_lead(2); internal_acc_lead(2); x_ego(2); v_ego(2); internal_acc_ego(2)];
init_set(i) = Star(lb, ub);


% reference input for neural network controller
% t_gap = 1.4; v_set = 30;

lb = [30; 1.4];
ub = [30; 1.4];

input_ref = Star(lb, ub); % input reference set

%% Start reachability analysis

N = 50; % number of control steps 
times = 0:0.01:N*0.1;

n_cores = 4; % number of cores 

% safety specification: relative distance > safe distance
% dis = x_lead - x_ego  
% safe distance between two cars, see here 
% https://www.mathworks.com/help/mpc/examples/design-an-adaptive-cruise-control-system-using-model-predictive-control.html
% dis_safe = D_default + t_gap * v_ego;

t_gap = 1.4;
D_default = 10;

% safety specification: x_lead - x_ego - t_gap * v_ego - D_default > 0
% unsafe region: x_lead - x_ego - t_gap * v_ego <= D_default 

unsafe_mat = [1 0 0 -1 -t_gap 0];
unsafe_vec = [D_default];

safe = cell(1, n); % safety status
VT = zeros(1,n); % verification time
counterExamples = cell(1,n); % counter examples
ReachSystem = [];
t = tic;
ReachSystem = ncs.reach('approx-star', init_set(i), input_ref, n_cores, N);
[safe1, ~] = ncs.check_safety(unsafe_mat, unsafe_vec, n_cores);
if safe1 == 1
    safe{i} = 'SAFE';
else
    Bi = init_set(i).getBox;
    Br = input_ref.getBox;
    [rs, ~, counterExamples, ~, ~, ~] = ncs.falsify(0.1, N, Bi, Br, unsafe_mat, unsafe_vec, 1000);

    if rs == 1
        safe{i} = 'UNSAFE';
    else
        safe{i} = 'UNKNOWN';
    end
end
VT(i) = toc(t);

alp = 1;
t_gap = 1.4;
D_default = 10;
dis = [];
safe_dis = [];
ego_vel = [];
lead_vel = [];

for i=1:length(ReachSystem)
    dis = [dis ReachSystem(i).affineMap([1 0 0 -1 0 0], [])];
    safe_dis = [safe_dis ReachSystem(i).affineMap([0 0 0 0 alp*t_gap 0], alp*D_default)];
    ego_vel = [ego_vel ReachSystem(i).affineMap([0 0 0 0 1 0], [])];
    lead_vel = [lead_vel ReachSystem(i).affineMap([0 1 0 0 0 0], [])];
end


%% Plot results
figure;
Star.plotBoxes_2D_noFill(ReachSystem, 1, 4,'b');
figure;
Star.plotBoxes_2D_noFill(ReachSystem, 2, 5,'b');
figure;
Star.plotBoxes_2D_noFill(ego_vel, 1, times,'b');
figure;
Star.plotBoxes_2D_noFill(dis, 2, times,'b');
figure;
Star.plotBoxes_2D_noFill(safe_dis, 1, times,'b');
figure;
Star.plotBoxes_2D_noFill(lead_vel, 2, times,'b');

