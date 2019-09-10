% This script does reachability analysis for mountain car benchmark and
% plot some reachable set

load MountainCar_ReluController.mat;
W = nnetwork.W; % weight matrices
b = nnetwork.b; % bias vectors

n = length(W);
Layers = [];
for i=1:n - 1
    L = Layer(W{1, i}, b{1, i}, 'ReLU');
    Layers = [Layers L];
end

L = Layer(W{1, n}, b{1, n}, 'Linear');

Layers = [Layers L];

Ts = 0.1; 
Controller = FFNN(Layers); % feedforward neural network controller
Plant = DNonLinearODE(2, 1, @discrete_car_dynamics, Ts); % two states and one input
output_mat = [1 0; 0 1];
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

feedbackMap = [0]; % feedback map, y[k] = [p[k]; v[k]], feedback both position and velocity with no delay 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial condition of the Plant (the car)
% x0 = [p0, v0], position and velocity

% This initial condition come from the paper: 
% Verisig: verifying safety properties of hybrid systems with neural network controllers, Radoslav Ivanov, HSCC2019, 

p0 = cell(13, 1);
p0{1,1} = [-0.41; -0.4];
p0{2,1} = [-0.415; -0.41];
p0{3,1} = [-0.42; -0.415];
p0{4,1} = [-0.43; -0.42];
p0{5,1} = [-0.45; -0.43];
p0{6,1} = [-0.48; -0.45];
p0{7,1} = [-0.50; -0.48];
p0{8,1} = [-0.53; -0.50];
p0{9,1} = [-0.55; -0.53];
p0{10,1} = [-0.57; -0.55];
p0{11,1} = [-0.58; -0.57];
p0{12,1} = [-0.59; -0.58];
p0{13,1} = [-0.6; -0.59];

v0 = [0; 0]; % range of initial velocity, v0 = 0

% there are 13 initial set of states 
reachTime = zeros(13, 1); % reachTime
safety_checkingTime = zeros(13, 1); % safety checking time 
verifyTime = zeros(13,1); % total verification time = reachTime + safety_checkingTime
reachSet = cell(13, 1);

input_ref = []; % empty input reference in this case study
N = 10;  % takes 20 seconds 
n_cores = 4; % number of cores 

for i=1:1 % test for first initial set
    init_pos = p0{i, 1};
    init_vel = v0;
    
    lb = [init_pos(1); init_vel(1)];
    ub = [init_pos(2); init_vel(2)];
    
    init_set = Star(lb, ub); 
    [reachSet{i, 1}, reachTime(i)] = ncs.reach('approx-star', init_set, input_ref, n_cores, N);
end

S = reachSet{1,1};
%Star.plots(S);
Star.plotBoxes_2D_noFill(S, 1, 2, 'blue');