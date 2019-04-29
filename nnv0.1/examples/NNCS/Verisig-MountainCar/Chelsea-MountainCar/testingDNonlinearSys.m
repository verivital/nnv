% This script does reachability analysis for toy discrete time nonlinear
% ODE and
% plots some reachable set

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

Ts = 0.01;
%Plant.set_timeStep(Ts)
% two inputs (state dim) and one output (control dim)
Controller = FFNN(Layers); % feedforward neural network controller
Plant = DNonLinearODE(2, 1, @discrete_car_dynamics, Ts);
output_mat = [1 0; 0 1];
Plant.set_output_mat(output_mat); % Defines the outputs that are fed back to the controller
feedbackMap = [0]; % feedback map, y[k] = [p[k]; v[k]], feedback both position and velocity with no delay 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial conditions
n = 2;
p0 = cell(n, 1);
p0{1,1} = [-.4; -.2];
p0{2,1} = [-.2; 0];

v0 = [0; 0]; % range of initial velocity, v0 = 0

reachTime = zeros(n, 1); % reachTime
safety_checkingTime = zeros(n, 1); % safety checking time 
verifyTime = zeros(n,1); % total verification time = reachTime + safety_checkingTime
reachSet = cell(n, 1);

input_ref = []; % empty input reference in this case study
N = 20;  % takes 20 seconds 
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

fig = figure;
Star.plots(S);
title("Reachable set of Chelsea-MountainCar");











