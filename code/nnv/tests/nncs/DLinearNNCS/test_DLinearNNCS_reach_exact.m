% Reachability analysis for Discrete Linear ACC model
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


% let x7 = -2*x(3) + 2 * a_lead -> x7(0) = -2*x(3)(0) + 2*alead
% -> dx7 = -2dx3


A = [0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 0 0 0 1; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 -2 0; 0 0 -2 0 0 0 0];
B = [0; 0; 0; 0; 0; 2; 0];
C = [1 0 0 -1 0 0 0; 0 1 0 0 -1 0 0; 0 0 0 0 1 0 0];  % feedback relative distance, relative velocity, longitudinal velocity
D = [0; 0; 0]; 

plant = LinearODE(A, B, C, D); % continuous plant model

plantd = plant.c2d(0.1); % discrete plant model

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

ncs = DLinearNNCS(Controller, plantd); % a discrete linear neural network control system


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


%% Reachability Analysis


reachPRM.init_set = init_set;
reachPRM.ref_input = ref_input;
reachPRM.numSteps = 10;
reachPRM.reachMethod = 'exact-star';
reachPRM.numCores = 1;
[R, reachTime] = ncs.reach(reachPRM);



%% Verification

% safety property: actual distance > alpha * safe distance <=> d = x1 - x4  > alpha * d_safe = alpha * (1.4 * v_ego + 10)

% usafe region: x1 - x4 <= alpha * (1.4 * v_ego + 10)
alpha = 1;
unsafe_mat = [1 0 0 -1 alpha*1.4 0 0];
unsafe_vec = alpha*10; 
unsafeRegion = HalfSpace(unsafe_mat, unsafe_vec);

[safe, counterExamples, verifyTime] = ncs.verify(reachPRM, unsafeRegion);

%% Plot reach sets

figure;
ncs.plotPlantReachSets('blue', 1); % plot position of lead car
hold on;
ncs.plotPlantReachSets('red',4); % plot position of ego car
title('Position reach sets of lead car (blue) and ego car (red)');

figure;
ncs.plotPlantReachSets('blue', 2); % plot velocity of lead car
hold on; 
ncs.plotPlantReachSets('red', 5); % plot velocity of ego car
title('Velocity reach sets of lead car (blue) and ego car (red)');

figure; 
ncs.plotControllerReachSets('green', 1); % plot control sets
title('Controller reach sets');

%% Plot output reach sets: actual distance vs. safe distance

% plot reachable set of the distance between two cars d = x1 - x4
figure; 
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
title('Actual Distance versus. Safe Distance');

