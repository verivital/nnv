%% Reachability analyisis of inverted pendulum NNCS

% inverted pendulum dynamics (linear dynamics, 4 states)
%
% - System states
%  - x1 = x position
%  - x2 = pendulum angle
%  - x3 = velocity (x)
%  - x4 = angular velocity (pendulum)

%% NNCS
load controller.mat;
load sysd.mat;

% Controller
L1 = LayerS(W{1,1}, b{1,1}, 'poslin');
L2 = LayerS(W{1,2}, b{2,1}, 'purelin');
Controller = NN({L1;L2});
Controller.InputSize = 4;
Controller.OutputSize = 1;

% Plant
Plant = DLinearODE(A, B, C, D, Ts);

% NNCS
feedbackMap = [0]; % feedback map, y[k] 
ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

%% Reachability analysis
% initial set of state of the plant x = [-0.1 <= x[1] <= -0.05, 0.85 <= x[2] <= 0.9, x[3] = 0; x[4] = 0]
lb = [-0.1; 0.85; 0; 0];
ub = [-0.05; 0.9; 0; 0];
init_set = Star(lb,ub);

n_cores = 1; % number of cores 
reachPRM.numSteps = 35; % number of control steps
reachPRM.numCores = n_cores;
reachPRM.ref_input = [];
reachPRM.reachMethod = 'approx-star';
reachPRM.init_set = init_set;

[P1, reachTime] = ncs.reach(reachPRM);

% plot output (position x[1] and velocity x[2])
maps = [1 0 0 0; 0 1 0 0]; 

Pos_Vel_ReachSet = [];
for i=1:length(P1)
    Pos_Vel_ReachSet = [Pos_Vel_ReachSet P1(i).affineMap(maps,[])];
end

% Visualize
t = tic;
figure;
Star.plots(Pos_Vel_ReachSet);
toc(t)

