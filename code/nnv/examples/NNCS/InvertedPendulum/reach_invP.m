% Reachability analyisis of inverted pendulum NNCS

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
lb = [0.08; 0.2; 0; 0];
ub = [0.1; 0.4; 0; 0];
% init_set = Polyhedron('lb', lb, 'ub', ub);
init_set = Star(lb,ub);

N = 10; % number of step
n_cores = 1; % number of cores 
reachPRM.numSteps = 25;
reachPRM.numCores = n_cores;
reachPRM.ref_input = [];
reachPRM.reachMethod = 'approx-star';
reachPRM.init_set = init_set;

[P1, reachTime1] = ncs.reach(reachPRM);

% plot output (position x[1] and velocity x[2])
maps = [1 0 0 0; 0 1 0 0]; 

Pos_Vel_ReachSet = [];
for i=1:length(P1)
    Pos_Vel_ReachSet = [Pos_Vel_ReachSet P1(i).affineMap(maps,[])];
end

% Visualize
figure;
Star.plots(Pos_Vel_ReachSet);
