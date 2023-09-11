% Reachability of spring electro control system

%% NNCS

% Controller
load controller.mat;
Layers = {};
Layers{1} = LayerS(W{1,1}, b{1,1}, 'poslin');
Layers{2} = LayerS(W{1,2}, b{1,2}, 'purelin');
Controller = NN(Layers); % feedforward neural network controller
Controller.InputSize = 3;
Controller.OutputSize = 1;

% Plant
load sysd.mat;
Plant = DLinearODE(sysd.A, sysd.B, sysd.C, sysd.D, 0.1);
feedbackMap = [0;1]; % feedback map, y[k] and y[k-1]

% NNCS
ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system


%% Reachability analysis
% Initial conditions
init_set = Star([0.8;-1],[1;-0.8]);
input_set = Star(0.99,1);
N = 5; % number of step

% reachability parameters
n_cores = 1; % number of cores 
reachPRM.numCores = n_cores;
reachPRM.init_set = init_set;
reachPRM.numSteps = N;
reachPRM.reachMethod = 'approx-star';
reachPRM.ref_input = input_set;
reachPRM2 = reachPRM;
reachPRM2.reachMethod = 'exact-star';

% Reachability computation
[S1, reachTime1] = ncs.reach(reachPRM);
[S2, reachTime2] = ncs.reach(reachPRM2);

%% Visualization
figure;
Star.plots(S1);
figure;
Star.plots(S2);