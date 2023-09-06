%% NNCS

% Controller
load controller.mat;

weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = {};
for i=1:n - 1
    Layers{i} = LayerS(weights{1, i}, bias{1, i}, 'poslin');
end
Layers{n} = LayerS(weights{1, n}, bias{1, n}, 'purelin');
Controller = NN(Layers); % feedforward neural network controller
Controller.InputSize = 2;
Controller.OutputSize = 1;

% Plant
controlPeriod = 0.2;
reachStep = 0.01;
Plant = NonLinearODE(2, 1, @dynamics, reachStep, controlPeriod, eye(2));

% NNCS
feedbackMap = [0]; % feedback map, y[k] 
ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

%% Reachability analysis

N = 50; % number of control steps   
n_cores = 1; % number of cores

lb = [0.7; 0.4];
ub = [0.75; 0.45];
init_set = Star(lb, ub);
input_ref = [];

[simTrace, controlTrace] = ncs.evaluate(0.2, N, [0.7; 0.6], []);

% Reachability parameters
reachPRM.numSteps = N;
reachPRM.numCores = 1;
reachPRM.reachMethod = 'approx-star';
reachPRM.ref_input = input_ref;
reachPRM.init_set = init_set;

% Reachability computation
[P, reachTime] = ncs.reach(reachPRM);

% Visualize
fig = figure;
Star.plotBoxes_2D_noFill(P, 1, 2, 'blue');

