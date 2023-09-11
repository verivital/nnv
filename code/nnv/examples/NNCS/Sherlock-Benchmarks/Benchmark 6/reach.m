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
Controller.InputSize = 3;
Controller.OutputSize = 1;
% offset = 3;
% scale_factor = 1;

% Plant
reachStep = 0.02;
controlPeriod = 0.2;
output_mat = eye(3); % feedback 
Plant = NonLinearODE(3, 1, @dynamics, reachStep, controlPeriod, output_mat);

% NNCS
feedbackMap = [0]; % feedback map, y[k] 
ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system


%% Analysis
N = 4; % number of control steps   
lb = [0.35; -0.35; 0.35];
ub = [0.36; -0.34; 0.36];
% ub = [0.4; -0.3; 0.4];
init_set = Star(lb, ub);
input_ref = [];

% Evaluate
[simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, []);

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

