% Reachability analysis of benchmark 8

%% NNCS

% Controller
load controller.mat;

weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = {};
for i=1:n-1
    Layers{i} = LayerS(weights{1, i}, bias{1, i}, 'poslin');
end
Layers{n} = LayerS(weights{1, n}, bias{1, n}, 'purelin');
Controller = NN(Layers); % feedforward neural network controller
Controller.InputSize = 4;
Controller.OutputSize = 1;

% Plant
reachStep = 0.02; % time step for reachability analysis of the plant
controlPeriod = 0.2; % control step
output_mat = eye(4); % feedback 
Plant = NonLinearODE(4, 1, @dynamics, reachStep, controlPeriod, output_mat);

% NNCS
feedbackMap = [0]; % feedback map, y[k] 
ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system


%% Analysis

% Initial conditions
N = 3;
lb = [0.5; 0.5; 0.5; 0.5];
ub = [0.6; 0.6; 0.6; 0.6];
init_set = Star(lb, ub);
input_ref = [];

% Evaluation
[simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, []);

% Reachability parameters
reachPRM.numSteps = N;
reachPRM.numCores = 1;
reachPRM.reachMethod = 'approx-star';
reachPRM.ref_input = input_ref;
reachPRM.init_set = init_set;

% Reach computation
[P, reachTime] = ncs.reach(reachPRM);

%% Visualization

figure;
Star.plotBoxes_2D_noFill(P, 1, 2, 'blue');
