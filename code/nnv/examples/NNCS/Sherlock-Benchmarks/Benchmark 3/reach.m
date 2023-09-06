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

%Plant
controlPeriod = 0.1;
reachStep = 0.01;
Plant = NonLinearODE(2, 1, @dynamics, reachStep, controlPeriod,eye(2));

% NNCS
feedbackMap = [0]; % feedback map, y[k] 
ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

%% Evaluation
N = 25; % number of control steps   

n_cores = 4; % number of cores

lb = [0.8; 0.4];
ub = [0.9; 0.5];
init_set = Star(lb, ub);
input_ref = [];

[simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, []);

plot(simTrace(1, :), simTrace(2, :));