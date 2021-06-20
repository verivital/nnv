load controller.mat;

weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = [];
for i=1:n - 1
    L = LayerS(weights{1, i}, bias{1, i}, 'poslin');
    Layers = [Layers L];
end

L = LayerS(weights{1, n}, bias{1, n}, 'purelin');

Layers = [Layers L];
Controller = FFNNS(Layers); % feedforward neural network controller
offset = 100;
scale_factor = 0.1;

reachStep = 0.02;
controlPeriod = 0.5;
output_mat = eye(3); % feedback 
Plant = NonLinearODE(3, 1, @dynamics, reachStep, controlPeriod, output_mat);

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% N = 50; % number of control steps   
N = 2;
n_cores = 4; % number of cores

% lb = [0.35; 0.45; 0.25];
lb = [0.445; 0.545; 0.345];
ub = [0.45; 0.55; 0.35];
init_set = Star(lb, ub);
input_ref = [];

[simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, []);

reachPRM.numSteps = N;
reachPRM.numCores = n_cores;
reachPRM.reachMethod = 'approx-star';
reachPRM.ref_input = input_ref;
reachPRM.init_set = init_set;

[P, reachTime] = ncs.reach(reachPRM);

fig = figure;
Star.plotBoxes_2D_noFill(P, 1, 2, 'blue');
saveas(fig, 'reachSet.pdf');
save result.mat reachTime ncs;
