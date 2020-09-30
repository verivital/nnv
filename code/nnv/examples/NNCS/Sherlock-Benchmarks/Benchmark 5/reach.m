load controller.mat;

weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = [];
for i=1:n - 1
    L = Layer(weights{1, i}, bias{1, i}, 'poslin');
    Layers = [Layers L];
end

L = Layer(weights{1, n}, bias{1, n}, 'purelin');

Layers = [Layers L];

Controller = FFNN(Layers); % feedforward neural network controller
offset = 10;
scale_factor = 1;
reachStep = 0.02;
controlPeriod = 0.2;
output_mat = eye(3);
Plant = NonLinearODE(2, 1, @dynamics, reachStep, controlPeriod, output_mat);

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial condition of the Plant


% reference input for neural network controller


N = 50; % number of control steps   

n_cores = 4; % number of cores

lb = [0.3; 0.3; -0.4];
ub = [0.4; 0.4; -0.3];
init_set = Star(lb, ub);
input_ref = [];

[simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, []);

[P, reachTime] = ncs.reach('approx-star', init_set, input_ref, n_cores, N);

fig = figure;
Star.plotBoxes_2D_noFill(P, 1, 2, 'blue');
saveas(fig, 'reachSet.pdf');
save result.mat reachTime ncs;
