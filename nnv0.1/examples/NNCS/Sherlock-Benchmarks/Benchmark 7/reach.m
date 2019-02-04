load controller.mat;

weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = [];
for i=1:n - 1
    L = Layer(weights{1, i}, bias{1, i}, 'ReLU');
    Layers = [Layers L];
end

L = Layer(weights{1, n}, bias{1, n}, 'Linear');

Layers = [Layers L];

Controller = FFNN(Layers); % feedforward neural network controller
Plant = NonLinearODE(2, 1, @dynamics);
Plant.set_timeStep(0.02); % time step for reachability analysis of the plant
Plant.set_tFinal(0.2); % Ts = 0.2, sampling time for control signal from neural network controller
output_mat = [1 0; 0 1]; % feedback 
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

% initial condition of the Plant


% reference input for neural network controller


N = 50; % number of control steps   

n_cores = 4; % number of cores

lb = [0.35; 0.45; 0.25];
ub = [0.45; 0.55; 0.35];
init_set = Star(lb, ub);
input_ref = [];

[simTrace, controlTrace] = ncs.evaluate(0.2, N, lb, []);

[P, reachTime] = ncs.reach('approx-star', init_set, input_ref, n_cores, N);

fig = figure;
Star.plotBoxes_2D_noFill(P, 1, 2, 'blue');
saveas(fig, 'reachSet.pdf');
save result.mat reachTime ncs;
