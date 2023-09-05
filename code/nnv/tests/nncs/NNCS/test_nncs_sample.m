load('../controller_3_20.mat');
weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = [];
for i=1:n - 1
    L = LayerS(weights{1, i}, bias{i, 1}, 'poslin');
    Layers{i} = L;
end

Layers{n} = LayerS(weights{1, n}, bias{n, 1}, 'purelin');

Ts = 0.1;
reachTimeStep = 0.01;
controlPeriod = 0.1;
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Controller = NN(Layers); % feedforward neural network controller
Controller.InputSize = 5;
Controller.OutputSize = 1;
Plant = NonLinearODE(6, 1, @car_dynamics, reachTimeStep, controlPeriod, output_mat);

ncs = NNCS(Controller, Plant); % the neural network control system


lb = [98; 32; 0; 10; 30; 0];
ub = [100; 32.2; 0; 11; 30.2; 0];

init_set = Box(lb, ub);
ref_input_set = Box([30; 1.4], [30; 1.4]);

n_steps = 50;
n_samples = 50; 

[sim_time, sim_traces, control_traces, sampled_init_states, sampled_ref_inputs] = ncs.sample(Ts, n_steps, init_set, ref_input_set, n_samples);


