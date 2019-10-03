load controller_3_20.mat;
weights = network.weights;
bias = network.bias;
n = length(weights);
Layers = [];
for i=1:n - 1
    L = Layer(weights{1, i}, bias{i, 1}, 'ReLU');
    Layers = [Layers L];
end

L = Layer(weights{1, n}, bias{n, 1}, 'Linear');

Layers = [Layers L];

Ts = 0.1;

Controller = FFNN(Layers); % feedforward neural network controller
Plant = NonLinearODE(6, 1, @car_dynamics);
Plant.set_timeStep(0.01); % time step for reachability analysis of the plant
Plant.set_tFinal(Ts); % Ts = 0.1, sampling time for control signal from neural network controller
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system


lb = [98; 32; 0; 10; 30; 0];
ub = [100; 32.2; 0; 11; 30.2; 0];

init_set = Box(lb, ub);
ref_input_set = Box([30; 1.4], [30; 1.4]);

n_steps = 50;
n_samples = 50; 

[sim_time, sim_traces, control_traces, sampled_init_states, sampled_ref_inputs] = ncs.sample(Ts, n_steps, init_set, ref_input_set, n_samples);


