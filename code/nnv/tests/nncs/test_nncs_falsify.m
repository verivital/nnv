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

Ts = 0.1; % control period

Controller = FFNN(Layers); % feedforward neural network controller
Plant = NonLinearODE(6, 1, @car_dynamics);
Plant.set_timeStep(0.01); % time step for reachability analysis of the plant
Plant.set_tFinal(Ts); % Ts = 0.1, sampling time for control signal from neural network controller, this is also control period of the controller
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Plant.set_output_mat(output_mat); % Define the outputs that is feedback to the controller

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system


lb = [98; 32; 0; 10; 30; 0];
ub = [100; 32.2; 0; 11; 30.2; 0];

init_set = Box(lb, ub);
ref_input_set = Box([30; 30], [1.4; 1.4]);

n_steps = 50;
n_samples = 50; 

% consider ACC case study to understand more
t_gap = 1.4;
D_default = 10;

% safety specification dis > safe_dis

unsafe_mat = [1 0 0 -1 -t_gap 0];
unsafe_vec = [D_default];

[falsify_result, falsify_time, counter_sim_traces, counter_control_traces, counter_init_states, counter_ref_inputs] = ncs.falsify(Ts, n_steps, init_set, ref_input_set, unsafe_mat, unsafe_vec, n_samples);

