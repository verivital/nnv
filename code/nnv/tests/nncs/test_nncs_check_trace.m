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
reachTimeStep = 0.01;
controlPeriod = 0.1;
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Controller = FFNN(Layers); % feedforward neural network controller
Plant = NonLinearODE(6, 1, @car_dynamics, reachTimeStep, controlPeriod, output_mat);

ncs = NNCS(Controller, Plant); % the neural network control system

x0 = [98; 32; 0; 10; 30; 0];
ref_input = [30; 1.4];

N = 50;
[simTrace, controlTrace] = ncs.evaluate(Ts, N, x0, ref_input);

% unsafe region: dis = x1 - x2 <= 80
unsafe_mat = [1 0 0 -1 0 0];
unsafe_vec = [78]; 

violate = NNCS.check_trace(simTrace, unsafe_mat, unsafe_vec);

fig = figure;
dis = [1 0 0 -1 0 0] * simTrace;
times = 0:Ts:N*Ts;
plot(times, dis, '-');
xlabel('t');
ylabel('distance');
