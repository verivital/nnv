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
timeStep = 0.01;
tFinal = Ts;
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
Plant = NonLinearODE(6, 1, @car_dynamics, timeStep, tFinal, output_mat);

feedbackMap = [0]; % feedback map, y[k] 

ncs = NNCS(Controller, Plant, feedbackMap); % the neural network control system

x0 = [98; 32; 0; 10; 30; 0];
ref_input = [30; 1.4];

N = 50;
[simTrace, controlTrace] = ncs.evaluate(Ts, N, x0, ref_input);