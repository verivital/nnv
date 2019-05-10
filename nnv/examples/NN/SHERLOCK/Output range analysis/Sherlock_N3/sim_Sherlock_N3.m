load neural_network_information_3.mat
n = length(W);
L = [];
for i=1:n
    L1 = Layer(W{1, i}, b{1, i}, 'ReLU');
    L = [L L1];
end

F = FFNN(L); % neural network

lb = [0; 0];
ub = [10; 10];

I = Box(lb, ub); % input set = Sherlock
num_samples = 5000; % use 5000 simulations to estimate ranges of the output

t = tic;
sim_range = F.estimate_ranges(I, num_samples); % using simulation to estimate output range
sim_time = toc(t);
save sim_range sim_range;
save sim_time sim_time;
