% Load network info
load network.mat;
% weights
W = network.weights;
% bias
b = network.bias;
n = length(b); % number of layers
Layers = cell(n,1);
for i=1:n - 1
    Li = LayerS(W{i}, b{i}, 'poslin'); % relu (hidden layers)
    Layers{i} = Li;
end
Layers{end} = LayerS(W{n}, b{n}, 'purelin'); % linear (output)

% Create network based on layer array
F = NN(Layers);

% Define input set
% x = [error, vref, vout, icurrent]
lb = [-3; 5; 5; 0];
ub = [3; 15; 15; 15];
I = Star(lb,ub);

% Define reachability parameters
reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
numCores = feature('numcores');
reachOptions.numCores = numCores;

% Compute reachability analsysis
t = tic;
R = F.reach(I, reachOptions); % exact reach set
toc(t);

% About 1-2 minutes to compute