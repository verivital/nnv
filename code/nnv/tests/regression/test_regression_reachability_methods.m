% test_regression_reachability_methods
% Regression test for comparing different reachability methods
% Tests exact-star vs approx-star methods
% Based on: examples/Tutorial/NN/compareReachability/reach_exact_vs_approx.m
% To run: results = runtests('test_regression_reachability_methods')

%% Test 1: Load pre-trained network
rng(42);

% Get path to neural network
net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'compareReachability', filesep, 'NeuralNetwork7_3.mat'];

net_data = load(net_path);

assert(isfield(net_data, 'W'), 'Should have weight matrices');
assert(isfield(net_data, 'b'), 'Should have bias vectors');

%% Test 2: Create NNV model from raw weights
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'compareReachability', filesep, 'NeuralNetwork7_3.mat'];
net_data = load(net_path);

W = net_data.W;
b = net_data.b;
n = length(b);

Layers = cell(n, 1);
% Hidden layers are all ReLUs
for i = 1:n-1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers{i} = Li;
end
% Output layer is linear
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');
Layers{end} = Ln;

% Create NN model
F = NN(Layers);

assert(~isempty(F), 'Network should be created');

%% Test 3: Define input set
rng(42);

lb = [0; 0; 0];
ub = [1; 1; 1];
I = Star(lb, ub);

assert(~isempty(I), 'Input Star should be created');
assert(I.dim == 3, 'Input dimension should be 3');

%% Test 4: Exact-star reachability
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'compareReachability', filesep, 'NeuralNetwork7_3.mat'];
net_data = load(net_path);
W = net_data.W;
b = net_data.b;
n = length(b);

Layers = cell(n, 1);
for i = 1:n-1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers{i} = Li;
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');
Layers{end} = Ln;
F = NN(Layers);

lb = [0; 0; 0];
ub = [1; 1; 1];
I = Star(lb, ub);

reachOptions = struct;
reachOptions.reachMethod = 'exact-star';

Re = F.reach(I, reachOptions);

assert(~isempty(Re), 'Exact reachable set should be computed');
assert(length(Re) >= 1, 'Should have at least one output set');

%% Test 5: Approx-star reachability
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'compareReachability', filesep, 'NeuralNetwork7_3.mat'];
net_data = load(net_path);
W = net_data.W;
b = net_data.b;
n = length(b);

Layers = cell(n, 1);
for i = 1:n-1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers{i} = Li;
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');
Layers{end} = Ln;
F = NN(Layers);

lb = [0; 0; 0];
ub = [1; 1; 1];
I = Star(lb, ub);

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';

Ra = F.reach(I, reachOptions);

assert(~isempty(Ra), 'Approx reachable set should be computed');
assert(length(Ra) == 1, 'Approx method should return single set');

%% Test 6: Approx-zono reachability
% NOTE: approx-zono requires Zono input, not Star
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'compareReachability', filesep, 'NeuralNetwork7_3.mat'];
net_data = load(net_path);
W = net_data.W;
b = net_data.b;
n = length(b);

Layers = cell(n, 1);
for i = 1:n-1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers{i} = Li;
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');
Layers{end} = Ln;
F = NN(Layers);

% Create Zono input (approx-zono requires Zono, not Star)
lb = [0; 0; 0];
ub = [1; 1; 1];
c = (lb + ub) / 2;  % center
V = diag((ub - lb) / 2);  % generators
I = Zono(c, V);

reachOptions = struct;
reachOptions.reachMethod = 'approx-zono';

Rz = F.reach(I, reachOptions);

assert(~isempty(Rz), 'Zono reachable set should be computed');

%% Test 7: Soundness - concrete outputs within exact bounds
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'compareReachability', filesep, 'NeuralNetwork7_3.mat'];
net_data = load(net_path);
W = net_data.W;
b = net_data.b;
n = length(b);

Layers = cell(n, 1);
for i = 1:n-1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers{i} = Li;
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');
Layers{end} = Ln;
F = NN(Layers);

lb = [0; 0; 0];
ub = [1; 1; 1];
I = Star(lb, ub);

reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
Re = F.reach(I, reachOptions);

% Compute bounds from union of exact sets
% Get output dimension from first reach set
[temp_lb, ~] = Re(1).getRanges;
out_dim = length(temp_lb);
all_lb = inf(out_dim, 1);
all_ub = -inf(out_dim, 1);
for i = 1:length(Re)
    [lb_i, ub_i] = Re(i).getRanges;
    all_lb = min(all_lb, lb_i);
    all_ub = max(all_ub, ub_i);
end

% Test corners and center
y_lb = F.evaluate(lb);
y_ub = F.evaluate(ub);
y_center = F.evaluate((lb + ub) / 2);

tol = 1e-5;
assert(all(y_lb >= all_lb - tol) && all(y_lb <= all_ub + tol), ...
    'LB output should be within exact bounds');
assert(all(y_ub >= all_lb - tol) && all(y_ub <= all_ub + tol), ...
    'UB output should be within exact bounds');
assert(all(y_center >= all_lb - tol) && all(y_center <= all_ub + tol), ...
    'Center output should be within exact bounds');

%% Test 8: Approx overapproximates exact
rng(42);

net_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'examples', filesep, 'Tutorial', filesep, ...
    'NN', filesep, 'compareReachability', filesep, 'NeuralNetwork7_3.mat'];
net_data = load(net_path);
W = net_data.W;
b = net_data.b;
n = length(b);

Layers = cell(n, 1);
for i = 1:n-1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers{i} = Li;
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');
Layers{end} = Ln;
F = NN(Layers);

lb = [0; 0; 0];
ub = [1; 1; 1];
I = Star(lb, ub);

% Compute exact
reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
Re = F.reach(I, reachOptions);

% Compute approx
reachOptions.reachMethod = 'approx-star';
Ra = F.reach(I, reachOptions);

% Get bounds from exact (union)
[temp_lb, ~] = Re(1).getRanges;
out_dim = length(temp_lb);
exact_lb = inf(out_dim, 1);
exact_ub = -inf(out_dim, 1);
for i = 1:length(Re)
    [lb_i, ub_i] = Re(i).getRanges;
    exact_lb = min(exact_lb, lb_i);
    exact_ub = max(exact_ub, ub_i);
end

% Get bounds from approx
[approx_lb, approx_ub] = Ra.getRanges;

% Approx should overapproximate: approx_lb <= exact_lb, approx_ub >= exact_ub
tol = 1e-5;
assert(all(approx_lb <= exact_lb + tol), ...
    'Approx lower bound should be <= exact lower bound (overapproximation)');
assert(all(approx_ub >= exact_ub - tol), ...
    'Approx upper bound should be >= exact upper bound (overapproximation)');

