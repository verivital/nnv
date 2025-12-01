% test_soundness_relax_star
% Tests for relax-star reachability method
% To run: results = runtests('test_soundness_relax_star')

%% Test 1: ReLU with relax-star-range method
rng(42);

L = ReluLayer();

% Create Star input
c = randn(5, 1);
V_basis = randn(5, 2) * 0.5;
V = [c V_basis];
C = [eye(2); -eye(2)];
d = ones(4, 1);

input_star = Star(V, C, d, [-1; -1], [1; 1]);

% Test relax-star-range method
output = L.reach(input_star, 'relax-star-range');
assert(~isempty(output), 'relax-star-range should produce output');

%% Test 2: ReLU with relax-star-area method
rng(42);

L2 = ReluLayer();

c = randn(5, 1);
V_basis = randn(5, 2) * 0.5;
V = [c V_basis];
C = [eye(2); -eye(2)];
d = ones(4, 1);

input_star = Star(V, C, d, [-1; -1], [1; 1]);

% Test relax-star-area method
output = L2.reach(input_star, 'relax-star-area');
assert(~isempty(output), 'relax-star-area should produce output');

%% Test 3: ReLU with relax-star-bound method
rng(42);

L3 = ReluLayer();

c = randn(5, 1);
V_basis = randn(5, 2) * 0.5;
V = [c V_basis];
C = [eye(2); -eye(2)];
d = ones(4, 1);

input_star = Star(V, C, d, [-1; -1], [1; 1]);

% Test relax-star-bound method
output = L3.reach(input_star, 'relax-star-bound');
assert(~isempty(output), 'relax-star-bound should produce output');

%% Test 4: LeakyReLU with relax-star-dis method
rng(42);

% LeakyReluLayer constructor: (name, numInputs, inputNames, numOutputs, outputNames, gamma)
L4 = LeakyReluLayer('leaky_relax', 1, {'in'}, 1, {'out'}, 0.01);

c = randn(5, 1);
V_basis = randn(5, 2) * 0.5;
V = [c V_basis];
C = [eye(2); -eye(2)];
d = ones(4, 1);

input_star = Star(V, C, d, [-1; -1], [1; 1]);

% Test relax-star-dis method for LeakyReLU
output = L4.reach(input_star, 'relax-star-dis');
assert(~isempty(output), 'LeakyReLU relax-star-dis should produce output');

%% Test 5: FullyConnected with relax-star-range
rng(42);

W = randn(3, 5) * 0.5;
b = randn(3, 1) * 0.1;
L5 = FullyConnectedLayer(W, b);

c = randn(5, 1);
V_basis = randn(5, 2) * 0.5;
V = [c V_basis];
C = [eye(2); -eye(2)];
d = ones(4, 1);

input_star = Star(V, C, d, [-1; -1], [1; 1]);

% FC layer should work with relax-star methods
output = L5.reach(input_star, 'relax-star-range');
assert(~isempty(output), 'FC with relax-star-range should produce output');

%% Test 6: Conv2D with relax-star-area
rng(42);

W = randn(3, 3, 1, 2) * 0.1;
b = zeros(1, 1, 2);
L6 = Conv2DLayer(W, b, [0 0 0 0], [1 1], [1 1]);

% V has 2 generators (indices 2 and 3 in 4th dim), so need 2 predicates
V = zeros(5, 5, 1, 3);
V(:,:,1,1) = rand(5, 5) * 0.5;  % center
V(:,:,1,2) = rand(5, 5) * 0.1;  % generator 1
V(:,:,1,3) = rand(5, 5) * 0.1;  % generator 2

% C and d define constraints on 2 predicates
C = [eye(2); -eye(2)];
d = ones(4, 1);
input_is = ImageStar(V, C, d, [-1; -1], [1; 1]);

% Conv2D should work with relax-star methods
output = L6.reach(input_is, 'relax-star-area');
assert(~isempty(output), 'Conv2D with relax-star-area should produce output');

