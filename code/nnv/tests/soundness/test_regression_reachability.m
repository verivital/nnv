% test_regression_reachability
% Regression tests for reachability analysis
% Tests that reachable sets have expected properties
% To run: results = runtests('test_regression_reachability')

%% Test 1: Star through identity FC layer
rng(42);

% Identity matrix - should preserve Star
W = eye(3);
b = zeros(3, 1);
L = FullyConnectedLayer(W, b);

c = [1; 2; 3];
V = [c eye(3)];
C = [eye(3); -eye(3)];
d = ones(6, 1);

input_star = Star(V, C, d);
output_star = L.reach(input_star, 'approx-star');

% Center should be preserved
assert(max(abs(output_star.V(:,1) - c)) < 1e-10, 'Identity FC should preserve center');

%% Test 2: ImageStar through identity Conv2D
rng(42);

% Identity 1x1 filter - use 2 filters to avoid MATLAB dimension collapse
W = zeros(1, 1, 1, 2);
W(1, 1, 1, 1) = 1.0;
W(1, 1, 1, 2) = 1.0;  % Both filters are identity
b = zeros(1, 1, 2);  % 3D bias: [1, 1, numFilters]

L2 = Conv2DLayer(W, b, [0 0 0 0], [1 1], [1 1]);

V = zeros(3, 3, 1, 2);
V(:,:,1,1) = [1 2 3; 4 5 6; 7 8 9];
V(:,:,1,2) = rand(3, 3) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L2.reach(input_is, 'approx-star');

% Center should be preserved
center_diff = max(abs(output_is.V(:,:,1,1) - V(:,:,1,1)), [], 'all');
assert(center_diff < 1e-10, 'Identity Conv2D should preserve center');

%% Test 3: ReLU on all-positive Star
rng(42);

L3 = ReluLayer();

% Star entirely in positive region
c = [2; 3; 4];
V = [c 0.1*eye(3)];
C = [eye(3); -eye(3)];
d = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

input_star = Star(V, C, d);
output_star = L3.reach(input_star, 'approx-star');

% All-positive input through ReLU should be unchanged
assert(max(abs(output_star.V(:,1) - c)) < 1e-10, 'ReLU should be identity on positive region');

%% Test 4: ReLU on all-negative Star
rng(42);

L4 = ReluLayer();

% Star entirely in negative region
c = [-2; -3; -4];
V = [c 0.1*eye(3)];
C = [eye(3); -eye(3)];
d = [0.5; 0.5; 0.5; 0.5; 0.5; 0.5];

input_star = Star(V, C, d);
output_star = L4.reach(input_star, 'approx-star');

% All-negative input through ReLU should be zero
assert(max(abs(output_star.V(:,1))) < 1e-10, 'ReLU should output zero on negative region');

%% Test 5: Zonotope bounds through ReLU
rng(42);

L5 = ReluLayer();

% Create Zonotope
c = [0; 1; -1];  % mixed signs
G = 0.5 * eye(3);
Z = Zono(c, G);

% Convert to Star and reach
input_star = Z.toStar();
output_star = L5.reach(input_star, 'approx-star');

% Check bounds are valid (non-negative for ReLU output)
[lb, ub] = output_star.getRanges();
assert(all(lb >= -1e-10), 'ReLU output lower bounds should be >= 0');

%% Test 6: Star dimension preservation through layers
rng(42);

% FC layer that changes dimension
W = randn(5, 3);
b = randn(5, 1);
L6 = FullyConnectedLayer(W, b);

c = randn(3, 1);
V = [c randn(3, 2)];
C = [eye(2); -eye(2)];
d = ones(4, 1);

input_star = Star(V, C, d);
output_star = L6.reach(input_star, 'approx-star');

% Output dimension should be 5
assert(size(output_star.V, 1) == 5, 'FC output dimension should be 5');

% Number of predicates should be preserved
assert(output_star.nVar == input_star.nVar, 'Number of predicates should be preserved');

%% Test 7: ImageStar spatial dimension change through pooling
rng(42);

L7 = MaxPooling2DLayer('maxpool_dim', [2 2], [2 2], [0 0 0 0]);

V = zeros(4, 4, 1, 2);
V(:,:,1,1) = rand(4, 4);
V(:,:,1,2) = rand(4, 4) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L7.reach(input_is, 'approx-star');

% Output should be 2x2
assert(size(output_is.V, 1) == 2, 'MaxPool output height should be 2');
assert(size(output_is.V, 2) == 2, 'MaxPool output width should be 2');

%% Test 8: Reachability with exact vs approx comparison
rng(42);

L8 = ReluLayer();

% Small Star for exact reachability
c = [0.5; -0.5];
V = [c 0.3*eye(2)];
C = [eye(2); -eye(2)];
d = ones(4, 1);

input_star = Star(V, C, d);

% Exact should give tighter bounds than approx
output_exact = L8.reach(input_star, 'exact-star');
output_approx = L8.reach(input_star, 'approx-star');

% Both should be non-empty
assert(~isempty(output_exact), 'Exact should produce output');
assert(~isempty(output_approx), 'Approx should produce output');

