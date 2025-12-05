% test_soundness_FullyConnectedLayer
% Soundness tests for FullyConnectedLayer
% Verifies that computed reachable sets contain all actual outputs
% To run: results = runtests('test_soundness_FullyConnectedLayer')

%% Test 1: Basic soundness with Star input
% Create a small fully connected layer
W = [1 2 -1; 0.5 -1 2; 1 1 1];  % 3x3 weight matrix
b = [0.1; -0.1; 0];             % 3x1 bias
L = FullyConnectedLayer(W, b);

% Create Star input with 2 predicate variables
c = [1; 0; -1];           % center
V_basis = [0.5 0; 0 0.3; 0.1 0.2];  % basis vectors
V = [c V_basis];
C = [1 0; -1 0; 0 1; 0 -1];  % box constraints
d = [1; 1; 1; 1];            % alpha in [-1, 1]^2
pred_lb = [-1; -1];
pred_ub = [1; 1];
input_star = Star(V, C, d, pred_lb, pred_ub);

% Verify soundness
[passed, msg] = soundness_test_utils.verify_layer_soundness_star(L, input_star, 'exact-star', 100);
assert(passed, 'Soundness test failed: %s', msg);

%% Test 2: Larger network soundness
% 5-input, 4-output layer
W2 = randn(4, 5);
b2 = randn(4, 1);
L2 = FullyConnectedLayer(W2, b2);

% Create Star input with 3 predicate variables
c2 = randn(5, 1);
V_basis2 = 0.1 * randn(5, 3);
V2 = [c2 V_basis2];
C2 = [eye(3); -eye(3)];
d2 = 0.5 * ones(6, 1);
pred_lb2 = -0.5 * ones(3, 1);
pred_ub2 = 0.5 * ones(3, 1);
input_star2 = Star(V2, C2, d2, pred_lb2, pred_ub2);

[passed, msg] = soundness_test_utils.verify_layer_soundness_star(L2, input_star2, 'exact-star', 100);
assert(passed, 'Soundness test failed for larger layer: %s', msg);

%% Test 3: Edge case - single predicate variable
W3 = [2; -1; 0.5];  % 3x1
b3 = [0; 0; 0];
L3 = FullyConnectedLayer(W3, b3);

c3 = [1];
V_basis3 = [0.2];
V3 = [c3 V_basis3];
C3 = [1; -1];
d3 = [1; 1];
pred_lb3 = -1;
pred_ub3 = 1;
input_star3 = Star(V3, C3, d3, pred_lb3, pred_ub3);

[passed, msg] = soundness_test_utils.verify_layer_soundness_star(L3, input_star3, 'exact-star', 50);
assert(passed, 'Soundness test failed for single predicate: %s', msg);

%% Test 4: Zero-centered input
W4 = [1 -1; -1 1];
b4 = [0; 0];
L4 = FullyConnectedLayer(W4, b4);

c4 = [0; 0];
V_basis4 = [1 0; 0 1];
V4 = [c4 V_basis4];
C4 = [eye(2); -eye(2)];
d4 = ones(4, 1);
pred_lb4 = [-1; -1];
pred_ub4 = [1; 1];
input_star4 = Star(V4, C4, d4, pred_lb4, pred_ub4);

[passed, msg] = soundness_test_utils.verify_layer_soundness_star(L4, input_star4, 'exact-star', 100);
assert(passed, 'Soundness test failed for zero-centered input: %s', msg);
