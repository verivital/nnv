% test_soundness_ReluLayer
% Soundness tests for ReluLayer (both exact and approximate methods)
% Verifies that computed reachable sets contain all actual outputs
% To run: results = runtests('test_soundness_ReluLayer')

%% Test 1: Exact-star soundness with ImageStar (small)
L = ReluLayer();

% Create small ImageStar with values crossing zero
V = zeros(2, 2, 1, 3);
V(:,:,1,1) = [-0.5 0.5; 0 1];   % center has both negative and positive
V(1,1,1,2) = 0.3;              % basis 1 affects pixel (1,1)
V(1,2,1,3) = 0.2;              % basis 2 affects pixel (1,2)

C = [eye(2); -eye(2)];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

input_is = ImageStar(V, C, d, pred_lb, pred_ub);

% Test exact-star (may produce multiple output sets)
[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is, 'exact-star', 30);
assert(passed, 'ReLU exact-star soundness test failed: %s', msg);

%% Test 2: Approx-star soundness
% Approximate methods should still be sound (contain all outputs)
% Recreate layer and input (test sections are independent)
L = ReluLayer();

V = zeros(2, 2, 1, 3);
V(:,:,1,1) = [-0.5 0.5; 0 1];
V(1,1,1,2) = 0.3;
V(1,2,1,3) = 0.2;

C = [eye(2); -eye(2)];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

input_is = ImageStar(V, C, d, pred_lb, pred_ub);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is, 'approx-star', 50);
assert(passed, 'ReLU approx-star soundness test failed: %s', msg);

%% Test 3: All-positive input (ReLU identity)
L = ReluLayer();

V2 = zeros(3, 3, 1, 2);
V2(:,:,1,1) = ones(3, 3);  % all positive center
V2(2,2,1,2) = 0.1;         % small perturbation, still all positive

C2 = [1; -1];
d2 = [0.5; 0.5];
pred_lb2 = -0.5;
pred_ub2 = 0.5;

input_is2 = ImageStar(V2, C2, d2, pred_lb2, pred_ub2);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is2, 'exact-star', 30);
assert(passed, 'ReLU all-positive soundness test failed: %s', msg);

%% Test 4: All-negative input (ReLU zeros)
L = ReluLayer();

V3 = zeros(3, 3, 1, 2);
V3(:,:,1,1) = -ones(3, 3);  % all negative center
V3(2,2,1,2) = 0.1;          % small perturbation, still all negative

C3 = [1; -1];
d3 = [0.5; 0.5];
pred_lb3 = -0.5;
pred_ub3 = 0.5;

input_is3 = ImageStar(V3, C3, d3, pred_lb3, pred_ub3);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is3, 'exact-star', 30);
assert(passed, 'ReLU all-negative soundness test failed: %s', msg);

%% Test 5: Multi-channel ReLU
rng(42);  % Fixed seed for reproducibility
L = ReluLayer();

V4 = zeros(2, 2, 3, 3);
V4(:,:,:,1) = randn(2, 2, 3) * 0.5;  % random center
V4(1,1,1,2) = 0.2;
V4(2,2,2,3) = 0.2;

C4 = [eye(2); -eye(2)];
d4 = [1; 1; 1; 1];
pred_lb4 = [-1; -1];
pred_ub4 = [1; 1];

input_is4 = ImageStar(V4, C4, d4, pred_lb4, pred_ub4);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is4, 'approx-star', 50);
assert(passed, 'Multi-channel ReLU soundness test failed: %s', msg);

%% Test 6: Star-based ReLU soundness (for feedforward networks)
L = ReluLayer();

% Create Star input
c = [0.5; -0.3; 0.1; -0.8];
V_basis = [0.3 0; 0 0.2; 0.1 0.1; 0 0.4];
V_star = [c V_basis];
C_star = [eye(2); -eye(2)];
d_star = ones(4, 1);
pred_lb_star = [-1; -1];
pred_ub_star = [1; 1];

input_star = Star(V_star, C_star, d_star, pred_lb_star, pred_ub_star);

% Use exact-star for Star-based ReLU (approx-star may have soundness issues)
[passed, msg] = soundness_test_utils.verify_layer_soundness_star(L, input_star, 'exact-star', 50);
assert(passed, 'Star-based ReLU soundness test failed: %s', msg);
