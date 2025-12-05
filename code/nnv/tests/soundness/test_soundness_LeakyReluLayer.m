% test_soundness_LeakyReluLayer
% Soundness tests for LeakyReluLayer (piecewise linear activation)
% Verifies that computed reachable sets contain all actual outputs
% To run: results = runtests('test_soundness_LeakyReluLayer')

%% Test 1: Basic LeakyReLU soundness with ImageStar
% LeakyReluLayer(name, numInputs, inputNames, numOutputs, outputNames, gamma)
L = LeakyReluLayer('leaky1', 1, {'in'}, 1, {'out'}, 0.1);

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

% Test with approx-star (exact-star may produce many sets)
[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is, 'approx-star', 50);
assert(passed, 'LeakyReLU approx-star soundness test failed: %s', msg);

%% Test 2: Different gamma values
gammas = [0.01, 0.1, 0.2];
for gamma_val = gammas
    L2 = LeakyReluLayer('leaky2', 1, {'in'}, 1, {'out'}, gamma_val);

    V2 = zeros(2, 2, 1, 2);
    V2(:,:,1,1) = [-0.3 0.5; 0.2 -0.4];
    V2(1,1,1,2) = 0.2;

    input_is2 = ImageStar(V2, [1; -1], [1; 1], -1, 1);

    [passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L2, input_is2, 'approx-star', 30);
    assert(passed, 'LeakyReLU (gamma=%.2f) soundness test failed: %s', gamma_val, msg);
end

%% Test 3: All-positive input (LeakyReLU identity)
L3 = LeakyReluLayer('leaky3', 1, {'in'}, 1, {'out'}, 0.1);

V3 = zeros(3, 3, 1, 2);
V3(:,:,1,1) = ones(3, 3);  % all positive center
V3(2,2,1,2) = 0.1;         % small perturbation, still all positive

input_is3 = ImageStar(V3, [1; -1], [0.5; 0.5], -0.5, 0.5);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L3, input_is3, 'approx-star', 30);
assert(passed, 'LeakyReLU all-positive soundness test failed: %s', msg);

%% Test 4: All-negative input (scaled by gamma)
L4 = LeakyReluLayer('leaky4', 1, {'in'}, 1, {'out'}, 0.2);

V4 = zeros(3, 3, 1, 2);
V4(:,:,1,1) = -ones(3, 3);  % all negative center
V4(2,2,1,2) = 0.1;          % small perturbation, still all negative

input_is4 = ImageStar(V4, [1; -1], [0.5; 0.5], -0.5, 0.5);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L4, input_is4, 'approx-star', 30);
assert(passed, 'LeakyReLU all-negative soundness test failed: %s', msg);

%% Test 5: Multi-channel LeakyReLU
rng(42);
L5 = LeakyReluLayer('leaky5', 1, {'in'}, 1, {'out'}, 0.1);

V5 = zeros(2, 2, 3, 3);
V5(:,:,:,1) = randn(2, 2, 3) * 0.5;  % random center
V5(1,1,1,2) = 0.2;
V5(2,2,2,3) = 0.2;

C5 = [eye(2); -eye(2)];
d5 = [1; 1; 1; 1];

input_is5 = ImageStar(V5, C5, d5, [-1; -1], [1; 1]);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L5, input_is5, 'approx-star', 50);
assert(passed, 'Multi-channel LeakyReLU soundness test failed: %s', msg);
