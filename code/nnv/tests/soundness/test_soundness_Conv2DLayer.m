% test_soundness_Conv2DLayer
% Soundness tests for Conv2DLayer
% Verifies that computed reachable sets contain all actual outputs
% To run: results = runtests('test_soundness_Conv2DLayer')

%% Test 1: Basic Conv2D soundness with ImageStar
% Create a simple 2x2 filter, 1 input channel, 1 output channel
W = ones(2, 2, 1, 1);  % [H, W, C_in, C_out]
b = 0;
L = Conv2DLayer(W, b);

% Create small ImageStar input (4x4x1 with 2 predicate variables)
V = zeros(4, 4, 1, 3);
V(:,:,1,1) = [1 2 3 4; 2 3 4 5; 3 4 5 6; 4 5 6 7] / 10;  % center
V(1,1,1,2) = 0.1;  % basis 1 affects pixel (1,1)
V(2,2,1,3) = 0.1;  % basis 2 affects pixel (2,2)

C = [1 0; -1 0; 0 1; 0 -1];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

input_is = ImageStar(V, C, d, pred_lb, pred_ub);

% Verify soundness
[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is, 'exact-star', 50);
assert(passed, 'Conv2D soundness test failed: %s', msg);

%% Test 2: Multi-channel Conv2D (single output channel)
% 3x3 filter, 2 input channels, 1 output channel
% Use fixed seed for reproducibility
rng(42);
W2 = randn(3, 3, 2, 1);  % Single output channel to avoid potential bugs
b2 = 0.1;
L2 = Conv2DLayer(W2, b2);

% Create ImageStar with 2 input channels
V2 = zeros(5, 5, 2, 3);
V2(:,:,:,1) = rand(5, 5, 2);  % center
V2(1,1,1,2) = 0.05;  % basis 1
V2(3,3,2,3) = 0.05;  % basis 2

C2 = [eye(2); -eye(2)];
d2 = 0.5 * ones(4, 1);
pred_lb2 = [-0.5; -0.5];
pred_ub2 = [0.5; 0.5];

input_is2 = ImageStar(V2, C2, d2, pred_lb2, pred_ub2);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L2, input_is2, 'exact-star', 50);
assert(passed, 'Multi-channel Conv2D soundness test failed: %s', msg);

%% Test 3: Conv2D with stride (set properties directly)
rng(42);
W3 = randn(2, 2, 1, 1);  % [H, W, C_in, C_out]
b3 = 0;
L3 = Conv2DLayer(W3, b3);
L3.Stride = [2 2];  % Set stride directly

V3 = zeros(6, 6, 1, 2);
V3(:,:,1,1) = rand(6, 6);
V3(2,2,1,2) = 0.1;

C3 = [1; -1];
d3 = [1; 1];
pred_lb3 = -1;
pred_ub3 = 1;

input_is3 = ImageStar(V3, C3, d3, pred_lb3, pred_ub3);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L3, input_is3, 'exact-star', 50);
assert(passed, 'Conv2D with stride soundness test failed: %s', msg);

%% Test 4: Conv2D with padding (set properties directly)
rng(42);
W4 = randn(3, 3, 1, 1);
b4 = 0.1;
L4 = Conv2DLayer(W4, b4);
L4.PaddingSize = [1 1 1 1];  % Set padding directly

V4 = zeros(4, 4, 1, 2);
V4(:,:,1,1) = rand(4, 4);
V4(2,3,1,2) = 0.1;

C4 = [1; -1];
d4 = [0.5; 0.5];
pred_lb4 = -0.5;
pred_ub4 = 0.5;

input_is4 = ImageStar(V4, C4, d4, pred_lb4, pred_ub4);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L4, input_is4, 'exact-star', 50);
assert(passed, 'Conv2D with padding soundness test failed: %s', msg);
