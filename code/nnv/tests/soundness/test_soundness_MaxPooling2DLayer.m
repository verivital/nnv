% test_soundness_MaxPooling2DLayer
% Soundness tests for MaxPooling2DLayer
% Verifies that computed reachable sets contain all actual outputs
% To run: results = runtests('test_soundness_MaxPooling2DLayer')

%% Test 1: Basic MaxPool soundness with exact-star
% Create MaxPooling layer: 2x2 pool, stride 2
L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

% Create 4x4 ImageStar input with 2 predicate variables
V = zeros(4, 4, 1, 3);
V(:,:,1,1) = [1 2 3 4; 2 3 4 5; 3 4 5 6; 4 5 6 7] / 10;  % center
V(1,1,1,2) = 0.15;  % basis 1 - can change max in top-left pool
V(2,2,1,3) = 0.15;  % basis 2 - can change max in top-left pool

C = [eye(2); -eye(2)];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

input_is = ImageStar(V, C, d, pred_lb, pred_ub);

% Verify soundness with exact-star (may produce multiple sets)
[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is, 'exact-star', 30);
assert(passed, 'MaxPool exact-star soundness test failed: %s', msg);

%% Test 2: MaxPool with approx-star (using stable max input)
% Recreate layer and input with smaller perturbations to ensure stable max
L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

V = zeros(4, 4, 1, 3);
V(:,:,1,1) = [0.1 0.2 0.3 0.4; 0.2 0.9 0.4 0.5; 0.3 0.4 0.5 0.6; 0.4 0.5 0.6 0.8];
V(2,2,1,2) = 0.01;  % small perturbation on max element (0.9)
V(4,4,1,3) = 0.01;  % small perturbation on other max element (0.8)

C = [eye(2); -eye(2)];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

input_is = ImageStar(V, C, d, pred_lb, pred_ub);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is, 'approx-star', 50);
assert(passed, 'MaxPool approx-star soundness test failed: %s', msg);

%% Test 3: Multi-channel MaxPool
L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

V2 = zeros(4, 4, 2, 3);
V2(:,:,1,1) = rand(4, 4);
V2(:,:,2,1) = rand(4, 4);
V2(1,1,1,2) = 0.2;
V2(3,3,2,3) = 0.2;

C2 = [eye(2); -eye(2)];
d2 = [0.5; 0.5; 0.5; 0.5];
pred_lb2 = [-0.5; -0.5];
pred_ub2 = [0.5; 0.5];

input_is2 = ImageStar(V2, C2, d2, pred_lb2, pred_ub2);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is2, 'approx-star', 50);
assert(passed, 'Multi-channel MaxPool soundness test failed: %s', msg);

%% Test 4: MaxPool with clear maximum (no crossing)
% When perturbations don't change which element is max
L = MaxPooling2DLayer([2 2], [2 2], [0 0 0 0]);

V3 = zeros(4, 4, 1, 2);
V3(:,:,1,1) = [0.1 0.2 0.3 0.4; 0.2 0.9 0.4 0.5; 0.3 0.4 0.5 0.6; 0.4 0.5 0.6 0.7];
V3(2,2,1,2) = 0.01;  % small perturbation on max element

C3 = [1; -1];
d3 = [1; 1];
pred_lb3 = -1;
pred_ub3 = 1;

input_is3 = ImageStar(V3, C3, d3, pred_lb3, pred_ub3);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is3, 'exact-star', 30);
assert(passed, 'Clear max MaxPool soundness test failed: %s', msg);

%% Test 5: MaxPool with larger pool size
% Test 3x3 pool on 6x6 input
L2 = MaxPooling2DLayer([3 3], [3 3], [0 0 0 0]);

V4 = zeros(6, 6, 1, 2);
V4(:,:,1,1) = rand(6, 6);
V4(2,3,1,2) = 0.1;

C4 = [1; -1];
d4 = [0.5; 0.5];
pred_lb4 = -0.5;
pred_ub4 = 0.5;

input_is4 = ImageStar(V4, C4, d4, pred_lb4, pred_ub4);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L2, input_is4, 'approx-star', 50);
assert(passed, 'Larger MaxPool soundness test failed: %s', msg);
