% test_soundness_AveragePooling2DLayer
% Soundness tests for AveragePooling2DLayer
% Verifies that computed reachable sets contain all actual outputs
% To run: results = runtests('test_soundness_AveragePooling2DLayer')

%% Test 1: Basic AvgPool soundness with exact-star
% Create AvgPooling layer: 2x2 pool, stride 2
L = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);

% Create 4x4 ImageStar input with 2 predicate variables
V = zeros(4, 4, 1, 3);
V(:,:,1,1) = [1 2 3 4; 2 3 4 5; 3 4 5 6; 4 5 6 7] / 10;  % center
V(1,1,1,2) = 0.1;  % basis 1
V(2,2,1,3) = 0.1;  % basis 2

C = [eye(2); -eye(2)];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

input_is = ImageStar(V, C, d, pred_lb, pred_ub);

% Verify soundness
[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L, input_is, 'exact-star', 50);
assert(passed, 'AvgPool exact-star soundness test failed: %s', msg);

%% Test 2: AvgPool with different pool sizes
rng(42);
L2 = AveragePooling2DLayer([3 3], [3 3], [0 0 0 0]);

V2 = zeros(6, 6, 1, 2);
V2(:,:,1,1) = rand(6, 6);
V2(2,2,1,2) = 0.2;

input_is2 = ImageStar(V2, [1; -1], [1; 1], -1, 1);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L2, input_is2, 'exact-star', 30);
assert(passed, 'AvgPool 3x3 soundness test failed: %s', msg);

%% Test 3: Multi-channel AvgPool
rng(42);
L3 = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);

V3 = zeros(4, 4, 3, 3);
V3(:,:,:,1) = rand(4, 4, 3);
V3(1,1,1,2) = 0.1;
V3(3,3,2,3) = 0.1;

input_is3 = ImageStar(V3, [eye(2); -eye(2)], ones(4,1), [-1;-1], [1;1]);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L3, input_is3, 'exact-star', 50);
assert(passed, 'Multi-channel AvgPool soundness test failed: %s', msg);

%% Test 4: AvgPool with padding
L4 = AveragePooling2DLayer([2 2], [2 2], [1 1 1 1]);

V4 = zeros(3, 3, 1, 2);
V4(:,:,1,1) = rand(3, 3);
V4(2,2,1,2) = 0.2;

input_is4 = ImageStar(V4, [1; -1], [0.5; 0.5], -0.5, 0.5);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L4, input_is4, 'exact-star', 30);
assert(passed, 'AvgPool with padding soundness test failed: %s', msg);

%% Test 5: Corner case sampling
L5 = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);

V5 = zeros(4, 4, 1, 3);
V5(:,:,1,1) = [0.1 0.2 0.3 0.4; 0.2 0.3 0.4 0.5; 0.3 0.4 0.5 0.6; 0.4 0.5 0.6 0.7];
V5(1,1,1,2) = 0.15;
V5(4,4,1,3) = 0.15;

C5 = [eye(2); -eye(2)];
d5 = [0.5; 0.5; 0.5; 0.5];
pred_lb5 = [-0.5; -0.5];
pred_ub5 = [0.5; 0.5];

input_is5 = ImageStar(V5, C5, d5, pred_lb5, pred_ub5);
output_is5 = L5.reach(input_is5, 'exact-star');

% Test all corners
corners = [-0.5 -0.5; -0.5 0.5; 0.5 -0.5; 0.5 0.5];
for i = 1:4
    alpha = corners(i, :)';
    input_concrete = V5(:,:,:,1) + alpha(1)*V5(:,:,:,2) + alpha(2)*V5(:,:,:,3);
    output_concrete = L5.evaluate(input_concrete);

    contained = soundness_test_utils.verify_imagestar_containment(output_is5, output_concrete);
    assert(contained, 'Corner %d not contained in AvgPool output', i);
end

%% Test 6: AvgPool preserves linearity (exact reachability)
% Since AvgPool is linear, the output ImageStar should exactly represent the set
rng(42);
L6 = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);

V6 = zeros(4, 4, 1, 2);
V6(:,:,1,1) = rand(4, 4);
V6(:,:,1,2) = rand(4, 4) * 0.1;

input_is6 = ImageStar(V6, [1; -1], [1; 1], -1, 1);
output_is6 = L6.reach(input_is6, 'exact-star');

% For 100 random samples, all should be contained
n_samples = 100;
for i = 1:n_samples
    alpha = -1 + 2*rand();
    input_concrete = V6(:,:,:,1) + alpha * V6(:,:,:,2);
    output_concrete = L6.evaluate(input_concrete);

    contained = soundness_test_utils.verify_imagestar_containment(output_is6, output_concrete);
    assert(contained, 'Sample %d not contained in linear AvgPool output', i);
end
