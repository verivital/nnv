% test_soundness_FlattenLayer
% Soundness tests for FlattenLayer
% Verifies that computed reachable sets contain all actual outputs
% Note: FlattenLayer.reach returns ImageStar with shape [1, 1, N, num_pred+1]
% To run: results = runtests('test_soundness_FlattenLayer')

%% Test 1: Basic Flatten soundness with ImageStar
L = FlattenLayer('flatten_test');
L.Type = 'nnet.cnn.layer.FlattenLayer';  % Required for evaluate

% Create 3x3x2 ImageStar (flattens to 18-element vector stored as [1,1,18])
V = zeros(3, 3, 2, 3);
V(:,:,1,1) = [1 2 3; 4 5 6; 7 8 9] / 10;
V(:,:,2,1) = [10 11 12; 13 14 15; 16 17 18] / 100;
V(1,1,1,2) = 0.1;  % basis 1
V(2,2,2,3) = 0.05; % basis 2

C = [eye(2); -eye(2)];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

input_is = ImageStar(V, C, d, pred_lb, pred_ub);

% Compute reachable set (Flatten returns ImageStar with shape [1,1,N,num_pred+1])
output_is = L.reach(input_is, 'exact-star');

% Test random samples - use ImageStar containment
n_samples = 50;
for i = 1:n_samples
    alpha = pred_lb + (pred_ub - pred_lb) .* rand(2, 1);
    if any(C * alpha > d)
        continue;
    end

    input_concrete = V(:,:,:,1) + alpha(1)*V(:,:,:,2) + alpha(2)*V(:,:,:,3);
    output_concrete = L.evaluate(input_concrete);

    % Reshape output to match ImageStar format [1,1,N]
    output_reshaped = reshape(output_concrete, [1, 1, numel(output_concrete)]);
    contained = soundness_test_utils.verify_imagestar_containment(output_is, output_reshaped);
    assert(contained, 'Flatten soundness violation at sample %d', i);
end

%% Test 2: Single channel Flatten
L2 = FlattenLayer('flatten2');
L2.Type = 'nnet.cnn.layer.FlattenLayer';

V2 = zeros(2, 2, 1, 2);
V2(:,:,1,1) = [0.5 0.6; 0.7 0.8];
V2(1,1,1,2) = 0.2;

input_is2 = ImageStar(V2, [1; -1], [1; 1], -1, 1);
output_is2 = L2.reach(input_is2, 'exact-star');

% Test corners
corners = [-1, 0, 1];
for i = 1:length(corners)
    alpha = corners(i);
    input_concrete = V2(:,:,:,1) + alpha * V2(:,:,:,2);
    output_concrete = L2.evaluate(input_concrete);

    output_reshaped = reshape(output_concrete, [1, 1, numel(output_concrete)]);
    contained = soundness_test_utils.verify_imagestar_containment(output_is2, output_reshaped);
    assert(contained, 'Single-channel Flatten corner %d not contained', i);
end

%% Test 3: Larger multi-channel image
rng(42);
L3 = FlattenLayer('flatten3');
L3.Type = 'nnet.cnn.layer.FlattenLayer';

V3 = zeros(4, 4, 3, 3);
V3(:,:,:,1) = rand(4, 4, 3);
V3(1,1,1,2) = 0.1;
V3(3,3,2,3) = 0.1;

input_is3 = ImageStar(V3, [eye(2); -eye(2)], ones(4,1), [-1;-1], [1;1]);
output_is3 = L3.reach(input_is3, 'exact-star');

% Test random samples
for i = 1:30
    alpha = -1 + 2*rand(2, 1);
    input_concrete = V3(:,:,:,1) + alpha(1)*V3(:,:,:,2) + alpha(2)*V3(:,:,:,3);
    output_concrete = L3.evaluate(input_concrete);

    output_reshaped = reshape(output_concrete, [1, 1, numel(output_concrete)]);
    contained = soundness_test_utils.verify_imagestar_containment(output_is3, output_reshaped);
    assert(contained, 'Multi-channel Flatten sample %d not contained', i);
end

%% Test 4: Flatten preserves set membership (linearity check)
L4 = FlattenLayer('flatten4');
L4.Type = 'nnet.cnn.layer.FlattenLayer';

V4 = zeros(2, 2, 1, 2);
V4(:,:,1,1) = [1 2; 3 4];
V4(:,:,1,2) = [0.1 0.2; 0.3 0.4];

input_is4 = ImageStar(V4, [1; -1], [1; 1], -1, 1);
output_is4 = L4.reach(input_is4, 'exact-star');

% Verify output dimensions: should be [1, 1, 4, 2] (ImageStar with 4 flattened elements)
expected_elements = 2 * 2 * 1;  % H * W * C = 4
assert(size(output_is4.V, 3) == expected_elements, ...
    'Flattened output should have %d elements, got %d', expected_elements, size(output_is4.V, 3));
