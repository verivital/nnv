% test_soundness_BatchNormalizationLayer
% Soundness tests for BatchNormalizationLayer
% Verifies that computed reachable sets contain all actual outputs
% To run: results = runtests('test_soundness_BatchNormalizationLayer')

%% Test 1: BatchNorm evaluate correctness
rng(42);
L = BatchNormalizationLayer('Name', 'bn1', ...
    'NumChannels', 1, ...
    'TrainedMean', 0.5, ...
    'TrainedVariance', 0.1, ...
    'Epsilon', 1e-5, ...
    'Scale', 1.2, ...
    'Offset', 0.1);

% Test that evaluate produces correct normalized values
input = rand(4, 4);
output = L.evaluate(input);

% Manual computation: y = Scale * (x - mean) / sqrt(var + eps) + Offset
expected = 1.2 * (input - 0.5) / sqrt(0.1 + 1e-5) + 0.1;
assert(max(abs(output(:) - expected(:))) < 1e-10, ...
    'BatchNorm evaluate should match manual computation');

%% Test 2: BatchNorm with zero mean and unit variance
rng(42);
L2 = BatchNormalizationLayer('Name', 'bn2', ...
    'NumChannels', 1, ...
    'TrainedMean', 0, ...
    'TrainedVariance', 1, ...
    'Epsilon', 1e-5, ...
    'Scale', 1.0, ...
    'Offset', 0);

V2 = zeros(3, 3, 1, 2);
V2(:,:,1,1) = rand(3, 3);
V2(2,2,1,2) = 0.2;

input_is2 = ImageStar(V2, [1; -1], [1; 1], -1, 1);

[passed, msg] = soundness_test_utils.verify_layer_soundness_imagestar(L2, input_is2, 'exact-star', 30);
assert(passed, 'Zero-mean BatchNorm soundness test failed: %s', msg);

%% Test 3: BatchNorm reach produces correct structure
rng(42);
L3 = BatchNormalizationLayer('Name', 'bn3', ...
    'NumChannels', 1, ...
    'TrainedMean', 0.5, ...
    'TrainedVariance', 0.25, ...
    'Epsilon', 1e-5, ...
    'Scale', 1.5, ...
    'Offset', 0.2);

V3 = zeros(2, 2, 1, 2);
V3(:,:,1,1) = [0.5 0.6; 0.7 0.8];
V3(1,1,1,2) = 0.2;

input_is3 = ImageStar(V3, [1; -1], [0.5; 0.5], -0.5, 0.5);
output_is3 = L3.reach(input_is3, 'exact-star');

% Verify structure
assert(isa(output_is3, 'ImageStar'), 'Output should be ImageStar');
assert(output_is3.numPred == input_is3.numPred, 'BatchNorm should preserve num predicates');

%% Test 4: BatchNorm random sampling
rng(42);
L4 = BatchNormalizationLayer('Name', 'bn4', ...
    'NumChannels', 1, ...
    'TrainedMean', 0.3, ...
    'TrainedVariance', 0.2, ...
    'Epsilon', 1e-5, ...
    'Scale', 2.0, ...
    'Offset', -0.5);

V4 = zeros(3, 3, 1, 2);
V4(:,:,1,1) = rand(3, 3);
V4(:,:,1,2) = rand(3, 3) * 0.1;

input_is4 = ImageStar(V4, [1; -1], [1; 1], -1, 1);
output_is4 = L4.reach(input_is4, 'exact-star');

% Sample random points and verify containment
for i = 1:30
    alpha = -1 + 2*rand();
    input_concrete = V4(:,:,:,1) + alpha * V4(:,:,:,2);
    output_concrete = L4.evaluate(input_concrete);
    output_expected = output_is4.V(:,:,:,1) + alpha * output_is4.V(:,:,:,2);

    % Verify direct computation matches
    assert(max(abs(output_concrete(:) - output_expected(:))) < 1e-10, ...
        'BatchNorm sample %d: output should match ImageStar formula', i);
end
