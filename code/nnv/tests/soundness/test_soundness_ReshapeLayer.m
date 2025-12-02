% test_soundness_ReshapeLayer
% Tests for Reshape layer
% To run: results = runtests('test_soundness_ReshapeLayer')

%% Test 1: ReshapeLayer constructor
rng(42);

L = ReshapeLayer('reshape_test', [4 4 1]);

assert(strcmp(L.Name, 'reshape_test'), 'Name should be set');
assert(isequal(L.targetDim, [4 4 1]), 'targetDim should match');

%% Test 2: ReshapeLayer evaluate
rng(42);

L2 = ReshapeLayer('reshape_eval', [2 8]);

input = rand(4, 4);  % 16 elements
output = L2.evaluate(input);

assert(isequal(size(output), [2 8]), 'Output shape should be [2 8]');
assert(numel(output) == numel(input), 'Number of elements should be preserved');

%% Test 3: ReshapeLayer 3D to 2D
rng(42);

L3 = ReshapeLayer('reshape_3dto2d', [16 1]);

input = rand(4, 4, 1);
output = L3.evaluate(input);

assert(isequal(size(output), [16 1]), 'Output should be column vector');

%% Test 4: ReshapeLayer reachability with ImageStar
rng(42);

L4 = ReshapeLayer('reshape_reach', [9 1]);

V = zeros(3, 3, 1, 2);
V(:,:,1,1) = rand(3, 3);
V(:,:,1,2) = rand(3, 3) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

output_is = L4.reach(input_is, 'approx-star');

assert(~isempty(output_is), 'Should produce output');

%% Test 5: ReshapeLayer soundness
rng(42);

L5 = ReshapeLayer('reshape_sound', [4 1]);

V = zeros(2, 2, 1, 2);
V(:,:,1,1) = rand(2, 2);
V(:,:,1,2) = rand(2, 2) * 0.2;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L5.reach(input_is, 'approx-star');

% Verify soundness
for i = 1:5
    alpha = -1 + 2*rand();
    input_concrete = V(:,:,:,1) + alpha * V(:,:,:,2);
    output_concrete = reshape(input_concrete, [4 1]);

    contained = soundness_test_utils.verify_imagestar_containment(output_is, output_concrete);
    assert(contained, 'Output should contain sample %d', i);
end

