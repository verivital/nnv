% test_soundness_ImageInputLayer
% Tests for Image Input layer (normalizes input images)
% To run: results = runtests('test_soundness_ImageInputLayer')

%% Test 1: ImageInputLayer constructor with InputSize
rng(42);

L = ImageInputLayer([28 28 1]);  % MNIST-like

assert(isequal(L.InputSize, [28 28 1]), 'InputSize should match');

%% Test 2: ImageInputLayer constructor with normalization
rng(42);

L2 = ImageInputLayer([32 32 3], 'none');

assert(isequal(L2.InputSize, [32 32 3]), 'InputSize should be [32 32 3]');
assert(strcmp(L2.Normalization, 'none'), 'Normalization should be none');

%% Test 3: ImageInputLayer evaluate with no normalization
rng(42);

L3 = ImageInputLayer([4 4 1], 'none');

input = rand(4, 4);
output = L3.evaluate(input);

% With no normalization, output should equal input
assert(max(abs(output(:) - input(:))) < 1e-10, 'No normalization should preserve input');

%% Test 4: ImageInputLayer reachability with ImageStar
rng(42);

L4 = ImageInputLayer([3 3 1], 'none');

V = zeros(3, 3, 1, 2);
V(:,:,1,1) = rand(3, 3);
V(:,:,1,2) = rand(3, 3) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L4.reach(input_is, 'approx-star');

assert(~isempty(output_is), 'Should produce output');

%% Test 5: ImageInputLayer soundness
rng(42);

L5 = ImageInputLayer([3 3 1], 'none');

V = zeros(3, 3, 1, 2);
V(:,:,1,1) = rand(3, 3);
V(:,:,1,2) = rand(3, 3) * 0.2;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L5.reach(input_is, 'approx-star');

% Verify soundness
for i = 1:5
    alpha = -1 + 2*rand();
    input_concrete = V(:,:,:,1) + alpha * V(:,:,:,2);
    output_concrete = L5.evaluate(input_concrete);

    contained = soundness_test_utils.verify_imagestar_containment(output_is, output_concrete);
    assert(contained, 'Output should contain sample %d', i);
end

