% test_soundness_TransposedConv2DLayer
% Tests for 2D Transposed Convolution layer
% To run: results = runtests('test_soundness_TransposedConv2DLayer')

%% Test 1: TransposedConv2DLayer constructor
rng(42);

% Weights: [H x W x NumFilters x NumChannels]
weights = randn(3, 3, 2, 1) * 0.1;
bias = zeros(1, 1, 2);
cropping = [0 0 0 0];
stride = [1 1];

L = TransposedConv2DLayer('tconv2d_test', weights, bias, cropping, stride);

assert(strcmp(L.Name, 'tconv2d_test'), 'Name should be set');
assert(L.NumFilters == 2, 'NumFilters should be 2');

%% Test 2: TransposedConv2DLayer with stride
rng(42);

weights = randn(3, 3, 4, 2) * 0.1;
bias = zeros(1, 1, 4);
cropping = [0 0 0 0];
stride = [2 2];

L2 = TransposedConv2DLayer('tconv2d_stride', weights, bias, cropping, stride);

assert(isequal(L2.Stride, stride), 'Stride should match');
assert(L2.NumFilters == 4, 'NumFilters should be 4');

%% Test 3: TransposedConv2DLayer evaluate
rng(42);

weights = randn(3, 3, 2, 1) * 0.1;
bias = zeros(1, 1, 2);

L3 = TransposedConv2DLayer('tconv2d_eval', weights, bias, [0 0 0 0], [1 1]);

input = rand(4, 4, 1);
output = L3.evaluate(input);

% Transposed conv increases spatial dimensions
assert(size(output, 1) >= size(input, 1), 'Output height should be >= input');
assert(size(output, 2) >= size(input, 2), 'Output width should be >= input');
assert(size(output, 3) == 2, 'Output channels should be NumFilters');

%% Test 4: TransposedConv2DLayer reachability
rng(42);

weights = randn(3, 3, 2, 1) * 0.1;
bias = zeros(1, 1, 2);

L4 = TransposedConv2DLayer('tconv2d_reach', weights, bias, [0 0 0 0], [1 1]);

V = zeros(4, 4, 1, 2);
V(:,:,1,1) = rand(4, 4);
V(:,:,1,2) = rand(4, 4) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

output_is = L4.reach(input_is, 'approx-star');

assert(~isempty(output_is), 'Should produce output');

%% Test 5: TransposedConv2DLayer soundness
rng(42);

weights = randn(2, 2, 2, 1) * 0.1;  % 2 filters to avoid single-filter edge case
bias = zeros(1, 1, 2);  % 3D bias: [1, 1, numFilters]

L5 = TransposedConv2DLayer('tconv2d_sound', weights, bias, [0 0 0 0], [1 1]);

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

