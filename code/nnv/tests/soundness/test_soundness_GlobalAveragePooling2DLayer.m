% test_soundness_GlobalAveragePooling2DLayer
% Tests for Global Average Pooling 2D layer
% To run: results = runtests('test_soundness_GlobalAveragePooling2DLayer')

%% Test 1: GlobalAveragePooling2DLayer constructor
rng(42);

L = GlobalAveragePooling2DLayer('gap2d_test');

assert(strcmp(L.Name, 'gap2d_test'), 'Name should be set');

%% Test 2: GlobalAveragePooling2DLayer evaluate
rng(42);

L2 = GlobalAveragePooling2DLayer('gap2d_eval');

% Create input: H x W x C
input = rand(4, 4, 3);

output = L2.evaluate(input);

% Output should be 1x1xC (global average per channel)
assert(size(output, 1) == 1, 'Output height should be 1');
assert(size(output, 2) == 1, 'Output width should be 1');
assert(size(output, 3) == 3, 'Output channels should be preserved');

% Verify values are means
for c = 1:3
    expected_mean = mean(mean(input(:,:,c)));
    assert(abs(output(1,1,c) - expected_mean) < 1e-10, 'Output should be mean of channel %d', c);
end

%% Test 3: GlobalAveragePooling2DLayer with different sizes
rng(42);

L3 = GlobalAveragePooling2DLayer('gap2d_sizes');

sizes = {[8 8 1], [16 16 4], [32 32 8]};
for i = 1:length(sizes)
    sz = sizes{i};
    input = rand(sz(1), sz(2), sz(3));
    output = L3.evaluate(input);

    assert(size(output, 1) == 1, 'Output height should be 1 for size %d', i);
    assert(size(output, 2) == 1, 'Output width should be 1 for size %d', i);
    assert(size(output, 3) == sz(3), 'Channels should be preserved for size %d', i);
end

%% Test 4: GlobalAveragePooling2DLayer reachability
rng(42);

L4 = GlobalAveragePooling2DLayer('gap2d_reach');

V = zeros(4, 4, 2, 2);
V(:,:,:,1) = rand(4, 4, 2);
V(:,:,:,2) = rand(4, 4, 2) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

output_is = L4.reach(input_is, 'approx-star');

assert(~isempty(output_is), 'Should produce output');

%% Test 5: GlobalAveragePooling2DLayer soundness
rng(42);

L5 = GlobalAveragePooling2DLayer('gap2d_sound');

V = zeros(4, 4, 1, 2);
V(:,:,1,1) = rand(4, 4);
V(:,:,1,2) = rand(4, 4) * 0.2;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L5.reach(input_is, 'approx-star');

% Verify soundness
for i = 1:5
    alpha = -1 + 2*rand();
    input_concrete = V(:,:,:,1) + alpha * V(:,:,:,2);
    output_concrete = mean(mean(input_concrete, 1), 2);

    contained = soundness_test_utils.verify_imagestar_containment(output_is, output_concrete);
    assert(contained, 'Output should contain sample %d', i);
end

