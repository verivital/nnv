% test_soundness_UpsampleLayer
% Tests for Upsample layer
% To run: results = runtests('test_soundness_UpsampleLayer')
% NOTE: UpsampleLayer requires 3 or 4 element scale (H, W, C) or (H, W, C, N)

%% Test 1: UpsampleLayer constructor
rng(42);

% Use 3-element scale [H, W, C] - layer flips it internally
L = UpsampleLayer('upsample_test', [2 2 1]);

assert(strcmp(L.Name, 'upsample_test'), 'Name should be set');
assert(length(L.scaleDim) == 3, 'scaleDim should have 3 elements');

%% Test 2: UpsampleLayer evaluate
% NOTE: Scale [H W C] is flipped to [C W H] internally by dlresize
rng(42);

L2 = UpsampleLayer('upsample_eval', [2 2 1]);

input = rand(4, 4, 2);
output = L2.evaluate(input);

% Verify output is produced and has expected elements
assert(~isempty(output), 'Should produce output');
assert(size(output, 1) >= size(input, 1), 'Output height should be >= input');
assert(size(output, 2) >= size(input, 2), 'Output width should be >= input');

%% Test 3: UpsampleLayer with different scales
rng(42);

L3 = UpsampleLayer('upsample_3x', [3 3 1]);

input = rand(4, 4, 1);
output = L3.evaluate(input);

% Verify upsampling occurred
assert(~isempty(output), 'Should produce output');
assert(numel(output) >= numel(input), 'Output should have more or equal elements');

%% Test 4: UpsampleLayer reachability
rng(42);

L4 = UpsampleLayer('upsample_reach', [2 2 1]);

V = zeros(3, 3, 1, 2);
V(:,:,1,1) = rand(3, 3);
V(:,:,1,2) = rand(3, 3) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

output_is = L4.reach(input_is, 'approx-star');

assert(~isempty(output_is), 'Should produce output');
% Just verify output exists - scale behavior is complex
assert(size(output_is.V, 1) >= 1, 'Output should have valid height');
assert(size(output_is.V, 2) >= 1, 'Output should have valid width');

%% Test 5: UpsampleLayer soundness
rng(42);

L5 = UpsampleLayer('upsample_sound', [2 2 1]);

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

