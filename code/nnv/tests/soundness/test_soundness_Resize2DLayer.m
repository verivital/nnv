% test_soundness_Resize2DLayer
% Tests for 2D Resize layer
% To run: results = runtests('test_soundness_Resize2DLayer')

%% Test 1: Resize2DLayer constructor with OutputSize
rng(42);

L = Resize2DLayer('resize_test', 1, 1, {'in1'}, {'out'}, [8 8], [], 'nearest', false, 'half-pixel', 'round');

assert(strcmp(L.Name, 'resize_test'), 'Name should be set');
assert(isequal(L.OutputSize, [8 8]), 'OutputSize should be [8 8]');
assert(strcmp(L.Method, 'nearest'), 'Method should be nearest');

%% Test 2: Resize2DLayer constructor with Scale
rng(42);

L2 = Resize2DLayer('resize_scale', 1, 1, {'in1'}, {'out'}, [], [2 2], 'bilinear', false, 'half-pixel', 'round');

assert(isequal(L2.Scale, [2 2]), 'Scale should be [2 2]');
assert(strcmp(L2.Method, 'bilinear'), 'Method should be bilinear');

%% Test 3: Resize2DLayer reachability with ImageStar
rng(42);

L3 = Resize2DLayer('resize_reach', 1, 1, {'in1'}, {'out'}, [6 6], [], 'nearest', false, 'half-pixel', 'round');

% Create small ImageStar
V = zeros(4, 4, 1, 2);
V(:,:,1,1) = rand(4, 4);
V(:,:,1,2) = rand(4, 4) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

output_is = L3.reach_single_input(input_is);

assert(~isempty(output_is), 'Should produce output');
assert(size(output_is.V, 1) == 6, 'Output height should be 6');
assert(size(output_is.V, 2) == 6, 'Output width should be 6');

%% Test 4: Resize2DLayer with scale factor
rng(42);

L4 = Resize2DLayer('resize_scale2', 1, 1, {'in1'}, {'out'}, [], [0.5 0.5], 'nearest', false, 'half-pixel', 'round');

V = zeros(8, 8, 1, 2);
V(:,:,1,1) = rand(8, 8);
V(:,:,1,2) = rand(8, 8) * 0.1;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);

output_is = L4.reach_single_input(input_is);

assert(~isempty(output_is), 'Should produce output');
assert(size(output_is.V, 1) == 4, 'Output height should be 4 (half)');
assert(size(output_is.V, 2) == 4, 'Output width should be 4 (half)');

%% Test 5: Resize2DLayer soundness check
rng(42);

L5 = Resize2DLayer('resize_sound', 1, 1, {'in1'}, {'out'}, [6 6], [], 'bilinear', false, 'half-pixel', 'round');

V = zeros(4, 4, 1, 2);
V(:,:,1,1) = rand(4, 4);
V(:,:,1,2) = rand(4, 4) * 0.2;

input_is = ImageStar(V, [1; -1], [1; 1], -1, 1);
output_is = L5.reach_single_input(input_is);

% Verify soundness
for i = 1:5
    alpha = -1 + 2*rand();
    input_concrete = V(:,:,:,1) + alpha * V(:,:,:,2);
    output_concrete = imresize(input_concrete, [6 6], 'bilinear');

    contained = soundness_test_utils.verify_imagestar_containment(output_is, output_concrete);
    assert(contained, 'Output should contain concrete sample %d', i);
end

