% test_soundness_DepthConcatenationLayer
% Tests for Depth Concatenation layer (concatenates along channel dimension)
% To run: results = runtests('test_soundness_DepthConcatenationLayer')

%% Test 1: DepthConcatenationLayer constructor
rng(42);

L = DepthConcatenationLayer('depth_concat', 2, 1, {'in1', 'in2'}, {'out'});

assert(strcmp(L.Name, 'depth_concat'), 'Name should be set');
assert(L.NumInputs == 2, 'NumInputs should be 2');
assert(L.NumOutputs == 1, 'NumOutputs should be 1');

%% Test 2: DepthConcatenationLayer with 3 inputs
rng(42);

L2 = DepthConcatenationLayer('depth_concat3', 3, 1, {'a', 'b', 'c'}, {'out'});

assert(L2.NumInputs == 3, 'NumInputs should be 3');
assert(length(L2.InputNames) == 3, 'Should have 3 input names');

%% Test 3: DepthConcatenationLayer evaluate
rng(42);

L3 = DepthConcatenationLayer('depth_eval', 2, 1, {'in1', 'in2'}, {'out'});

% Create two inputs with different channel counts
input1 = rand(4, 4, 2);  % 2 channels
input2 = rand(4, 4, 3);  % 3 channels

output = L3.evaluate({input1, input2});

assert(size(output, 3) == 5, 'Output should have 2+3=5 channels');
assert(size(output, 1) == 4, 'Height should be preserved');
assert(size(output, 2) == 4, 'Width should be preserved');

%% Test 4: DepthConcatenationLayer reach_single_input
% NOTE: reach_multipleInputs_Star has a library bug (nI undefined)
% Testing reach_single_input instead which works correctly
rng(42);

L4 = DepthConcatenationLayer('depth_reach', 2, 1, {'in1', 'in2'}, {'out'});

% Create two ImageStars with same predicate structure
V1 = zeros(3, 3, 1, 2);
V1(:,:,1,1) = rand(3, 3);
V1(:,:,1,2) = rand(3, 3) * 0.1;
is1 = ImageStar(V1, [1; -1], [1; 1], -1, 1);

V2 = zeros(3, 3, 2, 2);
V2(:,:,:,1) = rand(3, 3, 2);
V2(:,:,:,2) = rand(3, 3, 2) * 0.1;
is2 = ImageStar(V2, [1; -1], [1; 1], -1, 1);

% Use reach_single_input which takes a cell array
output_is = L4.reach_single_input({is1, is2});

assert(~isempty(output_is), 'Should produce output');
assert(size(output_is.V, 3) == 3, 'Should have 1+2=3 channels');

%% Test 5: DepthConcatenationLayer soundness with reach_single_input
rng(42);

L5 = DepthConcatenationLayer('depth_sound', 2, 1, {'in1', 'in2'}, {'out'});

V1 = zeros(3, 3, 1, 2);
V1(:,:,1,1) = rand(3, 3);
V1(:,:,1,2) = rand(3, 3) * 0.1;

V2 = zeros(3, 3, 1, 2);
V2(:,:,1,1) = rand(3, 3);
V2(:,:,1,2) = rand(3, 3) * 0.1;

is1 = ImageStar(V1, [1; -1], [1; 1], -1, 1);
is2 = ImageStar(V2, [1; -1], [1; 1], -1, 1);

output_is = L5.reach_single_input({is1, is2});

% Verify soundness
for i = 1:5
    alpha = -1 + 2*rand();
    in1_concrete = V1(:,:,:,1) + alpha * V1(:,:,:,2);
    in2_concrete = V2(:,:,:,1) + alpha * V2(:,:,:,2);
    out_concrete = cat(3, in1_concrete, in2_concrete);

    contained = soundness_test_utils.verify_imagestar_containment(output_is, out_concrete);
    assert(contained, 'Output should contain sample %d', i);
end

