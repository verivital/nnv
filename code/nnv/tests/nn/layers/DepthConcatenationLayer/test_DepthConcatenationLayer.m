% Test DepthConcatenationLayer functionality
% To run: results = runtests('test_DepthConcatenationLayer')

%% Test 1: DepthConcatenationLayer constructor - default
L = DepthConcatenationLayer();
assert(strcmp(L.Name, 'DepthConcatLayer'));
assert(L.NumInputs == 1);
assert(L.NumOutputs == 1);

%% Test 2: DepthConcatenationLayer constructor - with name
L = DepthConcatenationLayer('depthconcat1');
assert(strcmp(L.Name, 'depthconcat1'));

%% Test 3: DepthConcatenationLayer constructor - with all parameters
L = DepthConcatenationLayer('depthconcat2', 2, 1, {'in1', 'in2'}, {'out'});
assert(strcmp(L.Name, 'depthconcat2'));
assert(L.NumInputs == 2);

%% Test 4: DepthConcatenationLayer evaluate - two 3D inputs
L = DepthConcatenationLayer();

% Create two 3D arrays (H x W x C)
input1(:,:,1) = [1 2; 3 4];
input1(:,:,2) = [5 6; 7 8];

input2(:,:,1) = [9 10; 11 12];

inputs = {input1, input2};
output = L.evaluate(inputs);

% Should concatenate along 3rd dimension (channels)
assert(size(output, 1) == 2, 'Height should be unchanged');
assert(size(output, 2) == 2, 'Width should be unchanged');
assert(size(output, 3) == 3, 'Channels should be sum: 2 + 1 = 3');

%% Test 5: DepthConcatenationLayer evaluate - three inputs
L = DepthConcatenationLayer();

% Create three 2x2x1 arrays
input1(:,:,1) = [1 2; 3 4];
input2(:,:,1) = [5 6; 7 8];
input3(:,:,1) = [9 10; 11 12];

inputs = {input1, input2, input3};
output = L.evaluate(inputs);

% Should have 3 channels
assert(size(output, 3) == 3, 'Should have 3 channels after concatenation');

%% Test 6: DepthConcatenationLayer evaluate - different channel counts
L = DepthConcatenationLayer();

% Create inputs with different channel counts
input1 = rand(4, 4, 2);  % 2 channels
input2 = rand(4, 4, 3);  % 3 channels
input3 = rand(4, 4, 1);  % 1 channel

inputs = {input1, input2, input3};
output = L.evaluate(inputs);

% Total channels should be 2 + 3 + 1 = 6
assert(size(output, 3) == 6, 'Total channels should be 6');

%% Test 9: DepthConcatenationLayer toGPU
L = DepthConcatenationLayer();
L_gpu = L.toGPU();
assert(isa(L_gpu, 'DepthConcatenationLayer'));

%% Test 10: DepthConcatenationLayer changeParamsPrecision
L = DepthConcatenationLayer();
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'DepthConcatenationLayer'));
