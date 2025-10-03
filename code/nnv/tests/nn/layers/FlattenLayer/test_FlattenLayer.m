% Test FlattenLayer functionality
% To run: results = runtests('test_FlattenLayer')

%% Test 1: FlattenLayer constructor - default
L = FlattenLayer();
assert(strcmp(L.Name, 'FlattenLayer'));
assert(L.NumInputs == 1);
assert(L.NumOutputs == 1);

%% Test 2: FlattenLayer constructor - with name
L = FlattenLayer('flatten1');
assert(strcmp(L.Name, 'flatten1'));
assert(L.NumInputs == 1);

%% Test 3: FlattenLayer constructor - with all parameters
L = FlattenLayer('flatten2', 1, 1, {'input'}, {'output'});
assert(strcmp(L.Name, 'flatten2'));
assert(L.NumInputs == 1);
assert(L.NumOutputs == 1);
assert(strcmp(L.InputNames{1}, 'input'));
assert(strcmp(L.OutputNames{1}, 'output'));

%% Test 4: FlattenLayer evaluate - 2D input (MATLAB style)
L = FlattenLayer();
L.Type = 'nnet.cnn.layer.FlattenLayer';

% Create 2D input
input = [1 2 3; 4 5 6];
output = L.evaluate(input);

% Should flatten to [1 1 6]
expected_size = [1 6];
assert(all(size(output) == expected_size), 'FlattenLayer 2D evaluate failed');
assert(numel(output) == 6, 'FlattenLayer should preserve number of elements');

%% Test 5: FlattenLayer evaluate - 3D input (MATLAB style)
L = FlattenLayer();
L.Type = 'nnet.cnn.layer.FlattenLayer';

% Create 3D image
IM(:,:,1) = [1 2; 3 4];
IM(:,:,2) = [5 6; 7 8];
output = L.evaluate(IM);

% Should flatten to [1 1 8]
assert(size(output, 1) == 1);
assert(size(output, 2) == 1);
assert(size(output, 3) == 8);
assert(numel(output) == 8, 'FlattenLayer should preserve number of elements');

%% Test 6: FlattenLayer evaluate - 3D input (Keras C-style)
L = FlattenLayer();
L.Type = 'nnet.keras.layer.FlattenCStyleLayer';

% Create 3D image
IM(:,:,1) = [1 2; 3 4];
IM(:,:,2) = [5 6; 7 8];
output = L.evaluate(IM);

% Keras C-style flatten uses different permutation
assert(size(output, 1) == 1);
assert(size(output, 2) == 1);
assert(size(output, 3) == 8);

%% Test 7: FlattenLayer evaluate - 4D input (ONNX style)
L = FlattenLayer();
L.Type = 'nnet.onnx.layer.FlattenLayer';

% Create 4D input
input = rand(2, 2, 2, 2);
output = L.evaluate(input);

% Should flatten to [1 1 1 16]
assert(size(output, 1) == 1);
assert(size(output, 2) == 1);
assert(size(output, 3) == 16);
assert(size(output, 4) == 1);

%% Test 8: FlattenLayer evaluateSequence
L = FlattenLayer();

% Create sequence input
seq_input = rand(10, 5, 3);
output_seq = L.evaluateSequence(seq_input);

% Should squeeze dimensions
assert(ndims(output_seq) <= ndims(seq_input), 'evaluateSequence should reduce dimensions');

%% Test 9: FlattenLayer reach with ImageStar
L = FlattenLayer();
L.Type = 'nnet.cnn.layer.FlattenLayer';

% Create ImageStar input
IM(:,:,1) = [1 1; 0 1];
IM(:,:,2) = [0 1; 1 0];

LB(:,:,1) = [-0.1 -0.1; 0 0];
LB(:,:,2) = [-0.1 -0.1; 0 0];

UB(:,:,1) = [0.1 0.1; 0 0];
UB(:,:,2) = [0.1 0.1; 0 0];

image_star = ImageStar(IM, LB, UB);
output_star = L.reach(image_star, 'approx-star');

assert(isa(output_star, 'ImageStar'), 'reach should return ImageStar');
% Flattened output should have correct dimensions
assert(output_star.height == 1);
assert(output_star.width == 1);

%% Test 11: FlattenLayer toGPU
L = FlattenLayer();
L_gpu = L.toGPU();
assert(isa(L_gpu, 'FlattenLayer'));

%% Test 12: FlattenLayer changeParamsPrecision
L = FlattenLayer();
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'FlattenLayer'));
