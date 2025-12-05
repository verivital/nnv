% test_soundness_AveragePooling3DLayer
% Tests for 3D Average Pooling layer
% To run: results = runtests('test_soundness_AveragePooling3DLayer')

%% Test 1: AveragePooling3DLayer constructor
rng(42);

poolSize = [2 2 2];
stride = [2 2 2];
padding = [0 0 0; 0 0 0];  % 2x3 matrix

L = AveragePooling3DLayer('avgpool3d_test', poolSize, stride, padding);

assert(strcmp(L.Name, 'avgpool3d_test'), 'Name should be set correctly');
assert(isequal(L.PoolSize, poolSize), 'PoolSize should match');
assert(isequal(L.Stride, stride), 'Stride should match');

%% Test 2: AveragePooling3DLayer with different pool sizes
rng(42);

poolSize = [3 3 3];
stride = [1 1 1];
padding = [1 1 1; 1 1 1];  % Symmetric padding

L2 = AveragePooling3DLayer('avgpool3d_padded', poolSize, stride, padding);

assert(isequal(L2.PoolSize, poolSize), 'PoolSize should be 3x3x3');
assert(isequal(L2.PaddingSize, padding), 'Padding should be set correctly');

%% Test 3: AveragePooling3DLayer evaluate
rng(42);

poolSize = [2 2 2];
stride = [2 2 2];
padding = [0 0 0; 0 0 0];

L3 = AveragePooling3DLayer('avgpool3d_eval', poolSize, stride, padding);

% Create 4x4x4 input with 1 channel
input = rand(4, 4, 4, 1);

% Evaluate (check if method exists)
if ismethod(L3, 'evaluate')
    output = L3.evaluate(input);
    % Output should be 2x2x2x1
    assert(size(output, 1) == 2, 'Output height should be 2');
    assert(size(output, 2) == 2, 'Output width should be 2');
    assert(size(output, 3) == 2, 'Output depth should be 2');
end

%% Test 4: AveragePooling3DLayer properties validation
rng(42);

poolSize = [2 2 2];
stride = [1 1 1];
padding = [0 0 0; 0 0 0];

L4 = AveragePooling3DLayer('avgpool3d_props', poolSize, stride, padding);

assert(L4.NumInputs == 1, 'Should have 1 input');
assert(L4.NumOutputs == 1, 'Should have 1 output');

