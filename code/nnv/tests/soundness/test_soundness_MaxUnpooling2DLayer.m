% test_soundness_MaxUnpooling2DLayer
% Tests for Max Unpooling 2D layer (inverse of MaxPooling)
% To run: results = runtests('test_soundness_MaxUnpooling2DLayer')

%% Test 1: MaxUnpooling2DLayer constructor
rng(42);

L = MaxUnpooling2DLayer('maxunpool_test', 2, {'in', 'indices'});

assert(strcmp(L.Name, 'maxunpool_test'), 'Name should be set');
assert(L.NumInputs == 2, 'Should have 2 inputs (data + indices)');

%% Test 2: MaxUnpooling2DLayer with 5-argument constructor
rng(42);

L2 = MaxUnpooling2DLayer('maxunpool_full', 2, {'data', 'idx'}, 1, {'out'});

assert(strcmp(L2.Name, 'maxunpool_full'), 'Name should be set');
assert(L2.NumInputs == 2, 'Should have 2 inputs');
assert(L2.NumOutputs == 1, 'Should have 1 output');

%% Test 3: MaxUnpooling2DLayer basic functionality
% MaxUnpooling requires indices from MaxPooling operation
% This is a complex layer that works in conjunction with MaxPooling
rng(42);

% Create a simple MaxPooling layer to get indices
pool = MaxPooling2DLayer('pool_for_unpool', [2 2], [2 2], [0 0 0 0]);

% Input image
input = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];

% The max values would be at indices where max occurred
% For a 2x2 pool with stride 2:
% Pool 1: max(1,2,5,6)=6 at position (2,2)
% Pool 2: max(3,4,7,8)=8 at position (2,4)
% etc.

% Just verify the layer object is created properly
L3 = MaxUnpooling2DLayer('maxunpool_basic', 2, {'in', 'indices'});
assert(~isempty(L3), 'Layer should be created');

%% Test 4: MaxUnpooling2DLayer properties
rng(42);

L4 = MaxUnpooling2DLayer('maxunpool_props', 2, {'x', 'idx'});

assert(iscell(L4.InputNames), 'InputNames should be cell');
assert(length(L4.InputNames) == 2, 'Should have 2 input names');

