% test_soundness_ReshapeToConcatenationLayer
% Tests for ReshapeToConcatenation layer
% To run: results = runtests('test_soundness_ReshapeToConcatenationLayer')

%% Test 1: ReshapeToConcatenationLayer constructor
rng(42);

name = 'reshape_concat_test';
numInputs = 2;
numOutputs = 1;
inputNames = {'in1', 'in2'};
outputNames = {'out'};
targetShapes = {[4 4 1], [4 4 1]};
concatDim = 3;

L = ReshapeToConcatenationLayer(name, numInputs, numOutputs, inputNames, outputNames, targetShapes, concatDim);

assert(strcmp(L.Name, name), 'Name should be set');
assert(L.NumInputs == numInputs, 'NumInputs should match');
assert(L.NumOutputs == numOutputs, 'NumOutputs should match');
assert(L.ConcatDim == concatDim, 'ConcatDim should match');

%% Test 2: ReshapeToConcatenationLayer with 3 inputs
rng(42);

name = 'reshape_concat_3';
numInputs = 3;
numOutputs = 1;
inputNames = {'a', 'b', 'c'};
outputNames = {'output0'};
targetShapes = {[2 2 1], [2 2 1], [2 2 1]};
concatDim = 3;

L2 = ReshapeToConcatenationLayer(name, numInputs, numOutputs, inputNames, outputNames, targetShapes, concatDim);

assert(L2.NumInputs == 3, 'Should have 3 inputs');
assert(length(L2.InputNames) == 3, 'Should have 3 input names');
assert(length(L2.TargetShapes) == 3, 'Should have 3 target shapes');

%% Test 3: ReshapeToConcatenationLayer properties
rng(42);

L3 = ReshapeToConcatenationLayer('test', 2, 1, {'x', 'y'}, {'z'}, {[3 3 1], [3 3 1]}, 3);

assert(strcmp(L3.InputNames{1}, 'x'), 'First input name should be x');
assert(strcmp(L3.InputNames{2}, 'y'), 'Second input name should be y');
assert(strcmp(L3.OutputNames{1}, 'z'), 'Output name should be z');

