% Test SequenceInputLayer functionality
% To run: results = runtests('test_SequenceInputLayer')

%% Test 1: SequenceInputLayer constructor - default
L = SequenceInputLayer();
assert(strcmp(L.Name, 'SequenceInputLayer'));
assert(strcmp(L.Normalization, 'none'));
assert(L.MinLength == 1);

%% Test 2: SequenceInputLayer constructor - with Name and InputSize
L = SequenceInputLayer('Name', 'seq1', 'InputSize', 10);
assert(strcmp(L.Name, 'seq1'));
assert(L.InputSize == 10);

%% Test 3: SequenceInputLayer constructor - with multiple parameters
L = SequenceInputLayer('Name', 'seq2', 'InputSize', 20, 'MinLength', 5);
assert(strcmp(L.Name, 'seq2'));
assert(L.InputSize == 20);
assert(L.MinLength == 5);

%% Test 4: SequenceInputLayer constructor - with normalization
L = SequenceInputLayer('Name', 'seq3', 'Normalization', 'zerocenter', 'Mean', [0.5]);
assert(strcmp(L.Normalization, 'zerocenter'));
assert(L.Mean == 0.5);

%% Test 5: SequenceInputLayer evaluateSequence - no normalization
L = SequenceInputLayer('InputSize', 5, 'Normalization', 'none');

% Test with sequence input (features x time steps)
input = [1 2 3; 4 5 6; 7 8 9];  % 3 features x 3 time steps
output = L.evaluateSequence(input);

% Should return input unchanged (just converted to double)
assert(isequal(double(output), double(input)), 'No normalization should return input unchanged');

%% Test 6: SequenceInputLayer evaluateSequence - with normalization
L = SequenceInputLayer('InputSize', 3, 'Normalization', 'zerocenter', 'Mean', [1; 2; 3]);

% Test with sequence input
input = [2 3 4; 4 5 6; 6 7 8];  % 3 features x 3 time steps
output = L.evaluateSequence(input);

% Should subtract mean from each feature
expected = [1 2 3; 2 3 4; 3 4 5];
assert(all(abs(output(:) - expected(:)) < 1e-10), 'Zerocenter normalization failed');

%% Test 7: SequenceInputLayer reach - no normalization
L = SequenceInputLayer('InputSize', 2, 'Normalization', 'none');

% Create simple Star input
lb = [0 0; 0 0];
ub = [1 1; 1 1];
I_star = ImageStar(lb,ub);

output = L.reachSequence(I_star, 'exact-star');

%% Test 8: SequenceInputLayer toGPU
L = SequenceInputLayer('Name', 'seq_gpu', 'Mean', [1; 2; 3], 'StandardDeviation', [0.5; 0.5; 0.5]);

if gpuDeviceCount > 0
    L_gpu = L.toGPU();
    assert(isa(L_gpu.Mean, 'gpuArray'), 'Mean should be gpuArray after toGPU');
end

%% Test 9: SequenceInputLayer changeParamsPrecision - single
L = SequenceInputLayer('Mean', double([1; 2; 3]), 'StandardDeviation', double([0.5; 0.5; 0.5]));

L_single = L.changeParamsPrecision('single');
assert(isa(L_single.Mean, 'single'), 'Mean should be single precision');
assert(isa(L_single.StandardDeviation, 'single'), 'StandardDeviation should be single precision');

%% Test 10: SequenceInputLayer changeParamsPrecision - double
L = SequenceInputLayer('Mean', single([1; 2; 3]));

L_double = L.changeParamsPrecision('double');
assert(isa(L_double.Mean, 'double'), 'Mean should be double precision');
