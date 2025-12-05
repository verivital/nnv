% Test FeatureInputLayer functionality
% To run: results = runtests('test_FeatureInputLayer')

%% Test 1: FeatureInputLayer constructor - default
L = FeatureInputLayer();
assert(strcmp(L.Name, 'FeatureInputLayer'));
assert(strcmp(L.Normalization, 'none'));

%% Test 2: FeatureInputLayer constructor - with input size
L = FeatureInputLayer([10 1]);
assert(isequal(L.InputSize, [10 1]));
assert(strcmp(L.Normalization, 'none'));

%% Test 3: FeatureInputLayer constructor - with input size and normalization
L = FeatureInputLayer([10 1], 'zerocenter');
assert(isequal(L.InputSize, [10 1]));
assert(strcmp(L.Normalization, 'zerocenter'));

%% Test 4: FeatureInputLayer constructor - with all parameters
L = FeatureInputLayer('feature1', [5 1], 'zerocenter', 'auto', [0.5], [1.0], [0], [1]);
assert(strcmp(L.Name, 'feature1'));
assert(isequal(L.InputSize, [5 1]));
assert(strcmp(L.Normalization, 'zerocenter'));
assert(L.Mean == 0.5);

%% Test 5: FeatureInputLayer evaluate - no normalization
L = FeatureInputLayer([4 1], 'none');

% Test with simple input
input = [1; 2; 3; 4];
output = L.evaluate(input);

% Should return input unchanged (just converted to double)
assert(isequal(double(output), double(input)), 'No normalization should return input unchanged');

%% Test 6: FeatureInputLayer evaluate - zerocenter normalization
L = FeatureInputLayer([4 1], 'zerocenter');
L.Mean = [2; 2; 2; 2];

% Test with simple input
input = [3; 4; 5; 6];
output = L.evaluate(input);

% Should subtract mean
expected = [1; 2; 3; 4];
assert(all(abs(output - expected) < 1e-10), 'Zerocenter should subtract mean');

%% Test 7: FeatureInputLayer evaluate - zerocenter with single mean value
L = FeatureInputLayer([3 1], 'zerocenter');
L.Mean = 5;  % scalar mean

% Test with input
input = [5; 10; 15];
output = L.evaluate(input);

% Should subtract scalar mean from all values
expected = [0; 5; 10];
assert(all(abs(output - expected) < 1e-10), 'Zerocenter with scalar mean failed');

%% Test 8: FeatureInputLayer reach with Star - no normalization
L = FeatureInputLayer([2 1], 'none');

% Create Star input
lb = [0; 0];
ub = [1; 1];
B = Box(lb, ub);
I_star = B.toStar;

output = L.reach(I_star, 'exact-star');
assert(isa(output, 'Star'), 'reach should return Star');

%% Test 9: FeatureInputLayer reach with Star - zerocenter normalization
L = FeatureInputLayer([1 2], 'zerocenter');
L.Mean = [0.5; 0.5];

% Create Star input
lb = [0; 0];
ub = [1; 1];
I_star = Star(lb, ub);

output = L.reach(I_star, 'exact-star');
assert(isa(output, 'Star'), 'reach should return Star');

%% Test 10: FeatureInputLayer toGPU
L = FeatureInputLayer([4 1], 'zerocenter');
L.Mean = [1; 2; 3; 4];
L.StandardDeviation = [0.5; 0.5; 0.5; 0.5];
L.Min = [0; 0; 0; 0];
L.Max = [10; 10; 10; 10];

if gpuDeviceCount > 0
    L_gpu = L.toGPU();
    assert(isa(L_gpu.Mean, 'gpuArray'), 'Mean should be gpuArray');
end

%% Test 11: FeatureInputLayer changeParamsPrecision - single
L = FeatureInputLayer([4 1], 'zerocenter');
L.Mean = double([1; 2; 3; 4]);
L.StandardDeviation = double([0.5; 0.5; 0.5; 0.5]);

L_single = L.changeParamsPrecision('single');
assert(isa(L_single.Mean, 'single'), 'Mean should be single precision');

%% Test 12: FeatureInputLayer changeParamsPrecision - double
L = FeatureInputLayer([4 1], 'zerocenter');
L.Mean = single([1; 2; 3; 4]);

L_double = L.changeParamsPrecision('double');
assert(isa(L_double.Mean, 'double'), 'Mean should be double precision');
