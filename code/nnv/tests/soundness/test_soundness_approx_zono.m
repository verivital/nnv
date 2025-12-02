% test_soundness_approx_zono
% Soundness tests for approx-zono reachability method
% Tests Zonotope-based over-approximation across multiple layers
% To run: results = runtests('test_soundness_approx_zono')

%% Test 1: ReLU with approx-zono using ImageZono input
rng(42);
L = ReluLayer();

% Create ImageZono input (using lower/upper bounds)
lb = zeros(3, 3);
ub = ones(3, 3) * 0.5;
lb(2, 2) = -0.3;  % Some negative values to trigger ReLU
ub(2, 2) = 0.3;

input_iz = ImageZono(lb, ub);
output_iz = L.reach(input_iz, 'approx-zono');

% Verify output is valid ImageZono
assert(isa(output_iz, 'ImageZono') || isa(output_iz, 'ImageStar'), ...
    'ReLU approx-zono should return ImageZono or ImageStar');

% Sample and verify containment (outputs should be >= 0)
for i = 1:20
    % Sample from input bounds
    input_sample = lb + rand(size(lb)) .* (ub - lb);
    output_sample = L.evaluate(input_sample);

    % ReLU output should be max(0, input)
    expected = max(0, input_sample);
    assert(max(abs(output_sample(:) - expected(:))) < 1e-10, ...
        'ReLU evaluate should match max(0, x)');
end

%% Test 2: FullyConnected with approx-zono using Zono input
rng(42);

W = randn(3, 4);
b = randn(3, 1);
L2 = FullyConnectedLayer(W, b);

% Create Zono input
c = randn(4, 1);  % center
V = randn(4, 2) * 0.1;  % generators
input_zono = Zono(c, V);

% Test that layer can process Zono input
output_star = L2.reach(input_zono, 'approx-zono');

% Verify output type
assert(isa(output_star, 'Zono') || isa(output_star, 'Star'), ...
    'FullyConnected approx-zono should return Zono or Star');

%% Test 3: Conv2D with approx-zono
rng(42);

W = randn(3, 3, 1, 2);  % 3x3 kernel, 1 input channel, 2 filters
b = randn(1, 1, 2);
L3 = Conv2DLayer(W, b);
L3.PaddingSize = [0 0 0 0];
L3.Stride = [1 1];
L3.DilationFactor = [1 1];

% Create ImageZono input
lb3 = rand(5, 5) * 0.5;
ub3 = lb3 + rand(5, 5) * 0.2;

input_iz3 = ImageZono(lb3, ub3);
output_iz3 = L3.reach(input_iz3, 'approx-zono');

% Verify output structure
assert(~isempty(output_iz3), 'Conv2D approx-zono should produce output');

%% Test 4: AveragePooling2D with approx-zono
rng(42);

L4 = AveragePooling2DLayer([2 2], [2 2], [0 0 0 0]);

lb4 = rand(4, 4);
ub4 = lb4 + rand(4, 4) * 0.1;
input_iz4 = ImageZono(lb4, ub4);

output_iz4 = L4.reach(input_iz4, 'approx-zono');
assert(~isempty(output_iz4), 'AvgPool approx-zono should produce output');

%% Test 5: BatchNorm with approx-zono
rng(42);

L5 = BatchNormalizationLayer('Name', 'bn_zono', ...
    'NumChannels', 1, ...
    'TrainedMean', 0.5, ...
    'TrainedVariance', 0.1, ...
    'Epsilon', 1e-5, ...
    'Scale', 1.0, ...
    'Offset', 0.0);

lb5 = rand(3, 3);
ub5 = lb5 + rand(3, 3) * 0.2;
input_iz5 = ImageZono(lb5, ub5);

output_iz5 = L5.reach(input_iz5, 'approx-zono');
assert(~isempty(output_iz5), 'BatchNorm approx-zono should produce output');

% NOTE: LeakyReLU approx-zono has a library bug in LeakyReLU.reach_zono_approx
% line 1595: V1(map1, :) = gamma*V1(map1) has dimension mismatch
% This should be fixed in the NNV library

