%% Test Suite for RMSNormLayer
% Tests for RMS Normalization layer
%
% RMSNorm(x) = x / sqrt(mean(x^2) + eps) * gamma
%
% Author: NNV Team
% Date: November 2025

%% Test 1: Basic Construction - Dimension
fprintf('=== Test 1: Basic Construction - Dimension ===\n');

dim = 5;
layer = RMSNormLayer(dim);

assert(layer.Dim == dim, 'Dimension not set correctly');
assert(all(layer.Gamma == 1), 'Default gamma should be ones');
assert(layer.Epsilon == 1e-6, 'Default epsilon incorrect');

fprintf('Test 1 PASSED: Basic construction with dimension\n');

%% Test 2: Construction with Gamma
fprintf('\n=== Test 2: Construction with Gamma ===\n');

gamma = [1; 2; 0.5; 1.5; 0.8];
layer2 = RMSNormLayer(gamma);

assert(layer2.Dim == 5, 'Dimension from gamma incorrect');
assert(all(layer2.Gamma == gamma), 'Gamma not set correctly');

fprintf('Test 2 PASSED: Construction with gamma\n');

%% Test 3: Construction with Epsilon
fprintf('\n=== Test 3: Construction with Epsilon ===\n');

gamma3 = ones(4, 1);
eps3 = 1e-8;
layer3 = RMSNormLayer(gamma3, eps3);

assert(layer3.Epsilon == eps3, 'Epsilon not set correctly');

fprintf('Test 3 PASSED: Construction with epsilon\n');

%% Test 4: Basic Evaluation
fprintf('\n=== Test 4: Basic Evaluation ===\n');

gamma4 = ones(3, 1);
layer4 = RMSNormLayer(gamma4);

x = [3; 4; 0];  % RMS = sqrt((9+16+0)/3 + eps) = sqrt(25/3 + eps)
y = layer4.evaluate(x);

% Manual computation
rms = sqrt(mean(x.^2) + 1e-6);
expected = x / rms;

assert(all(abs(y - expected) < 1e-10), 'Basic evaluation failed');
fprintf('Test 4 PASSED: Basic evaluation\n');

%% Test 5: Evaluation with Custom Gamma
fprintf('\n=== Test 5: Evaluation with Custom Gamma ===\n');

gamma5 = [2; 0.5; 1];
layer5 = RMSNormLayer(gamma5);

x5 = [1; 2; 3];
y5 = layer5.evaluate(x5);

rms5 = sqrt(mean(x5.^2) + 1e-6);
expected5 = (x5 / rms5) .* gamma5;

assert(all(abs(y5 - expected5) < 1e-10), 'Evaluation with gamma failed');
fprintf('Test 5 PASSED: Evaluation with custom gamma\n');

%% Test 6: Normalization Property
fprintf('\n=== Test 6: Normalization Property ===\n');

% After RMSNorm with gamma=1, output should have RMS close to 1
gamma6 = ones(4, 1);
layer6 = RMSNormLayer(gamma6);

x6 = randn(4, 1) * 10;  % Random input with large variance
y6 = layer6.evaluate(x6);

% RMS of output should be close to 1
rms_out = sqrt(mean(y6.^2));
assert(abs(rms_out - 1) < 0.01, 'Output RMS should be ~1');

fprintf('Test 6 PASSED: Normalization property (RMS ≈ 1)\n');

%% Test 7: Scale Invariance
fprintf('\n=== Test 7: Scale Invariance ===\n');

% RMSNorm should be scale-invariant
gamma7 = ones(3, 1);
layer7 = RMSNormLayer(gamma7);

x7 = [1; 2; 3];
y7 = layer7.evaluate(x7);

% Scaled input
x7_scaled = 10 * x7;
y7_scaled = layer7.evaluate(x7_scaled);

% Outputs should be identical (up to numerical precision)
assert(all(abs(y7 - y7_scaled) < 1e-6), 'RMSNorm should be scale invariant');

fprintf('Test 7 PASSED: Scale invariance\n');

%% Test 8: Reachability - Star Set Basic
fprintf('\n=== Test 8: Reachability - Star Set Basic ===\n');

gamma8 = ones(2, 1);
layer8 = RMSNormLayer(gamma8);

lb8 = [-1; -1];
ub8 = [1; 1];
S_in = Star(lb8, ub8);

S_out = layer8.reach(S_in, 'approx-star');

assert(~isempty(S_out), 'Reachability returned empty');
assert(isa(S_out, 'Star'), 'Output should be a Star');
assert(S_out.dim == 2, 'Output dimension should match input');

fprintf('Test 8 PASSED: Star reachability basic\n');

%% Test 9: Reachability - Soundness Check
fprintf('\n=== Test 9: Reachability - Soundness Check ===\n');

gamma9 = [1.5; 0.8; 1.2];
layer9 = RMSNormLayer(gamma9);

lb9 = [0.5; -0.5; 0.1];
ub9 = [1.5; 0.5; 0.9];
S_in9 = Star(lb9, ub9);

S_out9 = layer9.reach(S_in9, 'approx-star');

% Get output bounds
[out_lb, out_ub] = S_out9.getRanges;

% Sample and verify soundness
n_samples = 100;
all_contained = true;

for i = 1:n_samples
    x_sample = lb9 + (ub9 - lb9) .* rand(3, 1);
    y_sample = layer9.evaluate(x_sample);

    tol = 1e-5;
    if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
        all_contained = false;
        fprintf('Sample %d outside bounds!\n', i);
    end
end

assert(all_contained, 'Soundness violation');
fprintf('Test 9 PASSED: Soundness check (%d samples)\n', n_samples);

%% Test 10: Zero Input
fprintf('\n=== Test 10: Near-Zero Input ===\n');

gamma10 = ones(2, 1);
eps10 = 1e-6;
layer10 = RMSNormLayer(gamma10, eps10);

% Near-zero input (can't be exactly zero due to division)
x10 = [1e-8; 1e-8];
y10 = layer10.evaluate(x10);

% Should not produce NaN or Inf
assert(~any(isnan(y10)), 'NaN in output');
assert(~any(isinf(y10)), 'Inf in output');

fprintf('Test 10 PASSED: Near-zero input handled\n');

%% Test 11: Negative Values
fprintf('\n=== Test 11: Negative Values ===\n');

gamma11 = ones(3, 1);
layer11 = RMSNormLayer(gamma11);

x11_neg = [-1; -2; -3];
y11_neg = layer11.evaluate(x11_neg);

x11_pos = [1; 2; 3];
y11_pos = layer11.evaluate(x11_pos);

% Output signs should match input signs
assert(all(sign(y11_neg) == sign(x11_neg)), 'Sign preservation failed for negative');
assert(all(sign(y11_pos) == sign(x11_pos)), 'Sign preservation failed for positive');

fprintf('Test 11 PASSED: Negative values handled correctly\n');

%% Test 12: Named Constructor
fprintf('\n=== Test 12: Named Constructor ===\n');

gamma12 = ones(4, 1);
layer12 = RMSNormLayer(gamma12, 1e-6, 'my_rmsnorm');

assert(strcmp(layer12.Name, 'my_rmsnorm'), 'Name not set correctly');

fprintf('Test 12 PASSED: Named constructor\n');

%% Test 13: set_gamma Method
fprintf('\n=== Test 13: set_gamma Method ===\n');

layer13 = RMSNormLayer(3);  % Start with dim=3
new_gamma = [0.5; 1.0; 1.5; 2.0];

layer13 = layer13.set_gamma(new_gamma);

assert(layer13.Dim == 4, 'Dimension not updated');
assert(all(layer13.Gamma == new_gamma), 'Gamma not updated');

fprintf('Test 13 PASSED: set_gamma method\n');

%% Test 14: Comparison with Manual RMSNorm
fprintf('\n=== Test 14: Comparison with Manual RMSNorm ===\n');

gamma14 = randn(5, 1);
eps14 = 1e-5;
layer14 = RMSNormLayer(gamma14, eps14);

x14 = randn(5, 1);
y14 = layer14.evaluate(x14);

% Manual computation
rms14 = sqrt(mean(x14.^2) + eps14);
y14_manual = (x14 / rms14) .* gamma14;

assert(all(abs(y14 - y14_manual) < 1e-10), 'Manual comparison failed');

fprintf('Test 14 PASSED: Matches manual computation\n');

%% Test 15: Layer in Network Context
fprintf('\n=== Test 15: Layer in Network Context ===\n');

% Mini network with RMSNorm
% NNV FullyConnectedLayer: W is [output_dim x input_dim]
W1 = randn(8, 4);  % 4 -> 8
b1 = randn(8, 1);
fc1 = FullyConnectedLayer('fc1', W1, b1);

gamma_norm = ones(8, 1);
rmsnorm = RMSNormLayer(gamma_norm);

W2 = randn(2, 8);  % 8 -> 2
b2 = randn(2, 1);
fc2 = FullyConnectedLayer('fc2', W2, b2);

x15 = randn(4, 1);
h1 = fc1.evaluate(x15);
h2 = rmsnorm.evaluate(h1);
y15 = fc2.evaluate(h2);

assert(all(size(y15) == [2, 1]), 'Network output shape wrong');

fprintf('Test 15 PASSED: Layer in network context\n');

%% Test 16: Visual Soundness Verification with Figure
fprintf('\n=== Test 16: Visual Soundness Verification ===\n');

% Create RMSNorm layer for visualization
input_dim16 = 4;
gamma16 = ones(input_dim16, 1);
eps16 = 1e-5;
layer16 = RMSNormLayer(gamma16, eps16, 'rmsnorm_visual');

% Create input Star with perturbation
center16 = 0.5 * ones(input_dim16, 1) + 0.1 * randn(input_dim16, 1);
perturbation = 0.1;
lb16 = center16 - perturbation;
ub16 = center16 + perturbation;
S_in16 = Star(lb16, ub16);

% Compute reachable set
S_out16 = layer16.reach(S_in16, 'approx-star');
[out_lb, out_ub] = S_out16.getRanges();
out_lb = out_lb(:);
out_ub = out_ub(:);

% Sample points and evaluate
n_samples = 50;
samples_out = zeros(input_dim16, n_samples);
for i = 1:n_samples
    x_sample = lb16 + (ub16 - lb16) .* rand(input_dim16, 1);
    samples_out(:, i) = layer16.evaluate(x_sample);
end

% Check containment
all_contained = true;
tol = 1e-5;
for i = 1:n_samples
    if any(samples_out(:,i) < out_lb - tol) || any(samples_out(:,i) > out_ub + tol)
        all_contained = false;
        break;
    end
end

% Create figure (AG News containment style)
fig = figure('Position', [100 100 800 500], 'Visible', 'off');
hold on;

% Define colors for each dimension
colors = [0.2 0.4 0.8;   % Blue
          0.8 0.2 0.2;   % Red
          0.2 0.7 0.3;   % Green
          0.7 0.4 0.9];  % Purple

% Plot bounds as boxes and samples as scatter points
for d = 1:input_dim16
    fill([d-0.4, d+0.4, d+0.4, d-0.4], ...
         [out_lb(d), out_lb(d), out_ub(d), out_ub(d)], ...
         colors(d,:), 'FaceAlpha', 0.3, 'EdgeColor', colors(d,:), 'LineWidth', 1.5);
    x_jitter = d + (rand(1, n_samples) - 0.5) * 0.6;
    scatter(x_jitter, samples_out(d, :), 30, colors(d,:), 'filled', 'MarkerFaceAlpha', 0.7);
end

% Status display
if all_contained
    title_str = sprintf('RMSNormLayer Soundness: SOUND (%d samples)', n_samples);
    title(title_str, 'Color', [0 0.6 0], 'FontWeight', 'bold', 'FontSize', 12);
else
    title_str = sprintf('RMSNormLayer Soundness: VIOLATION');
    title(title_str, 'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 12);
end

xlabel('Output Dimension', 'FontSize', 11);
ylabel('Value', 'FontSize', 11);
xlim([0.5, input_dim16 + 0.5]);
xticks(1:input_dim16);
grid on;
hold off;

% Save figure
try
    save_test_figure(fig, 'test_RMSNormLayer', 'soundness_containment', 16, 'subdir', 'nn/layers');
catch
    fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'tests', 'figures', 'nn', 'layers');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    saveas(fig, fullfile(fig_dir, 'test_RMSNormLayer_soundness_containment_16.png'));
    fprintf('  Figure saved to: %s\n', fig_dir);
end
if isvalid(fig)
    close(fig);
end

assert(all_contained, 'Visual soundness verification failed');
fprintf('Test 16 PASSED: Visual soundness verification\n');

%% Summary
fprintf('\n=== All RMSNormLayer Tests PASSED ===\n');
fprintf('Total: 16 tests\n');
