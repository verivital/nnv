%% Test Suite for SwiGLULayer
% Tests for SwiGLU (Gated Linear Unit with SiLU) layer
%
% SwiGLU(x) = SiLU(x @ W_gate + b_gate) .* (x @ W_up + b_up)
%
% Author: NNV Team
% Date: November 2025

%% Test 1: Basic Construction
fprintf('=== Test 1: Basic Construction ===\n');

% Create random weights
input_dim = 4;
hidden_dim = 8;

W_gate = randn(input_dim, hidden_dim);
W_up = randn(input_dim, hidden_dim);

layer = SwiGLULayer(W_gate, W_up);

assert(strcmp(layer.Name, 'swiglu'), 'Default name should be swiglu');
assert(layer.InputDim == input_dim, 'InputDim mismatch');
assert(layer.HiddenDim == hidden_dim, 'HiddenDim mismatch');
assert(all(layer.b_gate == 0), 'Default bias should be zero');
assert(all(layer.b_up == 0), 'Default bias should be zero');

fprintf('Test 1 PASSED: Basic construction\n');

%% Test 2: Construction with Biases
fprintf('\n=== Test 2: Construction with Biases ===\n');

input_dim = 4;
hidden_dim = 8;
W_gate = randn(input_dim, hidden_dim);
W_up = randn(input_dim, hidden_dim);
b_gate = randn(hidden_dim, 1);
b_up = randn(hidden_dim, 1);

layer2 = SwiGLULayer(W_gate, b_gate, W_up, b_up);

assert(all(layer2.b_gate == b_gate), 'Gate bias not set correctly');
assert(all(layer2.b_up == b_up), 'Up bias not set correctly');

fprintf('Test 2 PASSED: Construction with biases\n');

%% Test 3: Basic Evaluation
fprintf('\n=== Test 3: Basic Evaluation ===\n');

% Simple weights for predictable output
W_gate_simple = eye(3, 3);
W_up_simple = eye(3, 3);
b_gate_simple = zeros(3, 1);
b_up_simple = zeros(3, 1);

layer3 = SwiGLULayer(W_gate_simple, b_gate_simple, W_up_simple, b_up_simple);

x = [1; 0; -1];
y = layer3.evaluate(x);

% Expected: SiLU(x) .* x = x .* sigmoid(x) .* x = x^2 .* sigmoid(x)
expected = SiLU.evaluate(x) .* x;

assert(all(abs(y - expected) < 1e-10), 'Basic evaluation failed');
fprintf('Test 3 PASSED: Basic evaluation\n');

%% Test 4: Evaluation with Random Weights
fprintf('\n=== Test 4: Evaluation with Random Weights ===\n');

rng(42);
W_gate_rand = randn(5, 8);
b_gate_rand = randn(8, 1);
W_up_rand = randn(5, 8);
b_up_rand = randn(8, 1);

layer4 = SwiGLULayer(W_gate_rand, b_gate_rand, W_up_rand, b_up_rand);

x_rand = randn(5, 1);
y_rand = layer4.evaluate(x_rand);

% Manual computation
gate_linear = W_gate_rand' * x_rand + b_gate_rand;
gate = SiLU.evaluate(gate_linear);
up = W_up_rand' * x_rand + b_up_rand;
expected_rand = gate .* up;

assert(all(abs(y_rand - expected_rand) < 1e-10), 'Random weight evaluation failed');
fprintf('Test 4 PASSED: Evaluation with random weights\n');

%% Test 5: Output Shape
fprintf('\n=== Test 5: Output Shape ===\n');

input_dim5 = 10;
hidden_dim5 = 20;
W_gate5 = randn(input_dim5, hidden_dim5);
W_up5 = randn(input_dim5, hidden_dim5);

layer5 = SwiGLULayer(W_gate5, W_up5);

x5 = randn(input_dim5, 1);
y5 = layer5.evaluate(x5);

assert(length(y5) == hidden_dim5, 'Output dimension should be hidden_dim');
fprintf('Test 5 PASSED: Output shape correct\n');

%% Test 6: Reachability - Star Set Basic
fprintf('\n=== Test 6: Reachability - Star Set Basic ===\n');

% Small network for tractable reachability
W_gate6 = randn(2, 3) * 0.5;
W_up6 = randn(2, 3) * 0.5;
b_gate6 = zeros(3, 1);
b_up6 = zeros(3, 1);

layer6 = SwiGLULayer(W_gate6, b_gate6, W_up6, b_up6);

% Create input Star
lb6 = [-1; -1];
ub6 = [1; 1];
S_in = Star(lb6, ub6);

% Compute reachability
S_out = layer6.reach(S_in, 'approx-star');

assert(~isempty(S_out), 'Reachability returned empty');
assert(isa(S_out, 'Star'), 'Output should be a Star');
assert(S_out.dim == 3, 'Output dimension should be hidden_dim');

fprintf('Test 6 PASSED: Star reachability basic\n');

%% Test 7: Reachability - Soundness Check
fprintf('\n=== Test 7: Reachability - Soundness Check ===\n');

% Recreate layer from Test 6 (tests run independently)
W_gate6 = randn(2, 3) * 0.5;
W_up6 = randn(2, 3) * 0.5;
b_gate6 = zeros(3, 1);
b_up6 = zeros(3, 1);
layer6 = SwiGLULayer(W_gate6, b_gate6, W_up6, b_up6);

lb7 = [-0.5; -0.5];
ub7 = [0.5; 0.5];
S_in7 = Star(lb7, ub7);

S_out7 = layer6.reach(S_in7, 'approx-star');

% Get output bounds
[out_lb, out_ub] = S_out7.getRanges;

% Sample random points and verify soundness
n_samples = 100;
all_contained = true;

for i = 1:n_samples
    % Random input
    x_sample = lb7 + (ub7 - lb7) .* rand(2, 1);

    % Evaluate
    y_sample = layer6.evaluate(x_sample);

    % Check bounds
    tol = 1e-5;
    if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
        all_contained = false;
        fprintf('Sample %d outside bounds!\n', i);
        fprintf('  y: [%.6f, %.6f, %.6f]\n', y_sample(1), y_sample(2), y_sample(3));
        fprintf('  lb: [%.6f, %.6f, %.6f]\n', out_lb(1), out_lb(2), out_lb(3));
        fprintf('  ub: [%.6f, %.6f, %.6f]\n', out_ub(1), out_ub(2), out_ub(3));
    end
end

assert(all_contained, 'Soundness violation: samples outside bounds');
fprintf('Test 7 PASSED: Soundness check (%d samples)\n', n_samples);

%% Test 8: Zero Input
fprintf('\n=== Test 8: Zero Input ===\n');

% Recreate layer from Test 6 (tests run independently)
W_gate6 = randn(2, 3) * 0.5;
W_up6 = randn(2, 3) * 0.5;
b_gate6 = zeros(3, 1);
b_up6 = zeros(3, 1);
layer6 = SwiGLULayer(W_gate6, b_gate6, W_up6, b_up6);

x_zero = zeros(2, 1);
y_zero = layer6.evaluate(x_zero);

% SwiGLU(0) should be SiLU(b_gate) .* b_up
expected_zero = SiLU.evaluate(b_gate6) .* b_up6;

assert(all(abs(y_zero - expected_zero) < 1e-10), 'Zero input evaluation failed');
fprintf('Test 8 PASSED: Zero input\n');

%% Test 9: Named Constructor
fprintf('\n=== Test 9: Named Constructor ===\n');

% Define weights for this test
W_gate9 = randn(2, 3) * 0.5;
W_up9 = randn(2, 3) * 0.5;
b_gate9 = zeros(3, 1);
b_up9 = zeros(3, 1);

layer9 = SwiGLULayer(W_gate9, b_gate9, W_up9, b_up9, 'my_swiglu');
assert(strcmp(layer9.Name, 'my_swiglu'), 'Name not set correctly');

fprintf('Test 9 PASSED: Named constructor\n');

%% Test 10: set_weights Method
fprintf('\n=== Test 10: set_weights Method ===\n');

layer10 = SwiGLULayer();
W_g = randn(4, 6);
b_g = randn(6, 1);
W_u = randn(4, 6);
b_u = randn(6, 1);

layer10 = layer10.set_weights(W_g, b_g, W_u, b_u);

assert(layer10.InputDim == 4, 'InputDim not set');
assert(layer10.HiddenDim == 6, 'HiddenDim not set');
assert(all(all(layer10.W_gate == W_g)), 'W_gate not set');

fprintf('Test 10 PASSED: set_weights method\n');

%% Test 11: Interval Multiply Helper
fprintf('\n=== Test 11: Interval Multiply Helper ===\n');

% Test interval multiplication
[lb, ub] = SwiGLULayer.interval_multiply(-1, 2, -3, 4);
% All products: (-1)*(-3)=3, (-1)*4=-4, 2*(-3)=-6, 2*4=8
assert(lb == -6, 'Interval multiply lower bound wrong');
assert(ub == 8, 'Interval multiply upper bound wrong');

% Test with positive intervals
[lb2, ub2] = SwiGLULayer.interval_multiply(1, 2, 3, 4);
assert(lb2 == 3, 'Positive interval lb wrong');
assert(ub2 == 8, 'Positive interval ub wrong');

fprintf('Test 11 PASSED: Interval multiply helper\n');

%% Test 12: SiLU Bounds Helper
fprintf('\n=== Test 12: SiLU Bounds Helper ===\n');

% Test SiLU bounds
[lb_silu, ub_silu] = SwiGLULayer.silu_bounds(-2, 2);

% Sample points and verify bounds
x_test = linspace(-2, 2, 100);
y_test = SiLU.evaluate(x_test);

assert(lb_silu <= min(y_test) + 1e-6, 'SiLU lower bound too high');
assert(ub_silu >= max(y_test) - 1e-6, 'SiLU upper bound too low');

fprintf('Test 12 PASSED: SiLU bounds helper\n');

%% Test 13: Layer in Network Context
fprintf('\n=== Test 13: Layer in Network Context ===\n');

% Create a mini network with SwiGLU
% NNV FullyConnectedLayer: W is [output_dim x input_dim]
W1 = randn(3, 4);  % 4 -> 3
b1 = randn(3, 1);
fc1 = FullyConnectedLayer('fc1', W1, b1);

W_g = randn(3, 5);  % SwiGLU: input_dim=3, hidden_dim=5
W_u = randn(3, 5);
swiglu = SwiGLULayer(W_g, W_u);

W2 = randn(2, 5);  % 5 -> 2
b2 = randn(2, 1);
fc2 = FullyConnectedLayer('fc2', W2, b2);

% Forward pass
x = randn(4, 1);
h1 = fc1.evaluate(x);
h2 = swiglu.evaluate(h1);
y = fc2.evaluate(h2);

assert(all(size(y) == [2, 1]), 'Network output shape wrong');
fprintf('Test 13 PASSED: Layer in network context\n');

%% Test 14: Visual Soundness Verification with Figure
fprintf('\n=== Test 14: Visual Soundness Verification ===\n');

% Create SwiGLU layer with small dimensions for visualization
input_dim14 = 3;
hidden_dim14 = 4;
W_gate14 = randn(input_dim14, hidden_dim14) * 0.5;
W_up14 = randn(input_dim14, hidden_dim14) * 0.5;
b_gate14 = zeros(hidden_dim14, 1);
b_up14 = zeros(hidden_dim14, 1);
layer14 = SwiGLULayer(W_gate14, b_gate14, W_up14, b_up14, 'swiglu_visual');

% Create input Star with small perturbation
center14 = randn(input_dim14, 1);
eps14 = 0.05;
lb14 = center14 - eps14;
ub14 = center14 + eps14;
S_in14 = Star(lb14, ub14);

% Compute reachable set
S_out14 = layer14.reach(S_in14, 'approx-star');
[out_lb, out_ub] = S_out14.getRanges();
out_lb = out_lb(:);
out_ub = out_ub(:);

% Sample points and evaluate
n_samples = 50;
samples_out = zeros(hidden_dim14, n_samples);
for i = 1:n_samples
    x_sample = lb14 + (ub14 - lb14) .* rand(input_dim14, 1);
    samples_out(:, i) = layer14.evaluate(x_sample);
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
for d = 1:hidden_dim14
    fill([d-0.4, d+0.4, d+0.4, d-0.4], ...
         [out_lb(d), out_lb(d), out_ub(d), out_ub(d)], ...
         colors(d,:), 'FaceAlpha', 0.3, 'EdgeColor', colors(d,:), 'LineWidth', 1.5);
    x_jitter = d + (rand(1, n_samples) - 0.5) * 0.6;
    scatter(x_jitter, samples_out(d, :), 30, colors(d,:), 'filled', 'MarkerFaceAlpha', 0.7);
end

% Status display
if all_contained
    title_str = sprintf('SwiGLULayer Soundness: SOUND (%d samples)', n_samples);
    title(title_str, 'Color', [0 0.6 0], 'FontWeight', 'bold', 'FontSize', 12);
else
    title_str = sprintf('SwiGLULayer Soundness: VIOLATION');
    title(title_str, 'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 12);
end

xlabel('Output Dimension', 'FontSize', 11);
ylabel('Value', 'FontSize', 11);
xlim([0.5, hidden_dim14 + 0.5]);
xticks(1:hidden_dim14);
grid on;
hold off;

% Save figure
try
    save_test_figure(fig, 'test_SwiGLULayer', 'soundness_containment', 14, 'subdir', 'nn/layers');
catch
    fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'tests', 'figures', 'nn', 'layers');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    saveas(fig, fullfile(fig_dir, 'test_SwiGLULayer_soundness_containment_14.png'));
    fprintf('  Figure saved to: %s\n', fig_dir);
end
if isvalid(fig)
    close(fig);
end

assert(all_contained, 'Visual soundness verification failed');
fprintf('Test 14 PASSED: Visual soundness verification\n');

%% Summary
fprintf('\n=== All SwiGLULayer Tests PASSED ===\n');
fprintf('Total: 14 tests\n');
