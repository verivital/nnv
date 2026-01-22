%% Test Suite for SiLULayer
% Tests for SiLU (Sigmoid Linear Unit / Swish) activation layer
%
% SiLU(x) = x * sigmoid(x) = x / (1 + exp(-x))
%
% Author: NNV Team
% Date: November 2025

%% Test 1: Basic Evaluation - Scalar
fprintf('=== Test 1: Basic Evaluation - Scalar ===\n');

layer = SiLULayer('test_silu');
assert(strcmp(layer.Name, 'test_silu'), 'Layer name not set correctly');

% Test scalar evaluation
x = 0;
y = layer.evaluate(x);
expected = 0;  % SiLU(0) = 0 * sigmoid(0) = 0 * 0.5 = 0
assert(abs(y - expected) < 1e-10, 'SiLU(0) should be 0');

x = 1;
y = layer.evaluate(x);
expected = 1 / (1 + exp(-1));  % SiLU(1) = 1 * sigmoid(1) ≈ 0.7311
assert(abs(y - expected) < 1e-6, 'SiLU(1) incorrect');

x = -1;
y = layer.evaluate(x);
expected = -1 / (1 + exp(1));  % SiLU(-1) = -1 * sigmoid(-1) ≈ -0.2689
assert(abs(y - expected) < 1e-6, 'SiLU(-1) incorrect');

fprintf('Test 1 PASSED: Basic scalar evaluation\n');

%% Test 2: Evaluation - Vector
fprintf('\n=== Test 2: Evaluation - Vector ===\n');

layer = SiLULayer('test_silu');
x = [-2; -1; 0; 1; 2];
y = layer.evaluate(x);

% Compute expected values
expected = x ./ (1 + exp(-x));

assert(all(abs(y - expected) < 1e-10), 'Vector evaluation failed');
fprintf('Test 2 PASSED: Vector evaluation\n');

%% Test 3: Evaluation - Matrix
fprintf('\n=== Test 3: Evaluation - Matrix ===\n');

layer = SiLULayer('test_silu');
x = randn(5, 5);
y = layer.evaluate(x);
expected = x ./ (1 + exp(-x));

assert(all(abs(y(:) - expected(:)) < 1e-10), 'Matrix evaluation failed');
fprintf('Test 3 PASSED: Matrix evaluation\n');

%% Test 4: Evaluation - 3D Array (Image-like)
fprintf('\n=== Test 4: Evaluation - 3D Array ===\n');

layer = SiLULayer('test_silu');
x = randn(4, 4, 3);  % 4x4 image with 3 channels
y = layer.evaluate(x);
expected = x ./ (1 + exp(-x));

assert(all(abs(y(:) - expected(:)) < 1e-10), '3D array evaluation failed');
assert(isequal(size(y), size(x)), 'Output shape mismatch');
fprintf('Test 4 PASSED: 3D array evaluation\n');

%% Test 5: SiLU Properties - Minimum
fprintf('\n=== Test 5: SiLU Properties - Minimum ===\n');

% SiLU has a minimum at approximately x ≈ -1.278
[x_min, y_min] = SiLU.get_minimum();

% Verify this is indeed a minimum by checking nearby points
delta = 0.001;
y_left = SiLU.evaluate(x_min - delta);
y_right = SiLU.evaluate(x_min + delta);
y_center = SiLU.evaluate(x_min);

assert(y_center < y_left, 'Not a minimum: left point is lower');
assert(y_center < y_right, 'Not a minimum: right point is lower');
assert(abs(y_min - y_center) < 1e-6, 'Minimum value incorrect');

fprintf('Test 5 PASSED: SiLU minimum at x=%.4f, y=%.4f\n', x_min, y_min);

%% Test 6: SiLU Gradient
fprintf('\n=== Test 6: SiLU Gradient ===\n');

% Test gradient computation
x_test = [-2, -1, 0, 1, 2];
dy = SiLU.gradient(x_test);

% Verify gradient numerically
delta = 1e-6;
for i = 1:length(x_test)
    x = x_test(i);
    dy_numerical = (SiLU.evaluate(x + delta) - SiLU.evaluate(x - delta)) / (2 * delta);
    assert(abs(dy(i) - dy_numerical) < 1e-4, 'Gradient incorrect at x=%.1f', x);
end

fprintf('Test 6 PASSED: Gradient computation\n');

%% Test 7: Reachability - Star Set Basic
fprintf('\n=== Test 7: Reachability - Star Set Basic ===\n');

layer = SiLULayer('test_silu');

% Create a simple Star set from lower/upper bounds
lb_7 = [-1; -1];
ub_7 = [1; 1];
S_in = Star(lb_7, ub_7);

% Compute reachability
S_out = layer.reach_star_single_input(S_in, 'approx-star', 0, [], 'linprog');

assert(~isempty(S_out), 'Reachability returned empty');
assert(isa(S_out, 'Star'), 'Output should be a Star');
fprintf('Test 7 PASSED: Star reachability basic\n');

%% Test 8: Reachability - Soundness Check
fprintf('\n=== Test 8: Reachability - Soundness Check ===\n');

layer = SiLULayer('test_silu');

% Create input Star with known bounds
lb = [-1; -0.5];
ub = [1; 0.5];
S_in = Star(lb, ub);

% Compute reachability
S_out = layer.reach_star_single_input(S_in, 'approx-star', 0, [], 'linprog');

% Sample random points and verify they're within output bounds
n_samples = 100;
all_contained = true;

for i = 1:n_samples
    % Random input within bounds [lb, ub]
    x_sample = lb + (ub - lb) .* rand(2, 1);

    % Apply SiLU
    y_sample = SiLU.evaluate(x_sample);

    % Check if output is within reach set bounds
    [out_lb, out_ub] = S_out.getRanges;

    tol = 1e-4;
    if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
        all_contained = false;
        fprintf('Sample %d outside bounds!\n', i);
        fprintf('  y_sample: [%.6f, %.6f]\n', y_sample(1), y_sample(2));
        fprintf('  out_lb: [%.6f, %.6f]\n', out_lb(1), out_lb(2));
        fprintf('  out_ub: [%.6f, %.6f]\n', out_ub(1), out_ub(2));
    end
end

assert(all_contained, 'Soundness violation: some samples outside reach bounds');
fprintf('Test 8 PASSED: Soundness check (%d samples)\n', n_samples);

%% Test 9: Reachability - ImageStar
fprintf('\n=== Test 9: Reachability - ImageStar ===\n');

layer = SiLULayer('test_silu');

% Create ImageStar
h = 2; w = 2; c = 1;
lb_img = -ones(h, w, c);
ub_img = ones(h, w, c);
IS_in = ImageStar(lb_img, ub_img);

% Compute reachability
IS_out = layer.reach_star_single_input(IS_in, 'approx-star', 0, [], 'linprog');

assert(~isempty(IS_out), 'ImageStar reachability returned empty');
assert(isa(IS_out, 'ImageStar'), 'Output should be an ImageStar');
assert(IS_out.height == h, 'Output height mismatch');
assert(IS_out.width == w, 'Output width mismatch');

fprintf('Test 9 PASSED: ImageStar reachability\n');

%% Test 10: Reachability - ImageStar Soundness
fprintf('\n=== Test 10: Reachability - ImageStar Soundness ===\n');

layer = SiLULayer('test_silu');

% Create ImageStar with smaller bounds
lb_img = -0.5 * ones(2, 2, 1);
ub_img = 0.5 * ones(2, 2, 1);
IS_in = ImageStar(lb_img, ub_img);

% Compute reachability
IS_out = layer.reach_star_single_input(IS_in, 'approx-star', 0, [], 'linprog');

% Get output bounds
[out_lb, out_ub] = IS_out.getRanges;

% Sample and verify soundness
n_samples = 50;
all_contained = true;

for i = 1:n_samples
    % Random input
    x_sample = lb_img + (ub_img - lb_img) .* rand(size(lb_img));

    % Apply SiLU
    y_sample = SiLU.evaluate(x_sample);

    tol = 1e-4;
    if any(y_sample(:) < out_lb(:) - tol) || any(y_sample(:) > out_ub(:) + tol)
        all_contained = false;
    end
end

assert(all_contained, 'ImageStar soundness violation');
fprintf('Test 10 PASSED: ImageStar soundness (%d samples)\n', n_samples);

%% Test 11: Asymptotic Behavior - Large Positive
fprintf('\n=== Test 11: Asymptotic Behavior ===\n');

% For large positive x, SiLU(x) ≈ x
x_large = 10;
y_large = SiLU.evaluate(x_large);
assert(abs(y_large - x_large) < 0.001, 'SiLU should approach x for large positive x');

% For large negative x, SiLU(x) ≈ 0
x_neg = -10;
y_neg = SiLU.evaluate(x_neg);
assert(abs(y_neg) < 0.001, 'SiLU should approach 0 for large negative x');

fprintf('Test 11 PASSED: Asymptotic behavior\n');

%% Test 12: Default Constructor
fprintf('\n=== Test 12: Default Constructor ===\n');

layer_default = SiLULayer();
assert(strcmp(layer_default.Name, 'silu'), 'Default name should be "silu"');
assert(layer_default.NumInputs == 1, 'Default NumInputs should be 1');
assert(layer_default.NumOutputs == 1, 'Default NumOutputs should be 1');

fprintf('Test 12 PASSED: Default constructor\n');

%% Test 13: Comparison with Sigmoid Multiplication
fprintf('\n=== Test 13: Comparison with x * sigmoid(x) ===\n');

% Verify SiLU matches x * sigmoid(x) definition
x_test = linspace(-5, 5, 100)';
y_silu = SiLU.evaluate(x_test);
y_expected = x_test .* (1 ./ (1 + exp(-x_test)));

assert(all(abs(y_silu - y_expected) < 1e-10), 'SiLU should equal x * sigmoid(x)');
fprintf('Test 13 PASSED: SiLU = x * sigmoid(x)\n');

%% Test 14: Layer in Network Context
fprintf('\n=== Test 14: Layer in Network Context ===\n');

% Create a simple network with SiLU
W1 = randn(3, 4);
b1 = randn(3, 1);
fc1 = FullyConnectedLayer(W1, b1);
silu = SiLULayer('silu');
W2 = randn(2, 3);
b2 = randn(2, 1);
fc2 = FullyConnectedLayer(W2, b2);

% Test forward pass
x = randn(4, 1);
h1 = fc1.evaluate(x);
h2 = silu.evaluate(h1);
y = fc2.evaluate(h2);

assert(all(size(y) == [2, 1]), 'Output shape incorrect');
fprintf('Test 14 PASSED: Layer in network context\n');

%% Test 15: Visual Soundness Verification with Figure
fprintf('\n=== Test 15: Visual Soundness Verification ===\n');

layer15 = SiLULayer('test_silu_visual');

% Create 4D input for better visualization
input_dim = 4;
center = randn(input_dim, 1);
eps = 0.1;
lb15 = center - eps;
ub15 = center + eps;
S_in15 = Star(lb15, ub15);

% Compute reachable set
S_out15 = layer15.reach_star_single_input(S_in15, 'approx-star', 0, [], 'linprog');
[out_lb, out_ub] = S_out15.getRanges();
out_lb = out_lb(:);
out_ub = out_ub(:);

% Sample points and evaluate
n_samples = 50;
samples_out = zeros(input_dim, n_samples);
for i = 1:n_samples
    x_sample = lb15 + (ub15 - lb15) .* rand(input_dim, 1);
    samples_out(:, i) = layer15.evaluate(x_sample);
end

% Check containment
all_contained = true;
tol = 1e-6;
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
for d = 1:input_dim
    % Draw bound box using fill()
    fill([d-0.4, d+0.4, d+0.4, d-0.4], ...
         [out_lb(d), out_lb(d), out_ub(d), out_ub(d)], ...
         colors(d,:), 'FaceAlpha', 0.3, 'EdgeColor', colors(d,:), 'LineWidth', 1.5);

    % Overlay sample points using scatter()
    x_jitter = d + (rand(1, n_samples) - 0.5) * 0.6;
    scatter(x_jitter, samples_out(d, :), 30, colors(d,:), 'filled', 'MarkerFaceAlpha', 0.7);
end

% Status display
if all_contained
    title_str = sprintf('SiLULayer Soundness: SOUND (%d samples)', n_samples);
    title(title_str, 'Color', [0 0.6 0], 'FontWeight', 'bold', 'FontSize', 12);
else
    title_str = sprintf('SiLULayer Soundness: VIOLATION');
    title(title_str, 'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 12);
end

xlabel('Output Dimension', 'FontSize', 11);
ylabel('Value', 'FontSize', 11);
xlim([0.5, input_dim + 0.5]);
xticks(1:input_dim);
grid on;
hold off;

% Save figure using test infrastructure
try
    save_test_figure(fig, 'test_SiLULayer', 'soundness_containment', 15, 'subdir', 'nn/layers');
catch
    % If save_test_figure not available, save directly
    fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'results', 'tests', 'figures', 'nn', 'layers');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    saveas(fig, fullfile(fig_dir, 'test_SiLULayer_soundness_containment_15.png'));
    fprintf('  Figure saved to: %s\n', fig_dir);
end
if isvalid(fig)
    close(fig);
end

assert(all_contained, 'Visual soundness verification failed');
fprintf('Test 15 PASSED: Visual soundness verification\n');

%% Summary
fprintf('\n=== All SiLULayer Tests PASSED ===\n');
fprintf('Total: 15 tests\n');
