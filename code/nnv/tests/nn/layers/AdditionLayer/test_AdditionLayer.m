% Test AdditionLayer functionality
% To run: results = runtests('test_AdditionLayer')

%% Test 1: AdditionLayer constructor with single name argument
try
    L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});
    assert(strcmp(L.Name, 'add_layer'));
    assert(L.NumInputs == 2);
    assert(L.NumOutputs == 1);
    assert(length(L.InputNames) == 2);
    assert(strcmp(L.InputNames{1}, 'in1'));
    assert(strcmp(L.InputNames{2}, 'in2'));
catch e
    error('AdditionLayer constructor failed: %s', e.message);
end

%% Test 2: AdditionLayer evaluate with two inputs
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});

% Create two simple matrices
input1 = [1 2 3; 4 5 6];
input2 = [1 1 1; 2 2 2];
expected = [2 3 4; 6 7 8];

inputs = {input1, input2};
output = L.evaluate(inputs);

assert(isequal(output, expected), 'AdditionLayer evaluate failed for two inputs');

%% Test 3: AdditionLayer evaluate with three inputs
L = AdditionLayer('add_layer', 3, 1, {'in1', 'in2', 'in3'}, {'out'});

% Create three simple matrices
input1 = [1 2; 3 4];
input2 = [5 6; 7 8];
input3 = [9 10; 11 12];
expected = [15 18; 21 24];

inputs = {input1, input2, input3};
output = L.evaluate(inputs);

assert(isequal(output, expected), 'AdditionLayer evaluate failed for three inputs');

%% Test 4: AdditionLayer evaluate with image data
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});

% Create two 3D image arrays
IM1(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1];
IM1(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1];
IM1(:,:,3) = [1 1 1 1; 1 1 0 1; 0 1 1 0; 1 0 1 0];

IM2(:,:,1) = [0 1 0 0; 1 0 0 0; 0 1 0 1; 0 0 0 0];
IM2(:,:,2) = [1 0 1 1; 0 1 1 0; 1 0 0 1; 1 1 1 0];
IM2(:,:,3) = [0 0 0 0; 0 0 1 0; 1 0 0 1; 0 1 0 1];

inputs = {IM1, IM2};
output = L.evaluate(inputs);

expected(:,:,1) = IM1(:,:,1) + IM2(:,:,1);
expected(:,:,2) = IM1(:,:,2) + IM2(:,:,2);
expected(:,:,3) = IM1(:,:,3) + IM2(:,:,3);

assert(isequal(output, expected), 'AdditionLayer evaluate failed for image data');

%% Test 5: AdditionLayer reach with ImageStar inputs
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});

% Create first ImageStar
IM1(:,:,1) = [1 1 0 1; 0 0 1 1; 1 0 1 0; 1 1 1 1];
IM1(:,:,2) = [0 1 0 0; 1 0 0 1; 0 1 1 0; 0 0 0 1];

LB1(:,:,1) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB1(:,:,2) = [-0.1 -0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

UB1(:,:,1) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB1(:,:,2) = [0.1 0.1 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

image_star1 = ImageStar(IM1, LB1, UB1);

% Create second ImageStar
IM2(:,:,1) = [0 1 0 0; 1 0 0 0; 0 1 0 1; 0 0 0 0];
IM2(:,:,2) = [1 0 1 1; 0 1 1 0; 1 0 0 1; 1 1 1 0];

LB2(:,:,1) = [-0.05 -0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
LB2(:,:,2) = [-0.05 -0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

UB2(:,:,1) = [0.05 0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];
UB2(:,:,2) = [0.05 0.05 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

image_star2 = ImageStar(IM2, LB2, UB2);

% Test reach with approx-star
inputs = {image_star1, image_star2};
output_star = L.reach(inputs, 'approx-star');

assert(isa(output_star, 'ImageStar'), 'reach should return ImageStar');

%% Test 7: AdditionLayer toGPU
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});
L_gpu = L.toGPU();
assert(isa(L_gpu, 'AdditionLayer'));

%% Test 8: AdditionLayer changeParamsPrecision
L = AdditionLayer('add_layer', 2, 1, {'in1', 'in2'}, {'out'});
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'AdditionLayer'));

%% Test 9: Visual Soundness Verification with Star Sets
fprintf('\n=== Test 9: Visual Soundness Verification ===\n');

L9 = AdditionLayer('add_layer_visual', 2, 1, {'in1', 'in2'}, {'out'});

% Create two input Star sets with perturbations
input_dim = 4;
center1 = randn(input_dim, 1);
center2 = randn(input_dim, 1);
eps9 = 0.1;

lb1 = center1 - eps9;
ub1 = center1 + eps9;
S_in1 = Star(lb1, ub1);

lb2 = center2 - eps9;
ub2 = center2 + eps9;
S_in2 = Star(lb2, ub2);

% Compute reachable set
inputs9 = {S_in1, S_in2};
S_out9 = L9.reach(inputs9, 'approx-star');
[out_lb, out_ub] = S_out9.getRanges();
out_lb = out_lb(:);
out_ub = out_ub(:);

% Sample points and evaluate
n_samples = 50;
samples_out = zeros(input_dim, n_samples);
for i = 1:n_samples
    x1_sample = lb1 + (ub1 - lb1) .* rand(input_dim, 1);
    x2_sample = lb2 + (ub2 - lb2) .* rand(input_dim, 1);
    samples_out(:, i) = L9.evaluate({x1_sample, x2_sample});
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
    fill([d-0.4, d+0.4, d+0.4, d-0.4], ...
         [out_lb(d), out_lb(d), out_ub(d), out_ub(d)], ...
         colors(d,:), 'FaceAlpha', 0.3, 'EdgeColor', colors(d,:), 'LineWidth', 1.5);
    x_jitter = d + (rand(1, n_samples) - 0.5) * 0.6;
    scatter(x_jitter, samples_out(d, :), 30, colors(d,:), 'filled', 'MarkerFaceAlpha', 0.7);
end

% Status display
if all_contained
    title_str = sprintf('AdditionLayer Soundness: SOUND (%d samples)', n_samples);
    title(title_str, 'Color', [0 0.6 0], 'FontWeight', 'bold', 'FontSize', 12);
else
    title_str = sprintf('AdditionLayer Soundness: VIOLATION');
    title(title_str, 'Color', [0.8 0 0], 'FontWeight', 'bold', 'FontSize', 12);
end

xlabel('Output Dimension', 'FontSize', 11);
ylabel('Value', 'FontSize', 11);
xlim([0.5, input_dim + 0.5]);
xticks(1:input_dim);
grid on;
hold off;

% Save figure
try
    save_test_figure(fig, 'test_AdditionLayer', 'soundness_containment', 9, 'subdir', 'nn/layers/AdditionLayer');
catch
    fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'results', 'tests', 'figures', 'nn', 'layers', 'AdditionLayer');
    if ~exist(fig_dir, 'dir')
        mkdir(fig_dir);
    end
    saveas(fig, fullfile(fig_dir, 'test_AdditionLayer_soundness_containment_9.png'));
    fprintf('  Figure saved to: %s\n', fig_dir);
end
if isvalid(fig)
    close(fig);
end

assert(all_contained, 'Visual soundness verification failed');
fprintf('Test 9 PASSED: Visual soundness verification\n');

%% Summary
fprintf('\n=== All AdditionLayer Tests PASSED ===\n');
fprintf('Total: 9 tests\n');
