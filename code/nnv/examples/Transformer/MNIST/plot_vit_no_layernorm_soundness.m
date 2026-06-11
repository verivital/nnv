%% Soundness Plots for MNIST ViT without LayerNormalization
% Compares bound tightness between networks with and without LayerNorm
%
% Author: NNV Team
% Date: December 2025

fprintf('=== MNIST ViT No-LayerNorm Soundness Plots ===\n\n');

%% Setup
example_dir = fileparts(mfilename('fullpath'));
plots_dir = fullfile(example_dir, 'plots');
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end

%% Load Networks
fprintf('[Load] Loading networks...\n');

% Load no-LayerNorm model
no_ln_file = fullfile(example_dir, 'mnist_vit_no_layernorm_model.mat');
if exist(no_ln_file, 'file')
    no_ln_data = load(no_ln_file);
    no_ln_net = no_ln_data.net;
    fprintf('  Loaded No-LayerNorm ViT (Accuracy: %.2f%%)\n', no_ln_data.accuracy);
else
    error('No-LayerNorm model not found. Run train_mnist_vit_no_layernorm.m first.');
end

% Load original ReLU model (with LayerNorm) for comparison
relu_file = fullfile(example_dir, 'mnist_vit_relu_model.mat');
if exist(relu_file, 'file')
    relu_data = load(relu_file);
    relu_net = relu_data.net;
    fprintf('  Loaded ReLU ViT with LayerNorm (Accuracy: %.2f%%)\n', relu_data.accuracy);
    has_relu = true;
else
    fprintf('  ReLU model with LayerNorm not found (skipping comparison)\n');
    has_relu = false;
end

%% Convert to NNV
fprintf('[Convert] Converting to NNV...\n');

nnv_no_ln = matlab2nnv(no_ln_net);
fprintf('  No-LayerNorm NNV: %d layers\n', length(nnv_no_ln.Layers));

if has_relu
    nnv_relu = matlab2nnv(relu_net);
    fprintf('  With-LayerNorm NNV: %d layers\n', length(nnv_relu.Layers));
end

%% Load test data
fprintf('[Data] Loading test data...\n');
[XTest, YTest] = digitTest4DArrayData;
XTest = double(XTest);

%% Configuration
num_classes = 10;
epsilons = [1e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2];
n_samples = 50;
n_test_images = 1;

reachOptions = struct('reachMethod', 'approx-star');

%% Select test image
img_idx = 1;
img = XTest(:,:,:,img_idx);
true_label = double(YTest(img_idx));
fprintf('  Test image: digit %d\n', true_label);

%% Compute bounds for each epsilon
fprintf('\n[Verify] Computing reachability bounds...\n');

no_ln_ranges = zeros(length(epsilons), 1);
relu_ranges = zeros(length(epsilons), 1);
no_ln_sound = true(length(epsilons), 1);
relu_sound = true(length(epsilons), 1);

bounds_no_ln = cell(length(epsilons), 1);
bounds_relu = cell(length(epsilons), 1);
samples_all = cell(length(epsilons), 1);

for eps_idx = 1:length(epsilons)
    eps = epsilons(eps_idx);
    fprintf('  eps=%.0e: ', eps);

    lb = max(img - eps, 0);
    ub = min(img + eps, 1);
    input_set = ImageStar(lb, ub);

    % No-LayerNorm network
    try
        R_no_ln = nnv_no_ln.reach(input_set, reachOptions);
        if iscell(R_no_ln), R_no_ln = R_no_ln{1}; end
        if isa(R_no_ln, 'ImageStar')
            [out_lb, out_ub] = R_no_ln.estimateRanges();
        else
            [out_lb, out_ub] = R_no_ln.getRanges();
        end
        out_lb = out_lb(:); out_ub = out_ub(:);
        bounds_no_ln{eps_idx} = struct('lb', out_lb, 'ub', out_ub);
        no_ln_ranges(eps_idx) = max(out_ub - out_lb);
        fprintf('No-LN range=%.2f, ', no_ln_ranges(eps_idx));
    catch ME
        fprintf('No-LN FAILED, ');
        bounds_no_ln{eps_idx} = struct('lb', [], 'ub', []);
    end

    % With-LayerNorm network
    if has_relu
        try
            R_relu = nnv_relu.reach(input_set, reachOptions);
            if iscell(R_relu), R_relu = R_relu{1}; end
            if isa(R_relu, 'ImageStar')
                [out_lb, out_ub] = R_relu.estimateRanges();
            else
                [out_lb, out_ub] = R_relu.getRanges();
            end
            out_lb = out_lb(:); out_ub = out_ub(:);
            bounds_relu{eps_idx} = struct('lb', out_lb, 'ub', out_ub);
            relu_ranges(eps_idx) = max(out_ub - out_lb);
            fprintf('With-LN range=%.2f', relu_ranges(eps_idx));
        catch ME
            fprintf('With-LN FAILED');
            bounds_relu{eps_idx} = struct('lb', [], 'ub', []);
        end
    end

    % Sample outputs for soundness check - use LOGITS not softmax probabilities!
    % NNV computes bounds on pre-softmax logits, so we must compare to logits
    sampled = zeros(num_classes, n_samples);
    for s = 1:n_samples
        x_sample = lb + (ub - lb) .* rand(size(img));
        % Get LOGITS from classifier layer (pre-softmax)
        logits = activations(no_ln_net, single(x_sample), 'classifier');
        sampled(:, s) = double(logits(:));
    end
    samples_all{eps_idx} = sampled;

    % Check soundness
    if ~isempty(bounds_no_ln{eps_idx}.lb)
        for s = 1:n_samples
            if any(sampled(:,s) < bounds_no_ln{eps_idx}.lb - 1e-5) || ...
               any(sampled(:,s) > bounds_no_ln{eps_idx}.ub + 1e-5)
                no_ln_sound(eps_idx) = false;
                break;
            end
        end
    end

    fprintf('\n');
end

%% Plot 1: Bound width comparison
fprintf('\n[Plot] Creating comparison plots...\n');

figure('Position', [100, 100, 1200, 500]);

subplot(1, 2, 1);
if has_relu
    bar_data = [no_ln_ranges, relu_ranges];
    b = bar(bar_data);
    b(1).FaceColor = [0.2, 0.6, 0.2];
    b(2).FaceColor = [0.8, 0.2, 0.2];
    legend('Without LayerNorm', 'With LayerNorm', 'Location', 'northwest');
else
    bar(no_ln_ranges, 'FaceColor', [0.2, 0.6, 0.2]);
    legend('Without LayerNorm', 'Location', 'northwest');
end
xlabel('Epsilon Index');
ylabel('Max Output Range');
title('Output Bound Width: LayerNorm Comparison');
xticks(1:length(epsilons));
xticklabels(arrayfun(@(e) sprintf('%.0e', e), epsilons, 'UniformOutput', false));
grid on;

% Log scale comparison
subplot(1, 2, 2);
if has_relu
    semilogy(epsilons, no_ln_ranges, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    semilogy(epsilons, relu_ranges, 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
    legend('Without LayerNorm', 'With LayerNorm', 'Location', 'northwest');
else
    semilogy(epsilons, no_ln_ranges, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
    legend('Without LayerNorm', 'Location', 'northwest');
end
xlabel('Epsilon');
ylabel('Max Output Range (log scale)');
title('Bound Width vs Epsilon (Log Scale)');
grid on;

saveas(gcf, fullfile(plots_dir, 'layernorm_comparison.png'));
fprintf('  Saved: layernorm_comparison.png\n');

%% Plot 2: Soundness visualization for No-LayerNorm
figure('Position', [100, 100, 1400, 400]);

for eps_idx = 1:min(4, length(epsilons))
    subplot(1, 4, eps_idx);

    eps = epsilons(eps_idx);
    bounds = bounds_no_ln{eps_idx};
    samples = samples_all{eps_idx};

    if isempty(bounds.lb)
        text(0.5, 0.5, 'Failed', 'HorizontalAlignment', 'center');
        continue;
    end

    % Plot bounds as boxes
    for c = 1:num_classes
        x_pos = c - 0.3;
        width = 0.6;
        height = bounds.ub(c) - bounds.lb(c);
        rectangle('Position', [x_pos, bounds.lb(c), width, height], ...
            'FaceColor', [0.7, 0.9, 0.7], 'EdgeColor', [0.2, 0.6, 0.2], 'LineWidth', 1.5);
    end

    hold on;

    % Plot samples with jitter for visibility (20+ samples clearly visible)
    % Use larger markers, higher alpha, and jitter to show all samples
    sample_colors = lines(min(n_samples, 30));  % Different colors for samples
    for s = 1:n_samples
        % Add horizontal jitter to spread samples within each class
        jitter = (rand(1, num_classes) - 0.5) * 0.4;
        x_positions = (1:num_classes) + jitter;

        % Use larger markers with higher visibility
        if s <= 20
            % First 20 samples are most visible
            scatter(x_positions, samples(:, s), 40, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
        else
            % Additional samples slightly smaller
            scatter(x_positions, samples(:, s), 25, 'b', 'filled', 'MarkerFaceAlpha', 0.4);
        end
    end

    % Highlight true class
    xline(true_label + 1, 'k--', 'LineWidth', 2);

    xlabel('Class');
    ylabel('Output Value');
    xlim([0.5, num_classes + 0.5]);
    xticks(1:num_classes);
    xticklabels(0:9);

    status = 'SOUND';
    if ~no_ln_sound(eps_idx)
        status = 'VIOLATION';
    end
    title(sprintf('\\epsilon = %.0e | %s', eps, status));
    grid on;
end

sgtitle('No-LayerNorm ViT: Output Bounds and Samples', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, fullfile(plots_dir, 'no_layernorm_soundness.png'));
fprintf('  Saved: no_layernorm_soundness.png\n');

%% Plot 3: Improvement ratio
if has_relu
    figure('Position', [100, 100, 600, 400]);

    improvement = relu_ranges ./ no_ln_ranges;
    bar(improvement, 'FaceColor', [0.3, 0.5, 0.8]);
    xlabel('Epsilon');
    ylabel('Improvement Ratio (With-LN / Without-LN)');
    title('Bound Tightness Improvement by Removing LayerNorm');
    xticks(1:length(epsilons));
    xticklabels(arrayfun(@(e) sprintf('%.0e', e), epsilons, 'UniformOutput', false));
    grid on;

    % Add ratio labels
    for i = 1:length(improvement)
        text(i, improvement(i) + 5, sprintf('%.0fx', improvement(i)), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end

    saveas(gcf, fullfile(plots_dir, 'layernorm_improvement.png'));
    fprintf('  Saved: layernorm_improvement.png\n');
end

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('No-LayerNorm ViT:\n');
fprintf('  Accuracy: %.2f%%\n', no_ln_data.accuracy);
fprintf('  Output ranges by epsilon:\n');
for eps_idx = 1:length(epsilons)
    fprintf('    eps=%.0e: range=%.2f, sound=%d\n', ...
        epsilons(eps_idx), no_ln_ranges(eps_idx), no_ln_sound(eps_idx));
end

if has_relu
    fprintf('\nWith-LayerNorm ViT (for comparison):\n');
    fprintf('  Accuracy: %.2f%%\n', relu_data.accuracy);
    fprintf('  Output ranges by epsilon:\n');
    for eps_idx = 1:length(epsilons)
        fprintf('    eps=%.0e: range=%.2f (%.0fx larger)\n', ...
            epsilons(eps_idx), relu_ranges(eps_idx), relu_ranges(eps_idx)/no_ln_ranges(eps_idx));
    end
end

fprintf('\n=== Plots saved to: %s ===\n', plots_dir);
