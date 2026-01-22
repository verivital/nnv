%% Soundness Plots for MNIST ViT (ReLU and GELU versions)
% Visualizes output bounds with sampled points to demonstrate soundness
%
% Creates plots showing:
% - Output bounds (boxes) from reachability analysis
% - Sampled outputs (dots) contained within bounds
% - Class separation margins
%
% Author: NNV Team
% Date: December 2025

fprintf('=== MNIST ViT Soundness Plots ===\n\n');

%% Setup
example_dir = fileparts(mfilename('fullpath'));

% Create plots directory
plots_dir = fullfile(example_dir, 'plots');
if ~exist(plots_dir, 'dir')
    mkdir(plots_dir);
end

%% Load MNIST test data
fprintf('[Load] Loading MNIST test data...\n');
[XTest, YTest] = digitTest4DArrayData;
XTest = double(XTest);

%% Load Networks
fprintf('[Load] Loading trained networks...\n');

% Load ReLU model
relu_file = fullfile(example_dir, 'mnist_vit_relu_model.mat');
if exist(relu_file, 'file')
    relu_data = load(relu_file);
    relu_net = relu_data.net;
    fprintf('  Loaded ReLU ViT (Accuracy: %.2f%%)\n', relu_data.accuracy);
else
    error('ReLU model not found. Run train_mnist_vit_attention.m first.');
end

% Load GELU model
gelu_file = fullfile(example_dir, 'mnist_vit_gelu_model.mat');
if exist(gelu_file, 'file')
    gelu_data = load(gelu_file);
    gelu_net = gelu_data.net;
    fprintf('  Loaded GELU ViT (Accuracy: %.2f%%)\n', gelu_data.accuracy);
else
    error('GELU model not found. Run train_mnist_vit_gelu.m first.');
end

%% Convert to NNV
fprintf('[Convert] Converting networks to NNV...\n');

nnv_relu = matlab2nnv(relu_net);
fprintf('  ReLU NNV: %d layers\n', length(nnv_relu.Layers));

nnv_gelu = matlab2nnv(gelu_net);
fprintf('  GELU NNV: %d layers\n', length(nnv_gelu.Layers));

%% Configuration
num_classes = 10;
class_names = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
% Use smaller epsilons for tighter bounds visualization
% Note: Larger epsilons cause huge output ranges due to:
%   1. LayerNorm over-approximation
%   2. Large FC layers amplify uncertainty (~20x per layer)
epsilons = [1e-4, 5e-4, 1e-3, 2e-3];
n_samples_per_image = 50;  % Samples for soundness check
n_test_images = 3;         % Number of images to test

% Color map for classes
colors = lines(num_classes);

% Reach options
reachOptions = struct('reachMethod', 'approx-star');

%% Helper function to get bounds from reach output
function [out_lb, out_ub] = get_reach_bounds(R, num_classes)
    if iscell(R)
        out_lb = inf(num_classes, 1);
        out_ub = -inf(num_classes, 1);
        for r = 1:length(R)
            if isa(R{r}, 'ImageStar')
                if ~isempty(R{r}.im_lb)
                    lb_r = R{r}.im_lb(:);
                    ub_r = R{r}.im_ub(:);
                else
                    [lb_r, ub_r] = R{r}.estimateRanges();
                    lb_r = lb_r(:);
                    ub_r = ub_r(:);
                end
            else
                [lb_r, ub_r] = R{r}.getRanges();
                lb_r = lb_r(:);
                ub_r = ub_r(:);
            end
            out_lb = min(out_lb, lb_r);
            out_ub = max(out_ub, ub_r);
        end
    else
        if isa(R, 'ImageStar')
            if ~isempty(R.im_lb)
                out_lb = R.im_lb(:);
                out_ub = R.im_ub(:);
            else
                [out_lb, out_ub] = R.estimateRanges();
                out_lb = out_lb(:);
                out_ub = out_ub(:);
            end
        else
            [out_lb, out_ub] = R.getRanges();
            out_lb = out_lb(:);
            out_ub = out_ub(:);
        end
    end
end

%% Helper function to evaluate network on image - returns LOGITS (pre-softmax)
function y = eval_net(net, img)
    % Evaluate MATLAB network and return LOGITS (pre-softmax output)
    % NNV computes bounds on logits, so we must compare to logits
    % The classifier layer is the last FC layer before softmax
    if isa(net, 'DAGNetwork') || isa(net, 'SeriesNetwork')
        % Get the classifier layer (last FC before softmax)
        classifier_layer = 'classifier';  % Standard name for our ViT networks
        y = activations(net, single(img), classifier_layer);
        y = double(y(:));
    else
        % dlnetwork - get pre-softmax output
        img_dl = dlarray(single(img), 'SSCB');
        y_dl = predict(net, img_dl);
        y = double(extractdata(y_dl));
    end
end

%% Select test images (one from each correctly classified class)
fprintf('\n[Select] Selecting test images...\n');

test_indices = [];
classes_found = [];

for i = 1:size(XTest, 4)
    if length(test_indices) >= n_test_images
        break;
    end

    img = XTest(:,:,:,i);
    true_label = YTest(i);

    % Check if ReLU model classifies correctly (use original MATLAB net)
    y_relu = eval_net(relu_net, img);
    [~, pred_relu] = max(y_relu(:));
    pred_class_relu = pred_relu - 1;  % 0-indexed

    if pred_class_relu == double(true_label) && ~ismember(double(true_label), classes_found)
        test_indices = [test_indices, i];
        classes_found = [classes_found, double(true_label)];
        fprintf('  Image %d: digit %d (correctly classified)\n', i, double(true_label));
    end
end

fprintf('  Selected %d test images\n', length(test_indices));

%% Compute verification for each network and epsilon
fprintf('\n[Verify] Computing reachability for soundness plots...\n');

results_relu = struct();
results_gelu = struct();

for img_idx = 1:length(test_indices)
    i = test_indices(img_idx);
    img = XTest(:,:,:,i);
    true_label = double(YTest(i));

    fprintf('\n  Image %d (digit %d):\n', img_idx, true_label);

    for eps_idx = 1:length(epsilons)
        eps = epsilons(eps_idx);
        fprintf('    eps=%.3f: ', eps);

        % Create perturbation set
        lb = max(img - eps, 0);
        ub = min(img + eps, 1);
        input_set = ImageStar(lb, ub);

        % Sample points for soundness check
        sampled_relu = zeros(num_classes, n_samples_per_image);
        sampled_gelu = zeros(num_classes, n_samples_per_image);

        for s = 1:n_samples_per_image
            % Random sample within bounds
            x_sample = lb + (ub - lb) .* rand(size(img));

            % Evaluate networks (use original MATLAB nets for forward pass)
            y_relu = eval_net(relu_net, x_sample);
            y_gelu = eval_net(gelu_net, x_sample);

            sampled_relu(:, s) = y_relu(:);
            sampled_gelu(:, s) = y_gelu(:);
        end

        % Compute reachability for ReLU network
        try
            R_relu = nnv_relu.reach(input_set, reachOptions);
            [lb_relu, ub_relu] = get_reach_bounds(R_relu, num_classes);

            results_relu(img_idx, eps_idx).bounds.lb = lb_relu;
            results_relu(img_idx, eps_idx).bounds.ub = ub_relu;
            results_relu(img_idx, eps_idx).samples = sampled_relu;
            results_relu(img_idx, eps_idx).label = true_label;
            results_relu(img_idx, eps_idx).eps = eps;

            % Check soundness
            sound_relu = all(all(sampled_relu >= lb_relu - 1e-5)) && ...
                        all(all(sampled_relu <= ub_relu + 1e-5));
            results_relu(img_idx, eps_idx).sound = sound_relu;

            fprintf('ReLU %s, ', ifelse(sound_relu, 'SOUND', 'VIOLATION'));
        catch ME
            fprintf('ReLU ERROR, ');
            results_relu(img_idx, eps_idx).bounds.lb = [];
            results_relu(img_idx, eps_idx).bounds.ub = [];
            results_relu(img_idx, eps_idx).samples = sampled_relu;
            results_relu(img_idx, eps_idx).sound = false;
        end

        % Compute reachability for GELU network
        try
            R_gelu = nnv_gelu.reach(input_set, reachOptions);
            [lb_gelu, ub_gelu] = get_reach_bounds(R_gelu, num_classes);

            results_gelu(img_idx, eps_idx).bounds.lb = lb_gelu;
            results_gelu(img_idx, eps_idx).bounds.ub = ub_gelu;
            results_gelu(img_idx, eps_idx).samples = sampled_gelu;
            results_gelu(img_idx, eps_idx).label = true_label;
            results_gelu(img_idx, eps_idx).eps = eps;

            % Check soundness
            sound_gelu = all(all(sampled_gelu >= lb_gelu - 1e-5)) && ...
                        all(all(sampled_gelu <= ub_gelu + 1e-5));
            results_gelu(img_idx, eps_idx).sound = sound_gelu;

            fprintf('GELU %s\n', ifelse(sound_gelu, 'SOUND', 'VIOLATION'));
        catch ME
            fprintf('GELU ERROR\n');
            results_gelu(img_idx, eps_idx).bounds.lb = [];
            results_gelu(img_idx, eps_idx).bounds.ub = [];
            results_gelu(img_idx, eps_idx).samples = sampled_gelu;
            results_gelu(img_idx, eps_idx).sound = false;
        end
    end
end

%% Helper function for conditional string
function s = ifelse(cond, true_str, false_str)
    if cond
        s = true_str;
    else
        s = false_str;
    end
end

%% Plot 1: ReLU ViT Soundness - Output Bounds Containment
fprintf('\n[Plot 1] ReLU ViT - Output Bounds Containment\n');

figure('Position', [100 100 1400 900]);
sgtitle('ReLU ViT: Reachability Bounds Contain All Sampled Outputs', 'FontSize', 14, 'FontWeight', 'bold');

for eps_idx = 1:length(epsilons)
    subplot(2, 2, eps_idx);
    hold on;

    eps = epsilons(eps_idx);

    for img_idx = 1:min(n_test_images, length(test_indices))
        if isempty(results_relu(img_idx, eps_idx).bounds.lb)
            continue;
        end

        bounds = results_relu(img_idx, eps_idx).bounds;
        samples = results_relu(img_idx, eps_idx).samples;
        label = results_relu(img_idx, eps_idx).label;

        x_offset = (img_idx - 1) * 11;

        for c = 1:num_classes
            lb = bounds.lb(c);
            ub = bounds.ub(c);

            % Highlight true class
            if c == label + 1
                face_alpha = 0.5;
                edge_width = 2;
            else
                face_alpha = 0.2;
                edge_width = 1;
            end

            % Draw bound box
            fill([x_offset + c - 0.35, x_offset + c + 0.35, x_offset + c + 0.35, x_offset + c - 0.35], ...
                 [lb, lb, ub, ub], colors(c, :), 'FaceAlpha', face_alpha, 'EdgeColor', colors(c, :), 'LineWidth', edge_width);

            % Plot sample points (at least 20 visible samples with larger markers)
            x_jitter = x_offset + c + (rand(1, size(samples, 2)) - 0.5) * 0.5;
            scatter(x_jitter, samples(c, :), 30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
        end

        % Label
        text(x_offset + 5.5, min(ylim) + 0.1 * diff(ylim), sprintf('digit %d', label), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end

    xlabel('Output Class (grouped by image)');
    ylabel('Output Logit');

    n_sound = sum([results_relu(:, eps_idx).sound]);
    title(sprintf('\\epsilon = %.3f (Sound: %d/%d images)', eps, n_sound, n_test_images));
    grid on;
    hold off;
end

saveas(gcf, fullfile(plots_dir, 'relu_vit_soundness.png'));
saveas(gcf, fullfile(plots_dir, 'relu_vit_soundness.fig'));
fprintf('  Saved relu_vit_soundness.png\n');

%% Plot 2: GELU ViT Soundness - Output Bounds Containment
fprintf('\n[Plot 2] GELU ViT - Output Bounds Containment\n');

figure('Position', [100 100 1400 900]);
sgtitle('GELU ViT: Reachability Bounds Contain All Sampled Outputs', 'FontSize', 14, 'FontWeight', 'bold');

for eps_idx = 1:length(epsilons)
    subplot(2, 2, eps_idx);
    hold on;

    eps = epsilons(eps_idx);

    for img_idx = 1:min(n_test_images, length(test_indices))
        if isempty(results_gelu(img_idx, eps_idx).bounds.lb)
            continue;
        end

        bounds = results_gelu(img_idx, eps_idx).bounds;
        samples = results_gelu(img_idx, eps_idx).samples;
        label = results_gelu(img_idx, eps_idx).label;

        x_offset = (img_idx - 1) * 11;

        for c = 1:num_classes
            lb = bounds.lb(c);
            ub = bounds.ub(c);

            if c == label + 1
                face_alpha = 0.5;
                edge_width = 2;
            else
                face_alpha = 0.2;
                edge_width = 1;
            end

            fill([x_offset + c - 0.35, x_offset + c + 0.35, x_offset + c + 0.35, x_offset + c - 0.35], ...
                 [lb, lb, ub, ub], colors(c, :), 'FaceAlpha', face_alpha, 'EdgeColor', colors(c, :), 'LineWidth', edge_width);

            x_jitter = x_offset + c + (rand(1, size(samples, 2)) - 0.5) * 0.5;
            scatter(x_jitter, samples(c, :), 8, colors(c, :), 'filled', 'MarkerFaceAlpha', 0.5);
        end

        text(x_offset + 5.5, min(ylim) + 0.1 * diff(ylim), sprintf('digit %d', label), ...
            'HorizontalAlignment', 'center', 'FontSize', 9);
    end

    xlabel('Output Class (grouped by image)');
    ylabel('Output Logit');

    n_sound = sum([results_gelu(:, eps_idx).sound]);
    title(sprintf('\\epsilon = %.3f (Sound: %d/%d images)', eps, n_sound, n_test_images));
    grid on;
    hold off;
end

saveas(gcf, fullfile(plots_dir, 'gelu_vit_soundness.png'));
saveas(gcf, fullfile(plots_dir, 'gelu_vit_soundness.fig'));
fprintf('  Saved gelu_vit_soundness.png\n');

%% Plot 3: Detailed Single Image - ReLU vs GELU Comparison
fprintf('\n[Plot 3] Single Image Detail - ReLU vs GELU Comparison\n');

figure('Position', [100 100 1200 600]);
sgtitle('ReLU vs GELU ViT: Detailed Bound Containment (First Image)', 'FontSize', 14, 'FontWeight', 'bold');

for eps_idx = 1:length(epsilons)
    % ReLU subplot
    subplot(2, length(epsilons), eps_idx);
    hold on;

    if ~isempty(results_relu(1, eps_idx).bounds.lb)
        bounds = results_relu(1, eps_idx).bounds;
        samples = results_relu(1, eps_idx).samples;
        label = results_relu(1, eps_idx).label;

        for c = 1:num_classes
            lb = bounds.lb(c);
            ub = bounds.ub(c);

            if c == label + 1
                face_alpha = 0.4;
                edge_color = [0.8 0.2 0.2];
            else
                face_alpha = 0.15;
                edge_color = colors(c, :);
            end

            fill([c - 0.4, c + 0.4, c + 0.4, c - 0.4], ...
                 [lb, lb, ub, ub], colors(c, :), 'FaceAlpha', face_alpha, 'EdgeColor', edge_color, 'LineWidth', 1.5);

            x_jitter = c + (rand(1, size(samples, 2)) - 0.5) * 0.6;
            scatter(x_jitter, samples(c, :), 20, colors(c, :), 'filled', 'MarkerFaceAlpha', 0.6);
        end

        status = ifelse(results_relu(1, eps_idx).sound, 'SOUND', 'VIOLATION');
        status_color = ifelse(results_relu(1, eps_idx).sound, [0 0.6 0], [0.8 0 0]);
    else
        status = 'ERROR';
        status_color = [0.5 0.5 0.5];
    end

    xlabel('Class');
    ylabel('Logit');
    title(sprintf('ReLU | \\epsilon=%.3f | %s', epsilons(eps_idx), status), 'Color', status_color);
    xticks(1:num_classes);
    xticklabels(class_names);
    grid on;
    hold off;

    % GELU subplot
    subplot(2, length(epsilons), eps_idx + length(epsilons));
    hold on;

    if ~isempty(results_gelu(1, eps_idx).bounds.lb)
        bounds = results_gelu(1, eps_idx).bounds;
        samples = results_gelu(1, eps_idx).samples;
        label = results_gelu(1, eps_idx).label;

        for c = 1:num_classes
            lb = bounds.lb(c);
            ub = bounds.ub(c);

            if c == label + 1
                face_alpha = 0.4;
                edge_color = [0.8 0.2 0.2];
            else
                face_alpha = 0.15;
                edge_color = colors(c, :);
            end

            fill([c - 0.4, c + 0.4, c + 0.4, c - 0.4], ...
                 [lb, lb, ub, ub], colors(c, :), 'FaceAlpha', face_alpha, 'EdgeColor', edge_color, 'LineWidth', 1.5);

            x_jitter = c + (rand(1, size(samples, 2)) - 0.5) * 0.6;
            scatter(x_jitter, samples(c, :), 20, colors(c, :), 'filled', 'MarkerFaceAlpha', 0.6);
        end

        status = ifelse(results_gelu(1, eps_idx).sound, 'SOUND', 'VIOLATION');
        status_color = ifelse(results_gelu(1, eps_idx).sound, [0 0.6 0], [0.8 0 0]);
    else
        status = 'ERROR';
        status_color = [0.5 0.5 0.5];
    end

    xlabel('Class');
    ylabel('Logit');
    title(sprintf('GELU | \\epsilon=%.3f | %s', epsilons(eps_idx), status), 'Color', status_color);
    xticks(1:num_classes);
    xticklabels(class_names);
    grid on;
    hold off;
end

saveas(gcf, fullfile(plots_dir, 'relu_vs_gelu_detail.png'));
saveas(gcf, fullfile(plots_dir, 'relu_vs_gelu_detail.fig'));
fprintf('  Saved relu_vs_gelu_detail.png\n');

%% Plot 4: Bound Width Comparison
fprintf('\n[Plot 4] Bound Width Comparison\n');

figure('Position', [100 100 800 400]);

relu_widths = zeros(1, length(epsilons));
gelu_widths = zeros(1, length(epsilons));

for eps_idx = 1:length(epsilons)
    if ~isempty(results_relu(1, eps_idx).bounds.lb)
        relu_widths(eps_idx) = mean(results_relu(1, eps_idx).bounds.ub - results_relu(1, eps_idx).bounds.lb);
    end
    if ~isempty(results_gelu(1, eps_idx).bounds.lb)
        gelu_widths(eps_idx) = mean(results_gelu(1, eps_idx).bounds.ub - results_gelu(1, eps_idx).bounds.lb);
    end
end

bar_width = 0.35;
x = 1:length(epsilons);

bar(x - bar_width/2, relu_widths, bar_width, 'FaceColor', [0.2 0.4 0.8], 'EdgeColor', 'none');
hold on;
bar(x + bar_width/2, gelu_widths, bar_width, 'FaceColor', [0.8 0.4 0.2], 'EdgeColor', 'none');

xlabel('Perturbation Epsilon');
ylabel('Average Bound Width');
title('Output Bound Tightness: ReLU vs GELU ViT');
xticks(x);
xticklabels(arrayfun(@(e) sprintf('%.3f', e), epsilons, 'UniformOutput', false));
legend({'ReLU ViT', 'GELU ViT'}, 'Location', 'northwest');
grid on;
hold off;

saveas(gcf, fullfile(plots_dir, 'bound_width_comparison.png'));
fprintf('  Saved bound_width_comparison.png\n');

%% Save Results
save(fullfile(plots_dir, 'soundness_results.mat'), 'results_relu', 'results_gelu', 'epsilons', 'test_indices');
fprintf('\n  Saved soundness_results.mat\n');

%% Summary
fprintf('\n=== Soundness Plot Generation Complete ===\n');
fprintf('Generated plots in: %s\n', plots_dir);
fprintf('\nPlots:\n');
fprintf('  1. relu_vit_soundness.png - ReLU ViT bounds containment\n');
fprintf('  2. gelu_vit_soundness.png - GELU ViT bounds containment\n');
fprintf('  3. relu_vs_gelu_detail.png - Detailed comparison\n');
fprintf('  4. bound_width_comparison.png - Bound tightness comparison\n');
