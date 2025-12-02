%% Verify MNIST ViT Models - ReLU and GELU Versions
% This script performs robustness verification on both ViT models:
% - mnist_vit_relu_model.mat (ReLU activation)
% - mnist_vit_gelu_model.mat (GELU activation)
%
% Author: NNV Team
% Date: December 2025

fprintf('=== Verifying MNIST ViT Models (ReLU and GELU) ===\n\n');

%% Configuration
epsilons = [0.01, 0.02, 0.05];  % L-inf perturbation levels
num_test = 10;  % Number of test images

%% Load test data
fprintf('[Setup] Loading MNIST test data...\n');
[XTest, YTest] = digitTest4DArrayData;
XTest = double(XTest);

%% Load Networks
fprintf('[Setup] Loading trained networks...\n');

models = struct();

% ReLU model
relu_file = 'mnist_vit_relu_model.mat';
if exist(relu_file, 'file')
    relu_data = load(relu_file);
    models(1).name = 'ReLU ViT';
    models(1).net = relu_data.net;
    models(1).accuracy = relu_data.accuracy;
    fprintf('  Loaded ReLU ViT (Accuracy: %.2f%%)\n', relu_data.accuracy * 100);
else
    error('ReLU model not found. Run train_mnist_vit_attention.m first.');
end

% GELU model
gelu_file = 'mnist_vit_gelu_model.mat';
if exist(gelu_file, 'file')
    gelu_data = load(gelu_file);
    models(2).name = 'GELU ViT';
    models(2).net = gelu_data.net;
    models(2).accuracy = gelu_data.accuracy;
    fprintf('  Loaded GELU ViT (Accuracy: %.2f%%)\n', gelu_data.accuracy * 100);
else
    error('GELU model not found. Run train_mnist_vit_gelu.m first.');
end

%% Convert to NNV
fprintf('\n[Convert] Converting networks to NNV...\n');

for m = 1:length(models)
    models(m).nnv = matlab2nnv(models(m).net);
    fprintf('  %s: %d layers\n', models(m).name, length(models(m).nnv.Layers));

    % Identify attention layer
    for i = 1:length(models(m).nnv.Layers)
        if isa(models(m).nnv.Layers{i}, 'MultiHeadAttentionLayer')
            fprintf('    MultiHeadAttentionLayer at index %d\n', i);
            break;
        end
    end
end

%% Select test images
fprintf('\n[Select] Selecting %d test images...\n', num_test);

rng(42);  % Reproducibility
test_indices = randperm(size(XTest, 4), num_test);

% Filter to correctly classified images (by ReLU model)
correct_indices = [];
for i = test_indices
    img = XTest(:,:,:,i);
    true_label = double(YTest(i));

    pred = classify(models(1).net, single(img));
    pred_class = double(pred) - 1;

    if pred_class == true_label
        correct_indices = [correct_indices, i];
    end
end
test_indices = correct_indices(1:min(num_test, length(correct_indices)));
fprintf('  Selected %d correctly classified images\n', length(test_indices));

%% Initialize results
results = struct();
for m = 1:length(models)
    results(m).name = models(m).name;
    results(m).verified = zeros(length(test_indices), length(epsilons));
    results(m).times = zeros(length(test_indices), length(epsilons));
    results(m).errors = zeros(length(test_indices), length(epsilons));
end

%% Reachability options
reachOptions = struct('reachMethod', 'approx-star');

%% Run verification
fprintf('\n[Verify] Running robustness verification...\n\n');

for img_idx = 1:length(test_indices)
    i = test_indices(img_idx);
    img = XTest(:,:,:,i);
    true_label = double(YTest(i));

    fprintf('Image %d/%d (digit %d):\n', img_idx, length(test_indices), true_label);

    for eps_idx = 1:length(epsilons)
        eps = epsilons(eps_idx);

        % Create perturbation set
        lb = max(img - eps, 0);
        ub = min(img + eps, 1);
        input_set = ImageStar(lb, ub);

        fprintf('  eps=%.2f: ', eps);

        for m = 1:length(models)
            try
                tic;
                output_set = models(m).nnv.reach(input_set, reachOptions);
                reach_time = toc;
                results(m).times(img_idx, eps_idx) = reach_time;

                % Get output bounds
                if iscell(output_set)
                    output_set = output_set{1};
                end

                if isa(output_set, 'ImageStar')
                    if ~isempty(output_set.im_lb) && ~isempty(output_set.im_ub)
                        lb_out = output_set.im_lb(:);
                        ub_out = output_set.im_ub(:);
                    else
                        [lb_img, ub_img] = output_set.estimateRanges();
                        lb_out = lb_img(:);
                        ub_out = ub_img(:);
                    end
                elseif isa(output_set, 'Star')
                    [lb_out, ub_out] = output_set.getRanges();
                else
                    lb_out = zeros(10, 1);
                    ub_out = zeros(10, 1);
                end

                % Check robustness: true class lb > max other ub
                true_class = true_label + 1;  % 1-indexed
                true_class_lb = lb_out(true_class);
                other_ub = ub_out;
                other_ub(true_class) = -inf;
                max_other_ub = max(other_ub);

                if true_class_lb > max_other_ub
                    results(m).verified(img_idx, eps_idx) = 1;
                    fprintf('%s:V ', models(m).name(1:4));
                else
                    fprintf('%s:X ', models(m).name(1:4));
                end

            catch ME
                results(m).errors(img_idx, eps_idx) = 1;
                fprintf('%s:E ', models(m).name(1:4));
            end
        end
        fprintf('| ');
    end
    fprintf('\n');
end

%% Summary
fprintf('\n============================================\n');
fprintf('           Verification Summary\n');
fprintf('============================================\n\n');

for m = 1:length(models)
    fprintf('--- %s ---\n', models(m).name);
    fprintf('Training Accuracy: %.2f%%\n', models(m).accuracy * 100);
    fprintf('Layers: %d\n\n', length(models(m).nnv.Layers));

    for eps_idx = 1:length(epsilons)
        eps = epsilons(eps_idx);
        n_verified = sum(results(m).verified(:, eps_idx));
        n_errors = sum(results(m).errors(:, eps_idx));
        n_total = length(test_indices);
        avg_time = mean(results(m).times(:, eps_idx));

        fprintf('  eps=%.2f: %d/%d verified (%.1f%%), %d errors, avg time: %.2fs\n', ...
            eps, n_verified, n_total, 100*n_verified/n_total, n_errors, avg_time);
    end
    fprintf('\n');
end

fprintf('============================================\n');
fprintf('Legend: V=verified robust, X=not verified, E=error\n');
fprintf('============================================\n');

%% Save results
save('verification_results_both_vit.mat', 'results', 'epsilons', 'test_indices', 'models');
fprintf('\nResults saved to verification_results_both_vit.mat\n');

%% Create comparison plot
figure('Position', [100 100 800 400]);

bar_data = zeros(length(epsilons), length(models));
for m = 1:length(models)
    for e = 1:length(epsilons)
        bar_data(e, m) = sum(results(m).verified(:, e)) / length(test_indices) * 100;
    end
end

bar(bar_data);
xlabel('Perturbation Epsilon');
ylabel('Verified Robust (%)');
title('MNIST ViT Robustness Verification: ReLU vs GELU');
xticklabels(arrayfun(@(e) sprintf('%.2f', e), epsilons, 'UniformOutput', false));
legend({models.name}, 'Location', 'northeast');
ylim([0 110]);
grid on;

saveas(gcf, 'plots/verification_comparison.png');
saveas(gcf, 'plots/verification_comparison.fig');
fprintf('Saved verification_comparison.png\n');
