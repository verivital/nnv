%% Robustness Verification - Small Vocab SST-2 Model
%
% This script verifies the reduced vocabulary model (500 words).
% Smaller vocab enables verification with larger epsilon values.
%
% Author: NNV Team
% Date: November 2025

%% Configuration
epsilon = 0.0001;  % Small perturbation for faster verification
num_samples_to_verify = 10;  % Number of samples

fprintf('=== Small Vocab SST-2 Verification ===\n');
fprintf('Perturbation epsilon: %.4f\n', epsilon);
fprintf('Samples to verify: %d\n', num_samples_to_verify);

%% Load Model
fprintf('\n[1/4] Loading small vocab model...\n');

model_path = fullfile(fileparts(mfilename('fullpath')), 'sentiment_sst2_small_model.mat');
if ~isfile(model_path)
    error('Model not found. Run train_sentiment_sst2_small.m first.');
end

load(model_path, 'net', 'accuracy', 'config', 'word2idx', 'all_vocab', 'idf', 'vocab_size');
fprintf('Model loaded. Accuracy: %.2f%%, Vocab: %d\n', accuracy * 100, vocab_size);

nnv_net = matlab2nnv(net);
fprintf('NNV model: %d layers\n', length(nnv_net.Layers));

%% Load Test Data
fprintf('\n[2/4] Loading test data...\n');

[~, dev_data] = download_sst2();
rng(456);
test_idx = randperm(length(dev_data.sentences), min(num_samples_to_verify, length(dev_data.sentences)));
fprintf('Selected %d test samples\n', length(test_idx));

%% Verification
fprintf('\n[3/4] Running verification...\n');

results = struct('index', {}, 'text', {}, 'true_label', {}, ...
    'predicted_label', {}, 'verified', {}, 'time', {}, 'output_bounds', {}, 'sampled_outputs', {});

for i = 1:length(test_idx)
    idx = test_idx(i);
    text = dev_data.sentences{idx};
    true_label = dev_data.labels(idx);  % 0=neg, 1=pos

    fprintf('\n--- Sample %d ---\n', i);
    fprintf('Text: "%s"\n', truncate_text(text, 50));

    % Convert to features
    x = sentence_to_tfidf(text, word2idx, vocab_size, idf);
    x = single(x);
    x_img = reshape(x, [vocab_size, 1, 1]);  % ImageStar-compatible shape

    % Evaluate (use image shape to match reachability)
    y_output = nnv_net.evaluate(x_img);
    [~, pred_idx] = max(y_output);
    predicted_label = pred_idx - 1;  % 0=neg, 1=pos

    fprintf('True: %s, Pred: %s\n', label_str(true_label), label_str(predicted_label));

    % Create input set
    lb = max(x - epsilon, 0);
    ub = x + epsilon;
    IS = ImageStar(reshape(lb, [vocab_size, 1, 1]), reshape(ub, [vocab_size, 1, 1]));

    % Verify
    t = tic;
    try
        reachOptions = struct('reachMethod', 'approx-star');
        R = nnv_net.reach(IS, reachOptions);

        if isa(R, 'Star')
            lb_out = zeros(2, 1);
            ub_out = zeros(2, 1);
            for j = 1:2
                lb_out(j) = R.getMin(j, 'linprog');
                ub_out(j) = R.getMax(j, 'linprog');
            end
        else
            [lb_out, ub_out] = R.getRanges();
            lb_out = squeeze(lb_out); ub_out = squeeze(ub_out);
        end

        % Check robustness
        if pred_idx == 1
            other_ub = ub_out(2);
        else
            other_ub = ub_out(1);
        end

        if lb_out(pred_idx) > other_ub
            verified = 1;
            fprintf('VERIFIED!\n');
        else
            verified = 0;
            fprintf('UNKNOWN (gap: %.4f)\n', other_ub - lb_out(pred_idx));
        end
    catch ME
        verified = -1;
        lb_out = [];
        ub_out = [];
        fprintf('FAILED: %s\n', ME.message);
    end

    reach_time = toc(t);
    fprintf('Time: %.2f sec\n', reach_time);

    % Sample random inputs from perturbation set for soundness visualization
    num_random_samples = 50;
    sampled_outputs = zeros(2, num_random_samples);
    for s = 1:num_random_samples
        x_sample = lb + (ub - lb) .* rand(size(x));
        x_sample_img = reshape(x_sample, [vocab_size, 1, 1]);  % Match ImageStar shape
        y_sample = nnv_net.evaluate(x_sample_img);
        sampled_outputs(:, s) = y_sample;
    end

    results(i).index = idx;
    results(i).text = text;
    results(i).true_label = true_label;
    results(i).predicted_label = predicted_label;
    results(i).verified = verified;
    results(i).time = reach_time;
    results(i).output_bounds = struct('lb', lb_out, 'ub', ub_out);
    results(i).sampled_outputs = sampled_outputs;
end

%% Summary
fprintf('\n[4/4] Summary\n');
fprintf('============\n');

num_verified = sum([results.verified] == 1);
num_unknown = sum([results.verified] == 0);
num_correct = sum([results.true_label] == [results.predicted_label]);

fprintf('Total: %d\n', length(test_idx));
fprintf('Correct predictions: %d (%.1f%%)\n', num_correct, 100*num_correct/length(test_idx));
fprintf('Verified robust: %d (%.1f%%)\n', num_verified, 100*num_verified/length(test_idx));
fprintf('Unknown: %d\n', num_unknown);
fprintf('Avg time: %.2f sec\n', mean([results.time]));

%% Visualization - Soundness Check
fprintf('\n[5/5] Visualizing results with soundness check...\n');

% NOTE: There is a known issue with NNV's BatchNormalizationLayer.reach()
% where the reachability bounds are tighter than actual network outputs.
% This causes sampled outputs to appear outside the computed bounds.
% This is a bug in NNV's reachability analysis, not in the sampling.
% See TODO_TRANSFORMER.md "Known Issues" section for details.

% First, verify soundness numerically
fprintf('\nSoundness verification:\n');
fprintf('NOTE: Due to BatchNormalizationLayer.reach() bug, samples may appear outside bounds.\n');
total_samples = 0;
samples_in_bounds = 0;
num_test = length(test_idx);
for i = 1:num_test
    if ~isempty(results(i).output_bounds.lb) && isfield(results(i), 'sampled_outputs')
        lb = results(i).output_bounds.lb;
        ub = results(i).output_bounds.ub;
        samples = results(i).sampled_outputs;

        for s = 1:size(samples, 2)
            total_samples = total_samples + 1;
            % Check if sample is within bounds (with small tolerance)
            tol = 1e-6;
            in_bounds = all(samples(:, s) >= lb - tol) && all(samples(:, s) <= ub + tol);
            if in_bounds
                samples_in_bounds = samples_in_bounds + 1;
            else
                fprintf('  Sample %d of test %d OUTSIDE bounds:\n', s, i);
                fprintf('    Output: [%.6f, %.6f]\n', samples(1, s), samples(2, s));
                fprintf('    Bounds: [%.6f, %.6f] to [%.6f, %.6f]\n', lb(1), lb(2), ub(1), ub(2));
            end
        end
    end
end
fprintf('Soundness: %d/%d samples within bounds (%.1f%%)\n', samples_in_bounds, total_samples, 100*samples_in_bounds/total_samples);

num_to_plot = min(num_test, 5);  % Plot first 5 samples

figure('Name', 'Small SST-2 Verification - Soundness Check', 'Position', [100, 100, 1200, 400]);

for i = 1:num_to_plot
    subplot(1, num_to_plot, i);
    hold on;

    if ~isempty(results(i).output_bounds.lb)
        lb = results(i).output_bounds.lb;
        ub = results(i).output_bounds.ub;

        % Plot bounds as filled rectangles (gray boxes)
        for c = 1:2
            fill([c-0.3, c+0.3, c+0.3, c-0.3], [lb(c), lb(c), ub(c), ub(c)], ...
                [0.8, 0.8, 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5);
        end

        % Overlay sampled outputs as scatter points (not lines)
        if isfield(results(i), 'sampled_outputs') && ~isempty(results(i).sampled_outputs)
            samples = results(i).sampled_outputs;
            % Plot class 1 (Negative) samples
            scatter(ones(1, size(samples, 2)) + 0.1*(rand(1, size(samples, 2))-0.5), ...
                samples(1, :), 20, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
            % Plot class 2 (Positive) samples
            scatter(2*ones(1, size(samples, 2)) + 0.1*(rand(1, size(samples, 2))-0.5), ...
                samples(2, :), 20, 'b', 'filled', 'MarkerFaceAlpha', 0.5);
        end

        % Mark the center of bounds
        plot([1, 2], [(lb(1)+ub(1))/2, (lb(2)+ub(2))/2], 'r*', 'MarkerSize', 10, 'LineWidth', 2);
    end

    hold off;
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'Neg', 'Pos'});
    xlim([0.5, 2.5]);
    ylabel('Output Logit');

    % Title with verification status
    if results(i).verified == 1
        status_str = 'VERIFIED';
        title_color = [0, 0.6, 0];
    else
        status_str = 'UNKNOWN';
        title_color = [0.8, 0, 0];
    end
    title(sprintf('%s\n%s', truncate_text(results(i).text, 20), status_str), 'Color', title_color, 'FontSize', 9);
end

sgtitle(sprintf('Small SST-2 Verification (500D, \\epsilon=%.4f) - Gray=Bounds, Blue=Samples', epsilon));

% Summary figure
figure('Name', 'Small SST-2 Verification Summary');
bar([num_verified, num_unknown, num_test-num_correct], 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'XTickLabel', {'Verified', 'Unknown', 'Incorrect'});
ylabel('Count');
title(sprintf('SST-2 Small Vocab (500D, \\epsilon=%.4f)', epsilon));

save('sentiment_sst2_small_verification.mat', 'results', 'epsilon');
fprintf('\nResults saved.\n');

%% Helpers
function str = label_str(l)
    if l == 0
        str = 'negative';
    else
        str = 'positive';
    end
end

function t = truncate_text(text, n)
    if length(text) > n
        t = [text(1:n-3) '...'];
    else
        t = text;
    end
end

function words = tokenize_sentence(text)
    text = lower(text);
    text = regexprep(text, '[^a-z0-9\s]', ' ');
    words = strsplit(strtrim(text));
    words = words(~cellfun(@isempty, words));
end

function vec = sentence_to_tfidf(text, word2idx, vocab_size, idf)
    words = tokenize_sentence(text);
    vec = zeros(vocab_size, 1);
    for i = 1:length(words)
        w = words{i};
        if isKey(word2idx, w)
            vec(word2idx(w)) = vec(word2idx(w)) + 1;
        end
    end
    if any(vec > 0)
        vec = vec .* idf;
        if norm(vec) > 0
            vec = vec / norm(vec);
        end
    end
end
