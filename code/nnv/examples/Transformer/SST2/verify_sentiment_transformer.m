%% Robustness Verification of Sentiment Transformer
%
% This script demonstrates neural network verification of a
% Transformer-style sentiment classifier.
%
% Verification Property:
%   For an input text with sentiment s, verify that small perturbations
%   to the bag-of-words representation don't change the classification.
%
% Perturbation Model:
%   L-infinity perturbation in the normalized BoW feature space.
%   This models uncertainty in word importance weights.
%
% Author: NNV Team
% Date: November 2025

%% Configuration
epsilon = 0.01;  % L-infinity perturbation in feature space (small for tighter bounds)
num_samples_to_verify = 5;  % Number of samples to verify

fprintf('=== Sentiment Transformer Robustness Verification ===\n');
fprintf('Perturbation epsilon: %.2f (L-infinity in feature space)\n', epsilon);
fprintf('Samples to verify: %d\n', num_samples_to_verify);

%% 1) Load Trained Model
fprintf('\n[1/5] Loading trained sentiment model...\n');

% Check if model exists
model_path = fullfile(fileparts(mfilename('fullpath')), 'sentiment_transformer_model.mat');
if ~isfile(model_path)
    error('Model not found. Please run train_sentiment_transformer.m first.');
end

load(model_path, 'net', 'accuracy', 'config', 'word2idx', 'all_vocab');
fprintf('Model loaded. Training accuracy: %.2f%%\n', accuracy * 100);
fprintf('Vocabulary size: %d words\n', length(all_vocab));

% Convert to NNV
fprintf('Converting to NNV format...\n');
try
    nnv_net = matlab2nnv(net);
    fprintf('NNV model created with %d layers\n', length(nnv_net.Layers));
catch ME
    fprintf('NNV conversion failed: %s\n', ME.message);
    nnv_net = [];
end

%% 2) Create Test Samples
fprintf('\n[2/5] Creating test samples...\n');

% Define test sentences with known sentiment
test_sentences = {
    'great movie', 'positive'
    'excellent film', 'positive'
    'wonderful story', 'positive'
    'amazing performance', 'positive'
    'fantastic experience', 'positive'
    'terrible movie', 'negative'
    'awful film', 'negative'
    'boring story', 'negative'
    'bad performance', 'negative'
    'disappointing experience', 'negative'
};

% Convert to BoW representation
vocab_size = length(all_vocab);
num_test = min(num_samples_to_verify, size(test_sentences, 1));

X_test = zeros(vocab_size, 1, 1, num_test, 'single');
test_labels = cell(num_test, 1);
test_texts = cell(num_test, 1);

for i = 1:num_test
    bow = text_to_bow(test_sentences{i, 1}, word2idx, vocab_size);
    X_test(:, 1, 1, i) = single(bow);
    test_labels{i} = test_sentences{i, 2};
    test_texts{i} = test_sentences{i, 1};
end

fprintf('Created %d test samples\n', num_test);

%% 3) Verification Loop
fprintf('\n[3/5] Running verification...\n');

results = struct('index', {}, 'text', {}, 'true_label', {}, ...
    'predicted_label', {}, 'verified', {}, 'time', {}, 'output_bounds', {});

for i = 1:num_test
    fprintf('\n--- Sample %d ---\n', i);
    fprintf('Text: "%s"\n', test_texts{i});
    fprintf('True label: %s\n', test_labels{i});

    % Get input features
    x = X_test(:, 1, 1, i);
    x_img = reshape(x, [vocab_size, 1, 1]);  % ImageStar-compatible shape

    % Evaluate original input (use image shape to match reachability)
    if ~isempty(nnv_net)
        y_output = nnv_net.evaluate(x_img);
    else
        y_output = predict(net, x_img);
        y_output = y_output(:);
    end

    % Determine prediction (logits, so higher = more likely)
    [~, pred_idx] = max(y_output);
    if pred_idx == 1
        predicted_label = 'positive';
    else
        predicted_label = 'negative';
    end
    fprintf('Predicted label: %s\n', predicted_label);

    % Create input set with L-infinity perturbation
    % Only perturb non-zero features (words that appear in the text)
    lb = max(x - epsilon, 0);  % Features are non-negative
    ub = min(x + epsilon, 1);  % Normalized to [0, 1]

    % Reshape to image format [H x W x C] for ImageStar
    lb_img = reshape(lb, [vocab_size, 1, 1]);
    ub_img = reshape(ub, [vocab_size, 1, 1]);

    % Create ImageStar from bounds (network uses imageInputLayer)
    IS = ImageStar(lb_img, ub_img);

    % Run verification
    fprintf('Running reachability analysis...\n');
    t = tic;

    try
        if ~isempty(nnv_net)
            % Set up reach options
            reachOptions = struct();
            reachOptions.reachMethod = 'approx-star';

            % Compute reachable set
            R = nnv_net.reach(IS, reachOptions);

            % Get output bounds - check ImageStar first, then Star, then array
            R_out = R;
            if isa(R_out, 'ImageStar')
                [lb_out, ub_out] = R_out.getRanges();
                lb_out = squeeze(lb_out);
                ub_out = squeeze(ub_out);
            elseif isa(R_out, 'Star')
                n = R_out.dim;
                lb_out = zeros(n, 1);
                ub_out = zeros(n, 1);
                for j = 1:n
                    lb_out(j) = R_out.getMin(j, 'linprog');
                    ub_out(j) = R_out.getMax(j, 'linprog');
                end
            else
                % Array of stars
                lb_out = inf(2, 1);
                ub_out = -inf(2, 1);
                for j = 1:length(R_out)
                    [lb_temp, ub_temp] = R_out(j).getRanges();
                    lb_out = min(lb_out, lb_temp(:));
                    ub_out = max(ub_out, ub_temp(:));
                end
            end

            % Check verification
            % Robust if: predicted class has highest lower bound > all other upper bounds
            target_idx = pred_idx;
            if target_idx == 1
                other_max_ub = ub_out(2);
            else
                other_max_ub = ub_out(1);
            end

            if lb_out(target_idx) > other_max_ub
                verified = 1;  % Verified robust
                fprintf('VERIFIED: Network is robust!\n');
            else
                verified = 0;  % Unknown
                fprintf('UNKNOWN: Cannot verify robustness.\n');
                fprintf('  Target class (%s): [%.4f, %.4f]\n', predicted_label, lb_out(target_idx), ub_out(target_idx));
                fprintf('  Other class upper bound: %.4f\n', other_max_ub);
            end
        else
            % Fallback: Sample-based verification
            fprintf('Using sample-based verification (no NNV model)...\n');

            num_samples = 100;
            all_same = true;

            for s = 1:num_samples
                % Sample from perturbation set
                x_perturbed = lb + (ub - lb) .* rand(size(x));
                x_perturbed = reshape(x_perturbed, [vocab_size, 1, 1]);

                y_sample = predict(net, x_perturbed);
                [~, sample_pred] = max(y_sample);

                if sample_pred ~= pred_idx
                    all_same = false;
                    fprintf('Counterexample found at sample %d\n', s);
                    break;
                end
            end

            if all_same
                verified = 2;  % Unknown (no counterexample found)
                fprintf('No counterexample found in %d samples\n', num_samples);
            else
                verified = 0;
                fprintf('NOT ROBUST: Found counterexample\n');
            end

            lb_out = y_output;
            ub_out = y_output;
        end

        reach_time = toc(t);
        fprintf('Verification time: %.2f seconds\n', reach_time);

    catch ME
        fprintf('Verification failed: %s\n', ME.message);
        verified = -1;
        reach_time = toc(t);
        lb_out = [];
        ub_out = [];
    end

    % Sample random inputs from perturbation set for soundness visualization
    num_random_samples = 50;
    sampled_outputs = zeros(2, num_random_samples);
    for s = 1:num_random_samples
        x_sample = lb + (ub - lb) .* rand(size(x));
        x_sample_img = reshape(x_sample, [vocab_size, 1, 1]);  % Match ImageStar shape
        if ~isempty(nnv_net)
            y_sample = nnv_net.evaluate(x_sample_img);
        else
            y_sample = predict(net, x_sample_img);
            y_sample = y_sample(:);
        end
        sampled_outputs(:, s) = y_sample;
    end

    % Store results
    results(i).index = i;
    results(i).text = test_texts{i};
    results(i).true_label = test_labels{i};
    results(i).predicted_label = predicted_label;
    results(i).verified = verified;
    results(i).time = reach_time;
    results(i).output_bounds = struct('lb', lb_out, 'ub', ub_out);
    results(i).sampled_outputs = sampled_outputs;
end

%% 4) Summary
fprintf('\n[4/5] Verification Summary\n');
fprintf('========================\n');

num_verified = sum([results.verified] == 1);
num_unknown = sum([results.verified] == 0 | [results.verified] == 2);
num_failed = sum([results.verified] == -1);
num_correct_pred = sum(strcmp({results.true_label}, {results.predicted_label}));

fprintf('Total samples: %d\n', num_test);
fprintf('Correct predictions: %d\n', num_correct_pred);
fprintf('Verified robust: %d\n', num_verified);
fprintf('Unknown/Not robust: %d\n', num_unknown);
fprintf('Failed: %d\n', num_failed);
fprintf('Average time: %.2f seconds\n', mean([results.time]));

%% 5) Visualization - Soundness Check
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

figure('Name', 'Sentiment Transformer Verification - Soundness Check', 'Position', [100, 100, 1200, 400]);

for i = 1:num_test
    subplot(1, num_test, i);
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
    set(gca, 'XTick', [1, 2], 'XTickLabel', {'Negative', 'Positive'});
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
    title(sprintf('"%s"\n%s', results(i).text, status_str), 'Color', title_color, 'FontSize', 9);
end

sgtitle(sprintf('Toy SST Verification (\\epsilon=%.2f) - Gray=Bounds, Blue=Samples', epsilon));

% Summary figure
figure('Name', 'Sentiment Transformer Verification Summary');
status_counts = [num_verified, num_unknown, num_failed];
bar(status_counts, 'FaceColor', [0.2, 0.6, 0.8]);
set(gca, 'XTickLabel', {'Verified', 'Unknown', 'Failed'});
ylabel('Count');
title(sprintf('Verification Summary (\\epsilon=%.2f)', epsilon));
ylim([0, max(status_counts) + 1]);

% Save results
save('sentiment_verification_results.mat', 'results', 'epsilon', 'num_test');
fprintf('\nResults saved to: sentiment_verification_results.mat\n');

fprintf('\n=== Verification Complete ===\n');

%% Helper Functions
function bow = text_to_bow(text, word2idx, vocab_size)
    % Convert text to bag-of-words vector
    bow = zeros(vocab_size, 1);
    words = strsplit(lower(text));

    for i = 1:length(words)
        if isKey(word2idx, words{i})
            idx = word2idx(words{i});
            bow(idx) = bow(idx) + 1;
        end
    end

    % Normalize by L2 norm
    if norm(bow) > 0
        bow = bow / norm(bow);
    end
end
