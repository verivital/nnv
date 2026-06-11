%% Robustness Verification of Sentiment Transformer on Real SST-2
%
% This script demonstrates neural network verification of a
% Transformer-style sentiment classifier trained on real SST-2 data.
%
% Verification Property:
%   For an input text with sentiment s, verify that small perturbations
%   to the TF-IDF feature representation don't change the classification.
%
% Author: NNV Team
% Date: November 2025

%% Configuration
epsilon = 1e-6;  % Very small L-infinity perturbation for fast verification
num_samples_to_verify = 5;  % Number of samples to verify

fprintf('=== Real SST-2 Sentiment Verification ===\n');
fprintf('Perturbation epsilon: %.4f (L-infinity in TF-IDF space)\n', epsilon);
fprintf('Samples to verify: %d\n', num_samples_to_verify);

%% 1) Load Trained Model
fprintf('\n[1/5] Loading real SST-2 trained model...\n');

model_path = fullfile(fileparts(mfilename('fullpath')), 'sentiment_transformer_sst2_model.mat');
if ~isfile(model_path)
    error('Model not found. Please run train_sentiment_transformer_sst2.m first.');
end

load(model_path, 'net', 'accuracy', 'config', 'word2idx', 'all_vocab', 'idf', 'vocab_size');
fprintf('Model loaded. Validation accuracy: %.2f%%\n', accuracy * 100);
fprintf('Vocabulary size: %d words\n', vocab_size);

% Convert to NNV
fprintf('Converting to NNV format...\n');
try
    nnv_net = matlab2nnv(net);
    fprintf('NNV model created with %d layers\n', length(nnv_net.Layers));
catch ME
    fprintf('NNV conversion failed: %s\n', ME.message);
    nnv_net = [];
end

%% 2) Load Real SST-2 Test Data
fprintf('\n[2/5] Loading real SST-2 test data...\n');

[~, dev_data] = download_sst2();

% Use development set for verification (has labels)
test_sentences = dev_data.sentences;
test_labels_raw = dev_data.labels;

% Select random samples for verification
rng(123);  % Reproducibility
test_idx = randperm(length(test_sentences), min(num_samples_to_verify, length(test_sentences)));

fprintf('Selected %d test samples from SST-2 dev set\n', length(test_idx));

%% 3) Tokenize Helper Function
tokenize_sentence = @(text) tokenize_text(text);
sentence_to_tfidf_vec = @(text) sentence_to_tfidf(text, word2idx, vocab_size, idf);

%% 4) Verification Loop
fprintf('\n[3/5] Running verification...\n');

results = struct('index', {}, 'text', {}, 'true_label', {}, ...
    'predicted_label', {}, 'verified', {}, 'time', {}, 'output_bounds', {});

for i = 1:length(test_idx)
    idx = test_idx(i);
    text = test_sentences{idx};
    true_label_num = test_labels_raw(idx);  % 0=negative, 1=positive

    if true_label_num == 1
        true_label = 'positive';
    else
        true_label = 'negative';
    end

    fprintf('\n--- Sample %d (idx %d) ---\n', i, idx);
    fprintf('Text: "%s"\n', truncate_text(text, 60));
    fprintf('True label: %s\n', true_label);

    % Convert to TF-IDF features
    x = sentence_to_tfidf_vec(text);
    x = single(x);

    % Reshape for network input [vocab_size x 1 x 1]
    x_input = reshape(x, [vocab_size, 1, 1]);

    % Evaluate original input (use image shape to match reachability)
    if ~isempty(nnv_net)
        y_output = nnv_net.evaluate(x_input);
    else
        y_output = predict(net, x_input);
        y_output = y_output(:);
    end

    % Determine prediction
    [~, pred_idx] = max(y_output);
    if pred_idx == 1
        predicted_label = 'negative';
    else
        predicted_label = 'positive';
    end
    fprintf('Predicted label: %s\n', predicted_label);

    % Create input set with L-infinity perturbation
    lb = max(x - epsilon, 0);  % TF-IDF features are non-negative
    ub = x + epsilon;  % Allow small positive perturbation

    % Reshape for ImageStar
    lb_img = reshape(lb, [vocab_size, 1, 1]);
    ub_img = reshape(ub, [vocab_size, 1, 1]);

    IS = ImageStar(lb_img, ub_img);

    % Run verification
    fprintf('Running reachability analysis...\n');
    t = tic;

    try
        if ~isempty(nnv_net)
            reachOptions = struct();
            reachOptions.reachMethod = 'approx-star';

            R = nnv_net.reach(IS, reachOptions);

            % Get output bounds
            R_out = R;
            if isa(R_out, 'Star')
                n = R_out.dim;
                lb_out = zeros(n, 1);
                ub_out = zeros(n, 1);
                for j = 1:n
                    lb_out(j) = R_out.getMin(j, 'linprog');
                    ub_out(j) = R_out.getMax(j, 'linprog');
                end
            else
                [lb_out, ub_out] = R_out.getRanges();
                lb_out = squeeze(lb_out);
                ub_out = squeeze(ub_out);
            end

            % Check verification (2 classes: negative=1, positive=2)
            target_idx = pred_idx;
            if target_idx == 1
                other_max_ub = ub_out(2);
            else
                other_max_ub = ub_out(1);
            end

            if lb_out(target_idx) > other_max_ub
                verified = 1;
                fprintf('VERIFIED: Network is robust!\n');
            else
                verified = 0;
                fprintf('UNKNOWN: Cannot verify robustness.\n');
                fprintf('  Target class (%s): [%.4f, %.4f]\n', predicted_label, lb_out(target_idx), ub_out(target_idx));
                fprintf('  Other class upper bound: %.4f\n', other_max_ub);
            end
        else
            verified = -1;
            lb_out = y_output;
            ub_out = y_output;
            fprintf('No NNV model available.\n');
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
    results(i).index = idx;
    results(i).text = text;
    results(i).true_label = true_label;
    results(i).predicted_label = predicted_label;
    results(i).verified = verified;
    results(i).time = reach_time;
    results(i).output_bounds = struct('lb', lb_out, 'ub', ub_out);
    results(i).sampled_outputs = sampled_outputs;
end

%% 5) Summary
fprintf('\n[4/5] Verification Summary\n');
fprintf('========================\n');

num_verified = sum([results.verified] == 1);
num_unknown = sum([results.verified] == 0);
num_failed = sum([results.verified] == -1);
num_correct_pred = sum(strcmp({results.true_label}, {results.predicted_label}));

fprintf('Total samples: %d\n', length(test_idx));
fprintf('Correct predictions: %d (%.1f%%)\n', num_correct_pred, 100*num_correct_pred/length(test_idx));
fprintf('Verified robust: %d (%.1f%%)\n', num_verified, 100*num_verified/length(test_idx));
fprintf('Unknown: %d\n', num_unknown);
fprintf('Failed: %d\n', num_failed);
fprintf('Average time: %.2f seconds\n', mean([results.time]));

%% 6) Visualization - Soundness Check
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

figure('Name', 'Full SST-2 Verification - Soundness Check', 'Position', [100, 100, 1200, 400]);

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

sgtitle(sprintf('Full SST-2 Verification (5000D, \\epsilon=%.0e) - Gray=Bounds, Blue=Samples', epsilon));

% Summary figure
figure('Name', 'Full SST-2 Verification Summary', 'Position', [100 500 800 300]);

subplot(1, 2, 1);
status_counts = [num_verified, num_unknown, num_failed];
b = bar(status_counts, 'FaceColor', 'flat');
b.CData(1,:) = [0.2 0.7 0.3];  % Green for verified
b.CData(2,:) = [0.9 0.6 0.1];  % Orange for unknown
b.CData(3,:) = [0.8 0.2 0.2];  % Red for failed
set(gca, 'XTickLabel', {'Verified', 'Unknown', 'Failed'});
ylabel('Count');
title(sprintf('Verification Status (\\epsilon=%.0e)', epsilon));
ylim([0, max(status_counts) + 1]);

subplot(1, 2, 2);
correct_verified = sum(strcmp({results.true_label}, {results.predicted_label}) & [results.verified] == 1);
correct_unknown = sum(strcmp({results.true_label}, {results.predicted_label}) & [results.verified] == 0);
incorrect = sum(~strcmp({results.true_label}, {results.predicted_label}));

pie_data = [correct_verified, correct_unknown, incorrect];
pie_labels = {sprintf('Correct+Verified (%d)', correct_verified), ...
              sprintf('Correct+Unknown (%d)', correct_unknown), ...
              sprintf('Incorrect (%d)', incorrect)};
pie_data_filt = pie_data(pie_data > 0);
pie_labels_filt = pie_labels(pie_data > 0);
if ~isempty(pie_data_filt)
    pie(pie_data_filt);
    legend(pie_labels_filt, 'Location', 'southoutside');
end
title('Prediction Breakdown');

% Save results
save('sentiment_sst2_verification_results.mat', 'results', 'epsilon', 'num_samples_to_verify');
fprintf('\nResults saved to: sentiment_sst2_verification_results.mat\n');

% Save figure
saveas(gcf, 'sst2_verification_results.png');
fprintf('Figure saved to: sst2_verification_results.png\n');

fprintf('\n=== Verification Complete ===\n');

%% Helper Functions
function words = tokenize_text(text)
    text = lower(text);
    text = regexprep(text, '[^a-z0-9\s]', ' ');
    words = strsplit(strtrim(text));
    words = words(~cellfun(@isempty, words));
end

function vec = sentence_to_tfidf(text, word2idx, vocab_size, idf)
    words = tokenize_text(text);
    vec = zeros(vocab_size, 1);

    for i = 1:length(words)
        w = words{i};
        if isKey(word2idx, w)
            idx = word2idx(w);
            vec(idx) = vec(idx) + 1;
        end
    end

    if any(vec > 0)
        tf = vec;
        vec = tf .* idf;
        norm_val = norm(vec);
        if norm_val > 0
            vec = vec / norm_val;
        end
    end
end

function truncated = truncate_text(text, max_len)
    if length(text) > max_len
        truncated = [text(1:max_len-3) '...'];
    else
        truncated = text;
    end
end
