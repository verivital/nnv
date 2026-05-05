%% Robustness Verification of Vision Transformer on MNIST
%
% This script demonstrates neural network verification of a
% ViT-inspired model trained on MNIST digit classification.
%
% Verification Property:
%   For an input image x classified as digit d, verify that
%   all perturbed images x' within L-infinity ball ||x - x'||_inf <= eps
%   are also classified as digit d.
%
% Author: NNV Team
% Date: November 2025

%% Configuration
epsilon = 0.5;  % L-infinity perturbation (pixel values 0-255) - reduced for tighter bounds
num_images_to_verify = 5;  % Number of images to verify

fprintf('=== MNIST ViT Robustness Verification ===\n');
fprintf('Perturbation epsilon: %d (L-infinity)\n', epsilon);
fprintf('Images to verify: %d\n', num_images_to_verify);

%% 1) Load Trained Model
fprintf('\n[1/4] Loading trained ViT model...\n');

% Check if model exists
model_path = fullfile(fileparts(mfilename('fullpath')), 'mnist_vit_model.mat');
if ~isfile(model_path)
    error('Model not found. Please run train_mnist_vit.m first.');
end

load(model_path, 'net', 'accuracy', 'config');
fprintf('Model loaded. Training accuracy: %.2f%%\n', accuracy * 100);

% Convert to NNV
fprintf('Converting to NNV format...\n');
try
    nnv_net = matlab2nnv(net);
    fprintf('NNV model created with %d layers\n', length(nnv_net.Layers));
catch ME
    fprintf('Note: %s\n', ME.message);
    fprintf('Using MATLAB network for evaluation.\n');
    nnv_net = [];
end

%% 2) Load Test Data
fprintf('\n[2/4] Loading test data...\n');

% Use MATLAB's built-in digit dataset
digitDatasetPath = fullfile(matlabroot, 'toolbox', 'nnet', 'nndemos', ...
    'nndatasets', 'DigitDataset');

imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders', true, ...
    'LabelSource', 'foldernames');

% Get some test images
rng(42);  % For reproducibility
test_indices = randperm(numel(imds.Files), num_images_to_verify);

fprintf('Selected %d test images\n', num_images_to_verify);

%% 3) Verification Loop
fprintf('\n[3/4] Running verification...\n');

results = struct('index', {}, 'true_label', {}, 'predicted_label', {}, ...
    'verified', {}, 'time', {}, 'output_bounds', {});

for i = 1:num_images_to_verify
    idx = test_indices(i);

    % Load image
    [img, fileInfo] = readimage(imds, idx);
    img = single(img);
    if size(img, 3) == 3
        img = rgb2gray(img);
    end

    true_label = string(fileInfo.Label);

    fprintf('\n--- Image %d (index %d) ---\n', i, idx);
    fprintf('True label: %s\n', true_label);

    % Use raw image (0-255 range) - ImageInputLayer handles normalization
    img_raw = single(img);

    % Evaluate original image
    if ~isempty(nnv_net)
        % Pass raw image - NNV ImageInputLayer will apply rescale-zero-one normalization
        y_output = nnv_net.evaluate(img_raw);
    else
        y_output = predict(net, img);
        y_output = y_output(:);
    end

    [~, pred_idx] = max(y_output);
    predicted_label = num2str(pred_idx - 1);  % 0-indexed labels
    fprintf('Predicted label: %s\n', predicted_label);

    % Create input set with L-infinity perturbation (in 0-255 space)
    % ImageInputLayer will normalize during reachability
    lb_raw = max(img_raw - epsilon, 0);
    ub_raw = min(img_raw + epsilon, 255);

    % Create ImageStar (raw pixel values, normalization handled by ImageInputLayer)
    IS = ImageStar(lb_raw, ub_raw);

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

            % Get output bounds - R is the output set directly (not cell array)
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
                lb_out = inf(10, 1);
                ub_out = -inf(10, 1);
                for j = 1:length(R_out)
                    [lb_temp, ub_temp] = R_out(j).getRanges();
                    lb_out = min(lb_out, lb_temp(:));
                    ub_out = max(ub_out, ub_temp(:));
                end
            end

            % Check verification
            % Robust if: predicted class has highest lower bound > all other upper bounds
            target_idx = pred_idx;
            other_max_ub = max(ub_out([1:target_idx-1, target_idx+1:end]));

            if lb_out(target_idx) > other_max_ub
                verified = 1;  % Verified robust
                fprintf('VERIFIED: Network is robust!\n');
            else
                verified = 0;  % Unknown (could be robust or not)
                fprintf('UNKNOWN: Cannot verify robustness.\n');
                fprintf('  Target class %d: [%.4f, %.4f]\n', target_idx-1, lb_out(target_idx), ub_out(target_idx));
                fprintf('  Max other upper bound: %.4f\n', other_max_ub);
            end
        else
            % Fallback: Sample-based verification
            fprintf('Using sample-based verification (no NNV model)...\n');

            num_samples = 100;
            all_same = true;

            for s = 1:num_samples
                % Sample from perturbation set (already in [0,1] space)
                img_perturbed = lb_norm + (ub_norm - lb_norm) .* rand(size(img_normalized));

                y_sample = predict(net, img_perturbed);
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

    % Store results
    results(i).index = idx;
    results(i).true_label = true_label;
    results(i).predicted_label = predicted_label;
    results(i).verified = verified;
    results(i).time = reach_time;
    results(i).output_bounds = struct('lb', lb_out, 'ub', ub_out);
end

%% 4) Summary
fprintf('\n[4/4] Verification Summary\n');
fprintf('========================\n');

num_verified = sum([results.verified] == 1);
num_unknown = sum([results.verified] == 0 | [results.verified] == 2);
num_failed = sum([results.verified] == -1);

fprintf('Total images: %d\n', num_images_to_verify);
fprintf('Verified robust: %d\n', num_verified);
fprintf('Unknown/Not robust: %d\n', num_unknown);
fprintf('Failed: %d\n', num_failed);
fprintf('Average time: %.2f seconds\n', mean([results.time]));

%% 5) Visualize Results
fprintf('\n[5] Visualizing results...\n');

figure('Name', 'MNIST ViT Verification Results');

for i = 1:min(num_images_to_verify, 5)
    idx = test_indices(i);
    [img, ~] = readimage(imds, idx);

    subplot(2, num_images_to_verify, i);
    imshow(img, [0 255]);

    if results(i).verified == 1
        title_str = sprintf('%s (ROBUST)', results(i).predicted_label);
        title(title_str, 'Color', 'g');
    elseif results(i).verified == 0
        title_str = sprintf('%s (UNKNOWN)', results(i).predicted_label);
        title(title_str, 'Color', 'r');
    else
        title_str = sprintf('%s', results(i).predicted_label);
        title(title_str);
    end

    % Plot output bounds if available
    subplot(2, num_images_to_verify, num_images_to_verify + i);
    if ~isempty(results(i).output_bounds.lb)
        lb = results(i).output_bounds.lb;
        ub = results(i).output_bounds.ub;
        mid = (lb + ub) / 2;
        err = (ub - lb) / 2;

        errorbar(0:9, mid, err, '.');
        xlim([-0.5, 9.5]);
        xlabel('Digit');
        ylabel('Output');
        title('Output Bounds');
    end
end

% Save results
save('verification_results.mat', 'results', 'epsilon', 'num_images_to_verify');
fprintf('\nResults saved to: verification_results.mat\n');

fprintf('\n=== Verification Complete ===\n');
