%% Verify MNIST ViT with Real Self-Attention using NNV
% This script performs robustness verification on the ViT model
% that uses real MultiHeadAttentionLayer
%
% Author: NNV Team
% Date: December 2025

%% Load the trained model
fprintf('=== Verifying MNIST ViT with Real Attention ===\n\n');

model_file = 'mnist_vit_attention_model.mat';
if ~isfile(model_file)
    error('Model not found. Run train_mnist_vit_attention.m first.');
end

load(model_file, 'net', 'accuracy', 'config');
fprintf('Loaded model with %.2f%% validation accuracy\n', accuracy * 100);

%% Convert to NNV
fprintf('\n[1/4] Converting to NNV...\n');

nnv_net = matlab2nnv(net);
fprintf('Converted network has %d layers\n', length(nnv_net.Layers));

% Identify attention layer
attention_idx = -1;
for i = 1:length(nnv_net.Layers)
    if isa(nnv_net.Layers{i}, 'MultiHeadAttentionLayer')
        attention_idx = i;
        fprintf('MultiHeadAttentionLayer found at index %d\n', i);
        break;
    end
end

%% Load test images
fprintf('\n[2/4] Loading test images...\n');

digitDatasetPath = fullfile(matlabroot, 'toolbox', 'nnet', 'nndemos', ...
    'nndatasets', 'DigitDataset');
imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders', true, ...
    'LabelSource', 'foldernames');

% Use augmentedImageDatastore for proper preprocessing (same as training)
augimds = augmentedImageDatastore([28 28 1], imds);

% Get a few test images (one per class)
num_test = 5;
rng(42); % for reproducibility
test_indices = randperm(numel(imds.Files), num_test);

fprintf('Selected %d test images for verification\n', num_test);

%% Define perturbation
fprintf('\n[3/4] Setting up verification...\n');

% Small L-infinity perturbation (pixel values are normalized to [0,1])
epsilon = 0.02;  % 2% perturbation
fprintf('Perturbation epsilon: %.3f (L-inf norm)\n', epsilon);

%% Run verification
fprintf('\n[4/4] Running robustness verification...\n\n');

results = struct();
results.verified = 0;
results.not_verified = 0;
results.errors = 0;
results.times = [];

for i = 1:num_test
    idx = test_indices(i);

    % Load and preprocess image using augmentedImageDatastore (same as training)
    data = readByIndex(augimds, idx);
    img = data.input{1};  % Get the preprocessed image

    % Ensure correct size [28 28 1]
    if size(img, 3) ~= 1
        img = rgb2gray(img);
    end

    % Get true label
    true_label = imds.Labels(idx);
    true_class = double(true_label);

    % Predict with original network
    pred = classify(net, img);

    fprintf('Image %d: True label = %s, Predicted = %s\n', i, string(true_label), string(pred));

    if pred ~= true_label
        fprintf('  -> Skipping (incorrect prediction)\n');
        continue;
    end

    % Create input set with perturbation
    lb = max(0, img - epsilon);
    ub = min(1, img + epsilon);

    % Create ImageStar input set
    input_set = ImageStar(lb, ub);

    % Run reachability analysis
    fprintf('  Running reachability analysis...\n');

    try
        tic;

        % Use approximate reachability for speed
        reachOptions = struct('reachMethod', 'approx-star');
        output_set = nnv_net.reach(input_set, reachOptions);

        reach_time = toc;
        results.times = [results.times, reach_time];

        fprintf('  Reachability completed in %.2f seconds\n', reach_time);

        % Handle cell array output
        if iscell(output_set)
            output_set = output_set{1};
        end

        % Check robustness: verify that true class has highest output
        % Get output bounds based on output type
        if isa(output_set, 'ImageStar')
            if ~isempty(output_set.im_lb) && ~isempty(output_set.im_ub)
                lb_out = output_set.im_lb(:);
                ub_out = output_set.im_ub(:);
            else
                % Use V representation - estimate ranges
                [lb_img, ub_img] = output_set.estimateRanges();
                lb_out = lb_img(:);
                ub_out = ub_img(:);
            end
        elseif isa(output_set, 'Star')
            [lb_out, ub_out] = output_set.getRanges();
        else
            error('Unknown output type: %s', class(output_set));
        end

        % Check if true class lower bound > all other class upper bounds
        true_class_lb = lb_out(true_class);
        other_classes_ub = ub_out;
        other_classes_ub(true_class) = -inf;
        max_other_ub = max(other_classes_ub);

        if true_class_lb > max_other_ub
            fprintf('  -> VERIFIED ROBUST\n');
            results.verified = results.verified + 1;
        else
            fprintf('  -> NOT VERIFIED (bounds overlap)\n');
            fprintf('     True class [%d] lb: %.4f, Max other ub: %.4f\n', ...
                true_class, true_class_lb, max_other_ub);
            results.not_verified = results.not_verified + 1;
        end

    catch ME
        fprintf('  -> ERROR: %s\n', ME.message);
        results.errors = results.errors + 1;
    end

    fprintf('\n');
end

%% Summary
fprintf('============================================\n');
fprintf('           Verification Summary\n');
fprintf('============================================\n');
fprintf('Model: MNIST ViT with Real Attention\n');
fprintf('Perturbation: epsilon = %.3f\n', epsilon);
fprintf('Method: approx-star\n');
fprintf('--------------------------------------------\n');
fprintf('Total images tested: %d\n', num_test);
fprintf('Verified robust: %d\n', results.verified);
fprintf('Not verified: %d\n', results.not_verified);
fprintf('Errors: %d\n', results.errors);
if ~isempty(results.times)
    fprintf('Average reach time: %.2f seconds\n', mean(results.times));
end
fprintf('============================================\n');

%% Save results
save('verification_results_vit_attention.mat', 'results', 'epsilon', 'num_test');
fprintf('\nResults saved to verification_results_vit_attention.mat\n');
