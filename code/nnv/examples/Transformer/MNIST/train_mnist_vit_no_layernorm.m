%% Train MNIST ViT without LayerNormalization
% Vision Transformer for MNIST with real self-attention but NO LayerNorm
% This is to isolate the impact of LayerNorm on reachability bound looseness
%
% Architecture:
%   - Patch embedding via Conv2D
%   - Self-attention (without pre/post LayerNorm)
%   - FFN with ReLU (without LayerNorm)
%   - Linear classifier
%
% Author: NNV Team
% Date: December 2025

fprintf('=== Training MNIST ViT without LayerNormalization ===\n\n');

%% Configuration
config = struct();
config.patch_size = 4;          % 4x4 patches -> 7x7 = 49 patches
config.embed_dim = 64;          % Smaller embedding (no LayerNorm needs smaller network)
config.num_heads = 2;           % Attention heads
config.ffn_dim = 128;           % FFN hidden dimension
config.num_classes = 10;        % MNIST digits
config.epochs = 15;             % Training epochs
config.batch_size = 128;        % Batch size
config.learning_rate = 0.0005;  % Lower learning rate (no LayerNorm needs careful tuning)

fprintf('Configuration:\n');
fprintf('  Patch size: %dx%d\n', config.patch_size, config.patch_size);
fprintf('  Embedding dim: %d\n', config.embed_dim);
fprintf('  Attention heads: %d\n', config.num_heads);
fprintf('  FFN dim: %d\n', config.ffn_dim);
fprintf('  Learning rate: %.4f\n', config.learning_rate);
fprintf('  Epochs: %d\n', config.epochs);

%% Load MNIST Data
fprintf('\n[Data] Loading MNIST...\n');

digitDatasetPath = fullfile(matlabroot, 'toolbox', 'nnet', 'nndemos', ...
    'nndatasets', 'DigitDataset');

imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders', true, ...
    'LabelSource', 'foldernames');

% Split into training and validation
[imdsTrain, imdsVal] = splitEachLabel(imds, 0.8, 'randomized');

% Create augmented datastores
augmenter = imageDataAugmenter('RandXReflection', false);
augimdsTrain = augmentedImageDatastore([28 28 1], imdsTrain, ...
    'DataAugmentation', augmenter);
augimdsVal = augmentedImageDatastore([28 28 1], imdsVal);

fprintf('  Training: %d images\n', numel(imdsTrain.Labels));
fprintf('  Validation: %d images\n', numel(imdsVal.Labels));

%% Define Network Architecture (NO LayerNorm)
fprintf('\n[Model] Defining ViT without LayerNormalization...\n');

% Input: 28x28x1 MNIST image
img_size = [28, 28, 1];
num_patches = (img_size(1) / config.patch_size) * (img_size(2) / config.patch_size);
flattened_dim = num_patches * config.embed_dim;

fprintf('  Image size: %dx%dx%d\n', img_size(1), img_size(2), img_size(3));
fprintf('  Number of patches: %d\n', num_patches);
fprintf('  Flattened dimension: %d\n', flattened_dim);

% Build the network
layers = [
    % Input
    imageInputLayer(img_size, 'Name', 'input', 'Normalization', 'none')

    % Patch embedding: Conv2D with stride = patch_size
    % Output: [7, 7, embed_dim] for 28x28 image with 4x4 patches
    convolution2dLayer(config.patch_size, config.embed_dim, ...
        'Stride', config.patch_size, 'Padding', 'same', 'Name', 'patch_embed')

    % Flatten to sequence: [49, embed_dim] -> [49*embed_dim, 1]
    flattenLayer('Name', 'flatten_patches')

    % NO LayerNorm here (removed!)

    % Self-attention (treats flattened input as sequence)
    selfAttentionLayer(config.num_heads, config.embed_dim, 'Name', 'self_attention')

    % NO LayerNorm after attention (removed!)

    % FFN: expand then contract with ReLU
    fullyConnectedLayer(config.ffn_dim, 'Name', 'ff1')
    reluLayer('Name', 'relu')
    fullyConnectedLayer(flattened_dim, 'Name', 'ff2')

    % NO LayerNorm before classifier (removed!)

    % Classification head
    fullyConnectedLayer(config.num_classes, 'Name', 'classifier')
    softmaxLayer('Name', 'softmax')
    classificationLayer('Name', 'output')
];

lgraph = layerGraph(layers);

fprintf('  Network created with %d layers (no LayerNorm)\n', numel(lgraph.Layers));

% Display layer summary
fprintf('\n  Layer structure:\n');
for i = 1:numel(lgraph.Layers)
    L = lgraph.Layers(i);
    fprintf('    %2d: %s (%s)\n', i, L.Name, class(L));
end

%% Training Options
fprintf('\n[Training] Setting up training...\n');

options = trainingOptions('adam', ...
    'InitialLearnRate', config.learning_rate, ...
    'MaxEpochs', config.epochs, ...
    'MiniBatchSize', config.batch_size, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', augimdsVal, ...
    'ValidationFrequency', 50, ...
    'Plots', 'none', ...
    'Verbose', true, ...
    'VerboseFrequency', 100, ...
    'ExecutionEnvironment', 'cpu');

%% Train the Network
fprintf('\n[Training] Training ViT without LayerNorm...\n');
tic;
[net, trainInfo] = trainNetwork(augimdsTrain, lgraph, options);
train_time = toc;

fprintf('\n[Training] Completed in %.1f seconds\n', train_time);

%% Evaluate
fprintf('\n[Evaluate] Testing network...\n');
YPred = classify(net, augimdsVal);
YTest = imdsVal.Labels;
accuracy = sum(YPred == YTest) / numel(YTest) * 100;
fprintf('  Test Accuracy: %.2f%%\n', accuracy);

%% Analyze Attention Weights
fprintf('\n[Analysis] Checking attention weights...\n');
attn_layer = [];
for i = 1:numel(net.Layers)
    if isa(net.Layers(i), 'nnet.cnn.layer.SelfAttentionLayer')
        attn_layer = net.Layers(i);
        break;
    end
end

if ~isempty(attn_layer)
    q_max = max(abs(attn_layer.QueryWeights(:)));
    k_max = max(abs(attn_layer.KeyWeights(:)));
    v_max = max(abs(attn_layer.ValueWeights(:)));
    o_max = max(abs(attn_layer.OutputWeights(:)));

    fprintf('  QueryWeights max|val|: %.4e\n', q_max);
    fprintf('  KeyWeights max|val|: %.4e\n', k_max);
    fprintf('  ValueWeights max|val|: %.4e\n', v_max);
    fprintf('  OutputWeights max|val|: %.4e\n', o_max);

    if q_max < 1e-6 && k_max < 1e-6
        fprintf('  WARNING: Q/K weights near-zero (attention bypassed)\n');
    else
        fprintf('  Attention weights are non-trivial\n');
    end
end

%% Save the network
fprintf('\n[Save] Saving network...\n');
example_dir = fileparts(mfilename('fullpath'));

% Save as MAT
save(fullfile(example_dir, 'mnist_vit_no_layernorm_model.mat'), 'net', 'accuracy', 'config', 'trainInfo');
fprintf('  Saved: mnist_vit_no_layernorm_model.mat\n');

% Export to ONNX
try
    exportONNXNetwork(net, fullfile(example_dir, 'mnist_vit_no_layernorm.onnx'));
    fprintf('  Saved: mnist_vit_no_layernorm.onnx\n');
catch ME
    fprintf('  ONNX export failed: %s\n', ME.message);
end

%% Quick Verification Test
fprintf('\n[Verify] Quick reachability test...\n');

nnv_net = matlab2nnv(net);
fprintf('  NNV network: %d layers\n', length(nnv_net.Layers));

% Load a test image using digitTest4DArrayData
[XTestData, ~] = digitTest4DArrayData;
img = double(XTestData(:,:,:,1));
epsilon = 1e-4;
lb = max(img - epsilon, 0);
ub = min(img + epsilon, 1);

input_set = ImageStar(lb, ub);
reachOptions = struct('reachMethod', 'approx-star');

try
    tic;
    R = nnv_net.reach(input_set, reachOptions);
    reach_time = toc;

    if iscell(R)
        R = R{1};
    end

    if isa(R, 'ImageStar')
        [out_lb, out_ub] = R.estimateRanges();
    else
        [out_lb, out_ub] = R.getRanges();
    end

    out_range = max(out_ub(:) - out_lb(:));
    fprintf('  Reachability time: %.2f seconds\n', reach_time);
    fprintf('  Output max range (eps=%.0e): %.2f\n', epsilon, out_range);

    % Compare to LayerNorm version
    fprintf('\n  Comparison (with LayerNorm version at same eps):\n');
    fprintf('    With LayerNorm:    max_range ~ 280 (from analysis)\n');
    fprintf('    Without LayerNorm: max_range = %.2f\n', out_range);

catch ME
    fprintf('  Reachability failed: %s\n', ME.message);
end

fprintf('\n=== Training Complete ===\n');
fprintf('Model saved to: mnist_vit_no_layernorm_model.mat\n');
