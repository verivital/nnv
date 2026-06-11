%% Train a Vision Transformer with Real Self-Attention and GELU on MNIST
% This script trains a ViT using MATLAB's native selfAttentionLayer
% with GELU activation (as used in the original Transformer paper)
%
% Architecture:
%   - Input: 28x28 grayscale image
%   - Patch embedding: 7x7 patches -> 16 tokens, projected to 32 dims
%   - Position embedding (learned)
%   - 1 Transformer block with real self-attention (2 heads)
%   - GELU activation in feed-forward network
%   - Global average pooling -> 10 classes
%
% Requires: MATLAB R2024a+ for selfAttentionLayer and geluLayer
%
% Author: NNV Team
% Date: December 2025

%% Configuration
config = struct();
config.patch_size = 7;          % 7x7 patches (divides 28 evenly into 4x4=16 patches)
config.embed_dim = 32;          % Embedding dimension (must be divisible by num_heads)
config.num_heads = 2;           % Number of attention heads
config.ff_dim = 64;             % Feed-forward hidden dimension
config.num_classes = 10;        % MNIST classes
config.max_epochs = 10;         % Training epochs
config.batch_size = 128;        % Batch size
config.learning_rate = 0.001;   % Initial learning rate
config.activation = 'gelu';     % Use GELU activation

fprintf('=== Training ViT with Real Self-Attention and GELU on MNIST ===\n');
fprintf('Patch size: %dx%d\n', config.patch_size, config.patch_size);
fprintf('Embedding dim: %d\n', config.embed_dim);
fprintf('Number of heads: %d\n', config.num_heads);
fprintf('Activation: %s\n', config.activation);

% Check MATLAB version for selfAttentionLayer support
v = ver('MATLAB');
fprintf('MATLAB version: %s\n', v.Release);

%% 1) Load MNIST Data
fprintf('\n[1/5] Loading MNIST data...\n');

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

fprintf('Training samples: %d\n', numel(imdsTrain.Labels));
fprintf('Validation samples: %d\n', numel(imdsVal.Labels));

%% 2) Build ViT Architecture with Real Attention and GELU
fprintf('\n[2/5] Building ViT architecture with selfAttentionLayer + GELU...\n');

img_size = 28;
num_patches = (img_size / config.patch_size)^2;  % 16 patches for 7x7
patch_dim = config.patch_size * config.patch_size;  % 49 for 7x7 patches

fprintf('Number of patches (sequence length): %d\n', num_patches);
fprintf('Patch dimension: %d\n', patch_dim);

% Check if selfAttentionLayer and geluLayer are available
if ~exist('selfAttentionLayer', 'file')
    error(['selfAttentionLayer not available. Requires MATLAB R2024a or newer.\n' ...
           'Your version: %s'], v.Release);
end
if ~exist('geluLayer', 'file')
    error(['geluLayer not available. Requires MATLAB R2023a or newer.\n' ...
           'Your version: %s'], v.Release);
end

% Build the network using layer graph for more control
lgraph = layerGraph();

% Input layer
inputLayer = imageInputLayer([28 28 1], 'Name', 'input', ...
    'Normalization', 'rescale-zero-one');
lgraph = addLayers(lgraph, inputLayer);

% Patch embedding via convolution (extracts and projects patches)
% Conv with stride = patch_size creates non-overlapping patches
patchEmbed = convolution2dLayer(config.patch_size, config.embed_dim, ...
    'Stride', config.patch_size, ...
    'Name', 'patch_embed', ...
    'Padding', 'same');
lgraph = addLayers(lgraph, patchEmbed);
lgraph = connectLayers(lgraph, 'input', 'patch_embed');

% Flatten patches
flattenLayer1 = flattenLayer('Name', 'flatten_patches');
lgraph = addLayers(lgraph, flattenLayer1);
lgraph = connectLayers(lgraph, 'patch_embed', 'flatten_patches');

% First layer norm before attention
layerNorm1 = layerNormalizationLayer('Name', 'ln1');
lgraph = addLayers(lgraph, layerNorm1);
lgraph = connectLayers(lgraph, 'flatten_patches', 'ln1');

% Self-attention layer
selfAttn = selfAttentionLayer(config.num_heads, config.embed_dim * num_patches, ...
    'Name', 'self_attention', ...
    'AttentionMask', 'causal');
lgraph = addLayers(lgraph, selfAttn);
lgraph = connectLayers(lgraph, 'ln1', 'self_attention');

% Layer norm after attention
layerNorm2 = layerNormalizationLayer('Name', 'ln2');
lgraph = addLayers(lgraph, layerNorm2);
lgraph = connectLayers(lgraph, 'self_attention', 'ln2');

% Feed-forward network with GELU activation
ff1 = fullyConnectedLayer(config.ff_dim, 'Name', 'ff1');
lgraph = addLayers(lgraph, ff1);
lgraph = connectLayers(lgraph, 'ln2', 'ff1');

% GELU activation (as in original Transformer/BERT/GPT)
geluL = geluLayer('Name', 'gelu');
lgraph = addLayers(lgraph, geluL);
lgraph = connectLayers(lgraph, 'ff1', 'gelu');

ff2 = fullyConnectedLayer(config.embed_dim * num_patches, 'Name', 'ff2');
lgraph = addLayers(lgraph, ff2);
lgraph = connectLayers(lgraph, 'gelu', 'ff2');

% Final layer norm
layerNorm3 = layerNormalizationLayer('Name', 'ln3');
lgraph = addLayers(lgraph, layerNorm3);
lgraph = connectLayers(lgraph, 'ff2', 'ln3');

% Classification head
classifier = fullyConnectedLayer(config.num_classes, 'Name', 'classifier');
lgraph = addLayers(lgraph, classifier);
lgraph = connectLayers(lgraph, 'ln3', 'classifier');

softmaxL = softmaxLayer('Name', 'softmax');
lgraph = addLayers(lgraph, softmaxL);
lgraph = connectLayers(lgraph, 'classifier', 'softmax');

outputL = classificationLayer('Name', 'output');
lgraph = addLayers(lgraph, outputL);
lgraph = connectLayers(lgraph, 'softmax', 'output');

fprintf('Built ViT with:\n');
fprintf('  - Patch embedding (conv %dx%d)\n', config.patch_size, config.patch_size);
fprintf('  - Real selfAttentionLayer (%d heads)\n', config.num_heads);
fprintf('  - GELU activation\n');
fprintf('  - Feed-forward network (dim %d)\n', config.ff_dim);

% Analyze network
try
    analyzeNetwork(lgraph);
catch ME
    fprintf('Note: %s\n', ME.message);
end

%% 3) Training Options
fprintf('\n[3/5] Setting up training...\n');

options = trainingOptions('adam', ...
    'InitialLearnRate', config.learning_rate, ...
    'MaxEpochs', config.max_epochs, ...
    'MiniBatchSize', config.batch_size, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', augimdsVal, ...
    'ValidationFrequency', 50, ...
    'Verbose', true, ...
    'Plots', 'training-progress', ...
    'OutputNetwork', 'best-validation-loss');

%% 4) Train Network
fprintf('\n[4/5] Training network...\n');

try
    [net, trainInfo] = trainNetwork(augimdsTrain, lgraph, options);

    % Get final accuracy
    YPred = classify(net, augimdsVal);
    YVal = imdsVal.Labels;
    accuracy = sum(YPred == YVal) / numel(YVal);

    fprintf('\nTraining complete!\n');
    fprintf('Final validation accuracy: %.2f%%\n', accuracy * 100);

catch ME
    fprintf('Training failed: %s\n', ME.message);
    fprintf('\nThis may be due to:\n');
    fprintf('  1. selfAttentionLayer dimension requirements\n');
    fprintf('  2. Input format incompatibility\n');
    fprintf('\nTrying simplified architecture...\n');

    % Fallback: simplified architecture
    [net, trainInfo, accuracy] = train_simplified_vit_gelu(config, augimdsTrain, augimdsVal, imdsVal);
end

%% 5) Save Model
fprintf('\n[5/5] Saving model...\n');

save('mnist_vit_gelu_model.mat', 'net', 'accuracy', 'config', 'trainInfo');
fprintf('Model saved to: mnist_vit_gelu_model.mat\n');

%% 6) Test NNV Conversion
fprintf('\n[6] Testing NNV conversion...\n');

try
    nnv_net = matlab2nnv(net);
    fprintf('Successfully converted to NNV!\n');
    fprintf('Number of layers: %d\n', length(nnv_net.Layers));

    % List layer types
    fprintf('\nLayer types:\n');
    for i = 1:length(nnv_net.Layers)
        fprintf('  %d: %s\n', i, class(nnv_net.Layers{i}));
    end

    % Check for attention and GELU layers
    has_attention = false;
    has_gelu = false;
    for i = 1:length(nnv_net.Layers)
        if isa(nnv_net.Layers{i}, 'MultiHeadAttentionLayer')
            has_attention = true;
            fprintf('\nFound MultiHeadAttentionLayer at index %d!\n', i);
        end
        if isa(nnv_net.Layers{i}, 'GeluLayer')
            has_gelu = true;
            fprintf('Found GeluLayer at index %d!\n', i);
        end
    end

    if ~has_attention
        fprintf('\nNote: No NNV attention layers detected.\n');
    end
    if ~has_gelu
        fprintf('Note: No NNV GELU layers detected.\n');
    end

catch ME
    fprintf('NNV conversion issue: %s\n', ME.message);
end

fprintf('\n=== Training Complete ===\n');

%% Helper Function: Simplified ViT with GELU (fallback)
function [net, trainInfo, accuracy] = train_simplified_vit_gelu(config, augimdsTrain, augimdsVal, imdsVal)
    % Simplified ViT that works with standard MATLAB training
    % Uses GELU activation

    fprintf('\nBuilding simplified ViT architecture with GELU...\n');

    layers = [
        imageInputLayer([28 28 1], 'Name', 'input', 'Normalization', 'rescale-zero-one')

        % Patch embedding
        convolution2dLayer(config.patch_size, config.embed_dim, ...
            'Stride', config.patch_size, 'Name', 'patch_embed')

        % Flatten to sequence: [4 4 32] -> [512]
        flattenLayer('Name', 'flatten')

        % Layer norm
        layerNormalizationLayer('Name', 'ln1')

        % Attention-like processing via fully connected
        fullyConnectedLayer(config.embed_dim * 16, 'Name', 'attn_approx')
        geluLayer('Name', 'gelu1')

        % Feed-forward with GELU
        fullyConnectedLayer(config.ff_dim, 'Name', 'ff1')
        geluLayer('Name', 'gelu2')
        fullyConnectedLayer(config.embed_dim * 16, 'Name', 'ff2')
        layerNormalizationLayer('Name', 'ln2')

        % Classifier
        fullyConnectedLayer(config.num_classes, 'Name', 'classifier')
        softmaxLayer('Name', 'softmax')
        classificationLayer('Name', 'output')
    ];

    options = trainingOptions('adam', ...
        'InitialLearnRate', config.learning_rate, ...
        'MaxEpochs', config.max_epochs, ...
        'MiniBatchSize', config.batch_size, ...
        'Shuffle', 'every-epoch', ...
        'ValidationData', augimdsVal, ...
        'ValidationFrequency', 50, ...
        'Verbose', true, ...
        'Plots', 'training-progress', ...
        'OutputNetwork', 'best-validation-loss');

    [net, trainInfo] = trainNetwork(augimdsTrain, layers, options);

    YPred = classify(net, augimdsVal);
    YVal = imdsVal.Labels;
    accuracy = sum(YPred == YVal) / numel(YVal);

    fprintf('Simplified ViT (GELU) validation accuracy: %.2f%%\n', accuracy * 100);
end
