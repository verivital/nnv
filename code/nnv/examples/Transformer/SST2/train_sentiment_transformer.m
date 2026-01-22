%% Train a Transformer-based Sentiment Classifier
% This script trains a Transformer-style sentiment classifier for verification
%
% Architecture (Transformer-inspired for NLP):
%   - Bag-of-words input representation (fixed vocabulary)
%   - Learned "embedding" via FC layer
%   - Multi-layer FC network simulating attention (Q,K,V projections)
%   - Feed-forward network
%   - Binary classification head
%
% This architecture is designed to:
%   1. Work well with NNV verification
%   2. Achieve high accuracy on sentiment classification
%   3. Demonstrate transformer-style layer operations
%
% Author: NNV Team
% Date: November 2025

%% Configuration
config = struct();
config.vocab_size = 100;        % Vocabulary size (including special tokens)
config.embed_dim = 32;          % Embedding dimension
config.hidden_dim = 64;         % Hidden dimension
config.num_classes = 2;         % Binary classification
config.max_epochs = 100;        % Training epochs
config.batch_size = 32;         % Batch size
config.learning_rate = 0.001;   % Initial learning rate

fprintf('=== Training Transformer Sentiment Classifier ===\n');
fprintf('Vocabulary size: %d\n', config.vocab_size);
fprintf('Embedding dim: %d\n', config.embed_dim);

%% 1) Create Extended Dataset with Clear Sentiment Signals
fprintf('\n[1/6] Creating extended sentiment dataset...\n');

% Define clear sentiment vocabulary
positive_words = {'great', 'excellent', 'amazing', 'wonderful', 'fantastic', ...
    'brilliant', 'outstanding', 'superb', 'perfect', 'love', 'loved', ...
    'best', 'beautiful', 'impressive', 'incredible', 'marvelous', ...
    'terrific', 'awesome', 'delightful', 'enjoyable', 'masterpiece', ...
    'captivating', 'engaging', 'heartwarming', 'inspiring', 'remarkable'};

negative_words = {'terrible', 'awful', 'horrible', 'bad', 'poor', ...
    'disappointing', 'boring', 'dull', 'weak', 'hate', 'hated', ...
    'worst', 'ugly', 'unimpressive', 'mediocre', 'forgettable', ...
    'dreadful', 'annoying', 'tedious', 'painful', 'disaster', ...
    'confusing', 'predictable', 'depressing', 'uninspiring', 'bland'};

neutral_words = {'movie', 'film', 'story', 'acting', 'plot', 'direction', ...
    'performance', 'experience', 'work', 'cinematography', 'ending', ...
    'characters', 'script', 'visuals', 'music', 'scene', 'scenes'};

modifiers = {'very', 'really', 'truly', 'absolutely', 'so', 'quite', ''};

% Generate positive samples
positive_samples = {};
for i = 1:length(positive_words)
    for j = 1:min(3, length(neutral_words))
        for k = 1:length(modifiers)
            if isempty(modifiers{k})
                positive_samples{end+1} = sprintf('%s %s', positive_words{i}, neutral_words{j});
            else
                positive_samples{end+1} = sprintf('%s %s %s', modifiers{k}, positive_words{i}, neutral_words{j});
            end
        end
    end
end

% Generate negative samples
negative_samples = {};
for i = 1:length(negative_words)
    for j = 1:min(3, length(neutral_words))
        for k = 1:length(modifiers)
            if isempty(modifiers{k})
                negative_samples{end+1} = sprintf('%s %s', negative_words{i}, neutral_words{j});
            else
                negative_samples{end+1} = sprintf('%s %s %s', modifiers{k}, negative_words{i}, neutral_words{j});
            end
        end
    end
end

% Combine and balance
n_samples = min(length(positive_samples), length(negative_samples));
positive_samples = positive_samples(1:n_samples);
negative_samples = negative_samples(1:n_samples);

all_samples = [positive_samples'; negative_samples'];
labels = [ones(n_samples, 1); 2*ones(n_samples, 1)];
labels = categorical(labels, [1 2], {'positive', 'negative'});

% Shuffle
rng(42);
idx = randperm(length(all_samples));
all_samples = all_samples(idx);
labels = labels(idx);

% Split train/val (80/20)
n_train = round(0.8 * length(all_samples));
train_samples = all_samples(1:n_train);
train_labels = labels(1:n_train);
val_samples = all_samples(n_train+1:end);
val_labels = labels(n_train+1:end);

fprintf('Total samples: %d\n', length(all_samples));
fprintf('Training samples: %d\n', n_train);
fprintf('Validation samples: %d\n', length(val_samples));

%% 2) Build Vocabulary and Create Bag-of-Words Representation
fprintf('\n[2/6] Building vocabulary and creating input representation...\n');

% Build vocabulary from all sentiment and neutral words
all_vocab = [positive_words, negative_words, neutral_words, modifiers];
all_vocab = unique(lower(all_vocab));
all_vocab = all_vocab(~cellfun(@isempty, all_vocab));  % Remove empty strings
vocab_size = length(all_vocab);

% Create word-to-index mapping
word2idx = containers.Map();
for i = 1:vocab_size
    word2idx(all_vocab{i}) = i;
end

fprintf('Vocabulary size: %d words\n', vocab_size);

% Create bag-of-words representation
% Each sample is a vector of size [vocab_size x 1] with word presence/counts
X_train = zeros(vocab_size, 1, 1, n_train, 'single');
for i = 1:n_train
    bow = text_to_bow(train_samples{i}, word2idx, vocab_size);
    X_train(:, 1, 1, i) = single(bow);
end

X_val = zeros(vocab_size, 1, 1, length(val_samples), 'single');
for i = 1:length(val_samples)
    bow = text_to_bow(val_samples{i}, word2idx, vocab_size);
    X_val(:, 1, 1, i) = single(bow);
end

fprintf('Training data shape: %s\n', mat2str(size(X_train)));
fprintf('Validation data shape: %s\n', mat2str(size(X_val)));

%% 3) Build Transformer-style Architecture
fprintf('\n[3/6] Building Transformer-style architecture...\n');

% Architecture:
% 1. Input: Bag-of-words vector [vocab_size x 1 x 1]
% 2. "Embedding" projection: vocab_size -> embed_dim
% 3. Transformer-style processing (Q,K,V projections via FC)
% 4. Feed-forward network
% 5. Classification head

% NOTE: BatchNormalization layers removed due to known soundness issues
% in NNV's BatchNormalizationLayer.reach() method.
% See TODO_TRANSFORMER.md "Known Issues" section for details.
layers = [
    % Input layer - bag-of-words representation
    imageInputLayer([vocab_size 1 1], 'Name', 'input', 'Normalization', 'none')

    % Flatten to vector
    flattenLayer('Name', 'flatten')

    % "Embedding" projection - maps BoW to dense representation
    fullyConnectedLayer(config.embed_dim, 'Name', 'embed_proj')
    reluLayer('Name', 'relu1')

    % Transformer-style block 1: Q,K,V projections
    fullyConnectedLayer(config.hidden_dim * 3, 'Name', 'qkv_proj1')
    reluLayer('Name', 'relu_qkv1')
    fullyConnectedLayer(config.hidden_dim, 'Name', 'attn_out1')
    reluLayer('Name', 'relu_attn1')

    % Feed-forward block 1
    fullyConnectedLayer(config.hidden_dim * 2, 'Name', 'ff1')
    reluLayer('Name', 'relu_ff1')
    fullyConnectedLayer(config.hidden_dim, 'Name', 'ff1_out')
    reluLayer('Name', 'relu_ff1_out')

    % Transformer-style block 2: Another layer of processing
    fullyConnectedLayer(config.hidden_dim * 2, 'Name', 'qkv_proj2')
    reluLayer('Name', 'relu_qkv2')
    fullyConnectedLayer(config.hidden_dim, 'Name', 'attn_out2')
    reluLayer('Name', 'relu_attn2')

    % Final projection before classification
    fullyConnectedLayer(config.embed_dim, 'Name', 'final_proj')
    reluLayer('Name', 'relu_final')

    % Classification head
    fullyConnectedLayer(config.num_classes, 'Name', 'classifier')
    softmaxLayer('Name', 'softmax')
    classificationLayer('Name', 'output')
];

fprintf('Built Transformer-style architecture:\n');
fprintf('  - Input: Bag-of-words (%d dims)\n', vocab_size);
fprintf('  - Embedding projection (%d dims)\n', config.embed_dim);
fprintf('  - 2 Transformer-style blocks (FC attention simulation)\n');
fprintf('  - Binary classification head\n');

% Analyze network
analyzeNetwork(layers);

%% 4) Training Options
fprintf('\n[4/6] Setting up training...\n');

options = trainingOptions('adam', ...
    'InitialLearnRate', config.learning_rate, ...
    'MaxEpochs', config.max_epochs, ...
    'MiniBatchSize', config.batch_size, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {X_val, val_labels}, ...
    'ValidationFrequency', 50, ...
    'Verbose', true, ...
    'Plots', 'training-progress', ...
    'OutputNetwork', 'best-validation-loss');

%% 5) Train Network
fprintf('\n[5/6] Training network...\n');

try
    [net, trainInfo] = trainNetwork(X_train, train_labels, layers, options);

    % Get final accuracy
    YPred = classify(net, X_val);
    accuracy = sum(YPred == val_labels) / numel(val_labels);

    fprintf('\nTraining complete!\n');
    fprintf('Final validation accuracy: %.2f%%\n', accuracy * 100);

catch ME
    fprintf('Training failed: %s\n', ME.message);
    rethrow(ME);
end

%% 6) Save Model
fprintf('\n[6/6] Saving model...\n');

% Save the trained network and vocabulary
save('sentiment_transformer_model.mat', 'net', 'accuracy', 'config', 'trainInfo', 'word2idx', 'all_vocab');
fprintf('Model saved to: sentiment_transformer_model.mat\n');

% Export to ONNX
try
    exportONNXNetwork(net, 'sentiment_transformer_model.onnx');
    fprintf('ONNX model saved to: sentiment_transformer_model.onnx\n');
catch ME
    fprintf('ONNX export note: %s\n', ME.message);
end

%% 7) Test with NNV
fprintf('\n[7] Testing NNV compatibility...\n');

try
    % Convert to NNV
    nnv_net = matlab2nnv(net);
    fprintf('Successfully converted to NNV!\n');
    fprintf('Number of layers: %d\n', length(nnv_net.Layers));

    % Test evaluation on multiple samples
    fprintf('\nTesting NNV vs MATLAB outputs:\n');
    n_test = min(5, size(X_val, 4));
    max_diff = 0;

    for i = 1:n_test
        test_input = X_val(:,:,:,i);
        y_matlab = predict(net, test_input);
        y_nnv = nnv_net.evaluate(test_input);

        diff = max(abs(y_matlab(:) - y_nnv(:)));
        max_diff = max(max_diff, diff);

        fprintf('  Sample %d: MATLAB [%.4f, %.4f], NNV [%.4f, %.4f], diff=%.6f\n', ...
            i, y_matlab(1), y_matlab(2), y_nnv(1), y_nnv(2), diff);
    end

    if max_diff < 1e-4
        fprintf('All outputs match (max diff: %.6f)\n', max_diff);
    else
        fprintf('Warning: Max output difference: %.6f\n', max_diff);
    end

catch ME
    fprintf('NNV conversion note: %s\n', ME.message);
end

%% 8) Quick verification test
fprintf('\n[8] Quick verification test...\n');

try
    % Test on a known positive sample
    test_text = 'great movie';
    test_bow = text_to_bow(test_text, word2idx, vocab_size);
    test_input = reshape(single(test_bow), [vocab_size, 1, 1]);

    y_pred = predict(net, test_input);
    [~, pred_class] = max(y_pred);

    fprintf('Test: "%s"\n', test_text);
    fprintf('  Prediction: %s (%.2f%% positive, %.2f%% negative)\n', ...
        categorical([pred_class], [1 2], {'positive', 'negative'}), ...
        y_pred(1)*100, y_pred(2)*100);

    % Test on a known negative sample
    test_text2 = 'terrible film';
    test_bow2 = text_to_bow(test_text2, word2idx, vocab_size);
    test_input2 = reshape(single(test_bow2), [vocab_size, 1, 1]);

    y_pred2 = predict(net, test_input2);
    [~, pred_class2] = max(y_pred2);

    fprintf('Test: "%s"\n', test_text2);
    fprintf('  Prediction: %s (%.2f%% positive, %.2f%% negative)\n', ...
        categorical([pred_class2], [1 2], {'positive', 'negative'}), ...
        y_pred2(1)*100, y_pred2(2)*100);

catch ME
    fprintf('Verification test note: %s\n', ME.message);
end

fprintf('\n=== Training Complete ===\n');

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

    % Normalize by L2 norm (or leave as counts)
    if norm(bow) > 0
        bow = bow / norm(bow);
    end
end
