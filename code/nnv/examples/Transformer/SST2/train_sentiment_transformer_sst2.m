%% Train Sentiment Transformer on Real SST-2 Dataset
% This script trains a Transformer-style sentiment classifier on the
% actual Stanford Sentiment Treebank (SST-2) dataset.
%
% Architecture:
%   - Bag-of-Words input (TF-IDF weighted)
%   - Feedforward Transformer-style blocks
%   - Binary classification output
%
% Author: NNV Team
% Date: November 2025

%% Configuration
config = struct();
config.vocab_size = 5000;       % Top vocabulary words
config.embed_dim = 128;         % Embedding dimension
config.hidden_dim = 256;        % Hidden layer dimension
config.num_classes = 2;         % Binary classification
config.max_epochs = 15;         % Training epochs
config.batch_size = 64;         % Batch size
config.learning_rate = 0.001;   % Initial learning rate
config.max_train_samples = 10000;  % Limit training samples for faster training
config.use_tfidf = true;        % Use TF-IDF weighting

fprintf('=== Training Sentiment Transformer on SST-2 ===\n');
fprintf('Vocabulary size: %d\n', config.vocab_size);
fprintf('Embedding dim: %d\n', config.embed_dim);
fprintf('Max training samples: %d\n', config.max_train_samples);

%% 1) Load SST-2 Dataset
fprintf('\n[1/6] Loading SST-2 dataset...\n');

[train_data, dev_data] = download_sst2();

fprintf('Raw training samples: %d\n', length(train_data.sentences));
fprintf('Raw development samples: %d\n', length(dev_data.sentences));

% Subsample training data if needed
if config.max_train_samples < length(train_data.sentences)
    rng(42);
    idx = randperm(length(train_data.sentences), config.max_train_samples);
    train_data.sentences = train_data.sentences(idx);
    train_data.labels = train_data.labels(idx);
    fprintf('Subsampled to %d training samples\n', config.max_train_samples);
end

%% 2) Build Vocabulary
fprintf('\n[2/6] Building vocabulary...\n');

% Tokenize and count words
word_counts = containers.Map();

for i = 1:length(train_data.sentences)
    words = tokenize_sentence(train_data.sentences{i});
    for j = 1:length(words)
        w = words{j};
        if isKey(word_counts, w)
            word_counts(w) = word_counts(w) + 1;
        else
            word_counts(w) = 1;
        end
    end
end

% Sort by frequency and keep top vocab_size
all_words = keys(word_counts);
counts = zeros(length(all_words), 1);
for i = 1:length(all_words)
    counts(i) = word_counts(all_words{i});
end
[~, sort_idx] = sort(counts, 'descend');
all_words = all_words(sort_idx);

% Create word-to-index mapping
vocab_size = min(config.vocab_size, length(all_words));
word2idx = containers.Map();
idx2word = cell(vocab_size, 1);

for i = 1:vocab_size
    word2idx(all_words{i}) = i;
    idx2word{i} = all_words{i};
end

fprintf('Vocabulary size: %d unique words\n', vocab_size);
fprintf('Top words: %s, %s, %s, %s, %s\n', ...
    idx2word{1}, idx2word{2}, idx2word{3}, idx2word{4}, idx2word{5});

%% 3) Compute Document Frequencies for TF-IDF
if config.use_tfidf
    fprintf('\n[3/6] Computing TF-IDF weights...\n');

    doc_freq = zeros(vocab_size, 1);
    num_docs = length(train_data.sentences);

    for i = 1:num_docs
        words = tokenize_sentence(train_data.sentences{i});
        seen_words = containers.Map();

        for j = 1:length(words)
            w = words{j};
            if isKey(word2idx, w) && ~isKey(seen_words, w)
                idx = word2idx(w);
                doc_freq(idx) = doc_freq(idx) + 1;
                seen_words(w) = 1;
            end
        end
    end

    % IDF = log(N / df)
    idf = log(num_docs ./ (doc_freq + 1));
    fprintf('IDF range: [%.2f, %.2f]\n', min(idf), max(idf));
else
    idf = ones(vocab_size, 1);
end

%% 4) Convert to Feature Vectors
fprintf('\n[4/6] Converting sentences to feature vectors...\n');

% Training features
n_train = length(train_data.sentences);
X_train = zeros(vocab_size, 1, 1, n_train, 'single');

for i = 1:n_train
    bow = sentence_to_tfidf(train_data.sentences{i}, word2idx, vocab_size, idf);
    X_train(:, 1, 1, i) = single(bow);
end

% Training labels
Y_train = categorical(train_data.labels, [0, 1], {'negative', 'positive'});

% Validation features
n_val = length(dev_data.sentences);
X_val = zeros(vocab_size, 1, 1, n_val, 'single');

for i = 1:n_val
    bow = sentence_to_tfidf(dev_data.sentences{i}, word2idx, vocab_size, idf);
    X_val(:, 1, 1, i) = single(bow);
end

Y_val = categorical(dev_data.labels, [0, 1], {'negative', 'positive'});

fprintf('Training features shape: [%d, %d]\n', vocab_size, n_train);
fprintf('Validation features shape: [%d, %d]\n', vocab_size, n_val);

%% 5) Build Network Architecture
fprintf('\n[5/6] Building network architecture...\n');

% NOTE: BatchNormalization layers removed due to known soundness issues
% in NNV's BatchNormalizationLayer.reach() method.
% See TODO_TRANSFORMER.md "Known Issues" section for details.
layers = [
    % Input: [vocab_size x 1 x 1] bag-of-words/TF-IDF features
    imageInputLayer([vocab_size 1 1], 'Name', 'input', 'Normalization', 'none')

    % Flatten to vector
    flattenLayer('Name', 'flatten')

    % Embedding projection
    fullyConnectedLayer(config.embed_dim, 'Name', 'embed_proj')
    reluLayer('Name', 'relu1')

    % Transformer-style feedforward block 1
    fullyConnectedLayer(config.hidden_dim, 'Name', 'ff1_up')
    reluLayer('Name', 'ff1_relu')
    dropoutLayer(0.2, 'Name', 'ff1_drop')
    fullyConnectedLayer(config.embed_dim, 'Name', 'ff1_down')
    reluLayer('Name', 'relu2')

    % Transformer-style feedforward block 2
    fullyConnectedLayer(config.hidden_dim, 'Name', 'ff2_up')
    reluLayer('Name', 'ff2_relu')
    dropoutLayer(0.2, 'Name', 'ff2_drop')
    fullyConnectedLayer(config.embed_dim, 'Name', 'ff2_down')
    reluLayer('Name', 'relu3')

    % Classification head
    fullyConnectedLayer(config.num_classes, 'Name', 'classifier')
    softmaxLayer('Name', 'softmax')
    classificationLayer('Name', 'output')
];

fprintf('Network created with %d layers\n', length(layers));

%% 6) Train Network
fprintf('\n[6/6] Training network...\n');

options = trainingOptions('adam', ...
    'InitialLearnRate', config.learning_rate, ...
    'MaxEpochs', config.max_epochs, ...
    'MiniBatchSize', config.batch_size, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {X_val, Y_val}, ...
    'ValidationFrequency', floor(n_train / config.batch_size), ...
    'ValidationPatience', 5, ...
    'Verbose', true, ...
    'Plots', 'training-progress', ...
    'ExecutionEnvironment', 'auto');

[net, trainInfo] = trainNetwork(X_train, Y_train, layers, options);

%% Evaluate
fprintf('\n=== Evaluation ===\n');

% Training accuracy
Y_pred_train = classify(net, X_train);
train_accuracy = sum(Y_pred_train == Y_train) / n_train;
fprintf('Training accuracy: %.2f%%\n', train_accuracy * 100);

% Validation accuracy
Y_pred_val = classify(net, X_val);
val_accuracy = sum(Y_pred_val == Y_val) / n_val;
fprintf('Validation accuracy: %.2f%%\n', val_accuracy * 100);

% Per-class accuracy
pos_mask = Y_val == 'positive';
neg_mask = Y_val == 'negative';
pos_acc = sum(Y_pred_val(pos_mask) == 'positive') / sum(pos_mask);
neg_acc = sum(Y_pred_val(neg_mask) == 'negative') / sum(neg_mask);
fprintf('Positive class accuracy: %.2f%%\n', pos_acc * 100);
fprintf('Negative class accuracy: %.2f%%\n', neg_acc * 100);

%% Save Model
fprintf('\n=== Saving Model ===\n');

all_vocab = idx2word;
accuracy = val_accuracy;

save('sentiment_transformer_sst2_model.mat', 'net', 'accuracy', 'config', ...
    'word2idx', 'all_vocab', 'idf', 'vocab_size');
fprintf('Model saved to: sentiment_transformer_sst2_model.mat\n');

%% Test NNV Compatibility
fprintf('\n=== Testing NNV Compatibility ===\n');

try
    nnv_net = matlab2nnv(net);
    fprintf('Successfully converted to NNV!\n');
    fprintf('Number of NNV layers: %d\n', length(nnv_net.Layers));

    % Quick test
    test_input = X_val(:, 1, 1, 1);
    y_matlab = predict(net, reshape(test_input, [vocab_size, 1, 1]));
    y_nnv = nnv_net.evaluate(test_input);

    fprintf('MATLAB output: [%.4f, %.4f]\n', y_matlab(1), y_matlab(2));
    fprintf('NNV output: [%.4f, %.4f]\n', y_nnv(1), y_nnv(2));
    fprintf('Max difference: %.6f\n', max(abs(y_matlab(:) - y_nnv)));

catch ME
    fprintf('NNV conversion note: %s\n', ME.message);
end

fprintf('\n=== Training Complete ===\n');

%% Helper Functions
function words = tokenize_sentence(text)
    % Simple tokenization: lowercase, split on non-alphanumeric
    text = lower(text);
    text = regexprep(text, '[^a-z0-9\s]', ' ');
    words = strsplit(strtrim(text));
    words = words(~cellfun(@isempty, words));
end

function vec = sentence_to_tfidf(text, word2idx, vocab_size, idf)
    % Convert sentence to TF-IDF vector
    words = tokenize_sentence(text);
    vec = zeros(vocab_size, 1);

    % Count term frequencies
    for i = 1:length(words)
        w = words{i};
        if isKey(word2idx, w)
            idx = word2idx(w);
            vec(idx) = vec(idx) + 1;
        end
    end

    % Apply TF-IDF weighting
    if any(vec > 0)
        % TF = raw count (could use log(1+tf) for sublinear)
        tf = vec;
        vec = tf .* idf;

        % L2 normalize
        norm_val = norm(vec);
        if norm_val > 0
            vec = vec / norm_val;
        end
    end
end
