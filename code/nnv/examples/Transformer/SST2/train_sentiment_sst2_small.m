%% Train Sentiment Transformer on Real SST-2 with Reduced Vocabulary
% This script trains a model with smaller vocabulary for faster verification.
%
% Vocabulary: 500 words (instead of 5000)
% This enables verification with larger epsilon values.
%
% Author: NNV Team
% Date: November 2025

%% Configuration
config = struct();
config.vocab_size = 500;        % REDUCED vocabulary size
config.embed_dim = 64;          % Smaller embedding
config.hidden_dim = 128;        % Smaller hidden
config.num_classes = 2;         % Binary classification
config.max_epochs = 15;         % Training epochs
config.batch_size = 64;         % Batch size
config.learning_rate = 0.001;   % Initial learning rate
config.max_train_samples = 10000;  % Training samples
config.use_tfidf = true;        % Use TF-IDF weighting

fprintf('=== Training Sentiment Transformer (Small Vocab) ===\n');
fprintf('Vocabulary size: %d (reduced for faster verification)\n', config.vocab_size);
fprintf('Embedding dim: %d\n', config.embed_dim);

%% 1) Load SST-2 Dataset
fprintf('\n[1/6] Loading SST-2 dataset...\n');

[train_data, dev_data] = download_sst2();

% Subsample training data
if config.max_train_samples < length(train_data.sentences)
    rng(42);
    idx = randperm(length(train_data.sentences), config.max_train_samples);
    train_data.sentences = train_data.sentences(idx);
    train_data.labels = train_data.labels(idx);
end

fprintf('Training samples: %d\n', length(train_data.sentences));

%% 2) Build Vocabulary
fprintf('\n[2/6] Building vocabulary...\n');

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

% Sort by frequency
all_words = keys(word_counts);
counts = zeros(length(all_words), 1);
for i = 1:length(all_words)
    counts(i) = word_counts(all_words{i});
end
[~, sort_idx] = sort(counts, 'descend');
all_words = all_words(sort_idx);

% Keep top vocab_size words
vocab_size = min(config.vocab_size, length(all_words));
word2idx = containers.Map();
idx2word = cell(vocab_size, 1);

for i = 1:vocab_size
    word2idx(all_words{i}) = i;
    idx2word{i} = all_words{i};
end

fprintf('Vocabulary size: %d unique words\n', vocab_size);

%% 3) Compute TF-IDF
fprintf('\n[3/6] Computing TF-IDF weights...\n');

doc_freq = zeros(vocab_size, 1);
num_docs = length(train_data.sentences);

for i = 1:num_docs
    words = tokenize_sentence(train_data.sentences{i});
    seen = containers.Map();
    for j = 1:length(words)
        w = words{j};
        if isKey(word2idx, w) && ~isKey(seen, w)
            idx = word2idx(w);
            doc_freq(idx) = doc_freq(idx) + 1;
            seen(w) = 1;
        end
    end
end

idf = log(num_docs ./ (doc_freq + 1));
fprintf('IDF range: [%.2f, %.2f]\n', min(idf), max(idf));

%% 4) Convert to Features
fprintf('\n[4/6] Converting to feature vectors...\n');

n_train = length(train_data.sentences);
X_train = zeros(vocab_size, 1, 1, n_train, 'single');

for i = 1:n_train
    bow = sentence_to_tfidf(train_data.sentences{i}, word2idx, vocab_size, idf);
    X_train(:, 1, 1, i) = single(bow);
end

Y_train = categorical(train_data.labels, [0, 1], {'negative', 'positive'});

n_val = length(dev_data.sentences);
X_val = zeros(vocab_size, 1, 1, n_val, 'single');

for i = 1:n_val
    bow = sentence_to_tfidf(dev_data.sentences{i}, word2idx, vocab_size, idf);
    X_val(:, 1, 1, i) = single(bow);
end

Y_val = categorical(dev_data.labels, [0, 1], {'negative', 'positive'});

fprintf('Training: [%d, %d], Validation: [%d, %d]\n', vocab_size, n_train, vocab_size, n_val);

%% 5) Build Network
fprintf('\n[5/6] Building network...\n');

% NOTE: BatchNormalization layers removed due to known soundness issues
% in NNV's BatchNormalizationLayer.reach() method.
% See TODO_TRANSFORMER.md "Known Issues" section for details.
layers = [
    imageInputLayer([vocab_size 1 1], 'Name', 'input', 'Normalization', 'none')
    flattenLayer('Name', 'flatten')
    fullyConnectedLayer(config.embed_dim, 'Name', 'embed')
    reluLayer('Name', 'relu1')
    fullyConnectedLayer(config.hidden_dim, 'Name', 'ff1')
    reluLayer('Name', 'ff1_relu')
    dropoutLayer(0.2, 'Name', 'drop1')
    fullyConnectedLayer(config.embed_dim, 'Name', 'ff2')
    reluLayer('Name', 'relu2')
    fullyConnectedLayer(config.num_classes, 'Name', 'classifier')
    softmaxLayer('Name', 'softmax')
    classificationLayer('Name', 'output')
];

fprintf('Network: %d layers\n', length(layers));

%% 6) Train
fprintf('\n[6/6] Training...\n');

options = trainingOptions('adam', ...
    'InitialLearnRate', config.learning_rate, ...
    'MaxEpochs', config.max_epochs, ...
    'MiniBatchSize', config.batch_size, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {X_val, Y_val}, ...
    'ValidationFrequency', floor(n_train / config.batch_size), ...
    'ValidationPatience', 5, ...
    'Verbose', true, ...
    'Plots', 'training-progress');

[net, trainInfo] = trainNetwork(X_train, Y_train, layers, options);

%% Evaluate
Y_pred_val = classify(net, X_val);
val_accuracy = sum(Y_pred_val == Y_val) / n_val;
fprintf('\nValidation accuracy: %.2f%%\n', val_accuracy * 100);

%% Save
all_vocab = idx2word;
accuracy = val_accuracy;

save('sentiment_sst2_small_model.mat', 'net', 'accuracy', 'config', ...
    'word2idx', 'all_vocab', 'idf', 'vocab_size');
fprintf('Model saved to: sentiment_sst2_small_model.mat\n');

%% Test NNV
fprintf('\nTesting NNV compatibility...\n');
try
    nnv_net = matlab2nnv(net);
    fprintf('NNV conversion successful: %d layers\n', length(nnv_net.Layers));
catch ME
    fprintf('NNV note: %s\n', ME.message);
end

fprintf('\n=== Training Complete ===\n');

%% Helper Functions
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
            idx = word2idx(w);
            vec(idx) = vec(idx) + 1;
        end
    end
    if any(vec > 0)
        vec = vec .* idf;
        if norm(vec) > 0
            vec = vec / norm(vec);
        end
    end
end
