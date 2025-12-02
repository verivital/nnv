% test_soundness_LstmLayer
% Tests for LSTM layer functionality
% To run: results = runtests('test_soundness_LstmLayer')

%% Test 1: LSTM layer constructor and properties
rng(42);

inputSize = 3;
numHiddenUnits = 4;

% LSTM weights dimensions:
% inputWeights: [4*numHiddenUnits x inputSize]
% recurrentWeights: [4*numHiddenUnits x numHiddenUnits]
% bias: [4*numHiddenUnits x 1]
inputWeights = randn(4*numHiddenUnits, inputSize);
recurrentWeights = randn(4*numHiddenUnits, numHiddenUnits);
bias = randn(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_test', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

% Check properties are set (Name may default to 'lstm_layer' based on implementation)
assert(L.numHiddenUnits == numHiddenUnits, 'NumHiddenUnits should match');
assert(L.InputSize == inputSize, 'InputSize should match');
assert(~isempty(L.inputWeights), 'InputWeights should be set');
assert(~isempty(L.recurrentWeights), 'RecurrentWeights should be set');

%% Test 2: LSTM layer evaluate
rng(42);

inputSize = 2;
numHiddenUnits = 3;
seqLength = 5;

inputWeights = randn(4*numHiddenUnits, inputSize) * 0.1;
recurrentWeights = randn(4*numHiddenUnits, numHiddenUnits) * 0.1;
bias = zeros(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L2 = LstmLayer('Name', 'lstm_eval', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

% Create input sequence: inputSize x seqLength
input_seq = randn(inputSize, seqLength);
output = L2.evaluate(input_seq);

assert(~isempty(output), 'LSTM evaluate should produce output');
assert(size(output, 1) == numHiddenUnits, 'Output should have numHiddenUnits rows');

%% Test 3: LSTM deterministic output
rng(42);

inputSize = 2;
numHiddenUnits = 2;

inputWeights = ones(4*numHiddenUnits, inputSize) * 0.1;
recurrentWeights = ones(4*numHiddenUnits, numHiddenUnits) * 0.1;
bias = zeros(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L3 = LstmLayer('Name', 'lstm_det', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

% Same input should give same output
input_seq = ones(inputSize, 3);
output1 = L3.evaluate(input_seq);
output2 = L3.evaluate(input_seq);

assert(max(abs(output1(:) - output2(:))) < 1e-10, ...
    'Same input should produce same output');

%% Test 4: LSTM with different sequence lengths
rng(42);

inputSize = 3;
numHiddenUnits = 4;

inputWeights = randn(4*numHiddenUnits, inputSize) * 0.1;
recurrentWeights = randn(4*numHiddenUnits, numHiddenUnits) * 0.1;
bias = randn(4*numHiddenUnits, 1) * 0.1;
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L4 = LstmLayer('Name', 'lstm_seqlen', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

% Test with different sequence lengths
for seqLen = [1, 3, 5, 10]
    input_seq = randn(inputSize, seqLen);
    output = L4.evaluate(input_seq);
    assert(size(output, 1) == numHiddenUnits, ...
        'Output should have numHiddenUnits rows for seqLen=%d', seqLen);
end

