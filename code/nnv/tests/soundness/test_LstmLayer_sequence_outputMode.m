% test_LstmLayer_sequence_outputMode
% Regression test: with OutputMode = 'sequence', star reachability must
% return one set of ranges per time step. Previously the assembly loop
% in reach_star_single_input iterated `for i = length(ht)` (instead of
% `for i = 1:length(ht)`), so the output contained only the last time
% step regardless of the input sequence length.
% To run: results = runtests('test_LstmLayer_sequence_outputMode')

%% Test 1: sequence-mode output has one column per time step
rng(42);

inputSize = 3;
numHiddenUnits = 4;
seqLength = 5;

inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
bias = 0.1*randn(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_seq_mode', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'sequence');

x = randn(inputSize, seqLength);
ep = 0.05;
IM = ImageStar(x - ep, x + ep);

R = L.reach_star_single_input(IM);

assert(isa(R, 'ImageStar'), 'sequence-mode output must be an ImageStar');
assert(R.height == numHiddenUnits, ...
    'sequence-mode output must have numHiddenUnits rows');
assert(R.width == seqLength, ...
    'sequence-mode output must have one column per time step');

[lb, ub] = R.getRanges;
assert(all(lb(:) <= ub(:) + 1e-12), ...
    'lower bounds must not exceed upper bounds');
