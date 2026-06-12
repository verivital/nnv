% test_LstmLayer_reach_dispatch
% Regression test: LstmLayer.reach() with a star-based method dispatches
% to reach_star_multipleInputs. That method was previously missing, so
% reach() raised "Unrecognized method ... 'reach_star_multipleInputs'"
% for every star-based method. 'exact-star' now fails closed with a
% clear message, since the LSTM star set is an over-approximation.
% To run: results = runtests('test_LstmLayer_reach_dispatch')

%% Test 1: reach() works for approx-star on a single input
rng(42);

inputSize = 3;
numHiddenUnits = 4;
seqLength = 5;

inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
bias = 0.1*randn(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_dispatch', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

x = randn(inputSize, seqLength);
ep = 0.05;
IM = ImageStar(x - ep, x + ep);

R1 = L.reach(IM, 'approx-star', []);
assert(isa(R1, 'ImageStar'), 'reach output must be an ImageStar');
assert(R1.height == numHiddenUnits, 'output must have numHiddenUnits rows');

% 'exact-star' must fail closed: the LSTM star set is an
% over-approximation, and NN.verify treats exact-star results as
% complete, so silently aliasing it to approx-star would be unsound
threw = false;
try
    L.reach(IM, 'exact-star', 'single');
catch
    threw = true;
end
assert(threw, "reach() must reject 'exact-star' for LstmLayer");

%% Test 2: reach() works for an array of input sets
rng(43);

inputSize = 3;
numHiddenUnits = 4;
seqLength = 4;

inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
bias = 0.1*randn(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_dispatch_multi', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

x1 = randn(inputSize, seqLength);
x2 = randn(inputSize, seqLength);
ep = 0.05;
IMs = [ImageStar(x1 - ep, x1 + ep), ImageStar(x2 - ep, x2 + ep)];

R = L.reach(IMs, 'approx-star', 'single');
assert(length(R) == 2, 'one output set per input set');
assert(all(arrayfun(@(s) isa(s, 'ImageStar'), R)), 'outputs must be ImageStars');
assert(all(arrayfun(@(s) s.height == numHiddenUnits, R)), ...
    'each output must have numHiddenUnits rows');

%% Test 3: sequence output mode through reach() has one column per step
rng(44);

inputSize = 3;
numHiddenUnits = 4;
seqLength = 5;

inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
bias = 0.1*randn(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_dispatch_seq', ...
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

R = L.reach(IM, 'approx-star', []);
assert(isa(R, 'ImageStar'), 'sequence-mode output must be an ImageStar');
assert(R.height == numHiddenUnits, 'output must have numHiddenUnits rows');
assert(R.width == seqLength, 'sequence-mode output must have one column per time step');

%% Test 4: predicate-only ImageStar input (no im_lb/im_ub) is accepted
rng(45);

inputSize = 3;
numHiddenUnits = 4;
seqLength = 4;
ep = 0.05;

inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
bias = 0.1*randn(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_dispatch_pred', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

% one shared predicate variable: input = x + a*ep, a in [-1, 1]; this
% is the form produced by upstream layers (im_lb/im_ub empty), which
% the previous size(im_lb) check rejected
x = randn(inputSize, seqLength);
V = zeros(inputSize, seqLength, 1, 2);
V(:,:,1,1) = x;
V(:,:,1,2) = ep*ones(inputSize, seqLength);
C = [1; -1];
d = [1; 1];
IM = ImageStar(V, C, d, -1, 1);

R = L.reach(IM, 'approx-star', []);
assert(isa(R, 'ImageStar'), 'reach output must be an ImageStar');
assert(R.height == numHiddenUnits, 'output must have numHiddenUnits rows');
