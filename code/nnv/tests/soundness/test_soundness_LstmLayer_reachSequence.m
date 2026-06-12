% test_soundness_LstmLayer_reachSequence
% Soundness (sandwich) tests for LstmLayer.reachSequence: concrete
% trajectories evaluated from inside the input set must lie inside the
% computed reachable set. The previous implementation of
% reach_star_single_input_Sequence propagated the im_lb/im_ub bound
% images concretely through evaluateSequence (or mapped each basis
% vector of V through the nonlinear network), and sampled trajectories
% routinely left the returned set.
% To run: results = runtests('test_soundness_LstmLayer_reachSequence')

%% Test 1: reachSequence encloses concrete trajectories (bounded input)
rng(42);

inputSize = 3;
numHiddenUnits = 4;
seqLength = 5;
ep = 0.05;
tol = 1e-6;

inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
bias = 0.1*randn(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_seq_sound', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

x = randn(inputSize, seqLength);
IM = ImageStar(x - ep, x + ep);

R = L.reachSequence(IM, 'approx-star');
assert(isa(R, 'ImageStar'), 'reachSequence output must be an ImageStar');
[lb, ub] = R.getRanges;

% bounds of the output star along random directions (catches
% unsoundness that per-dimension ranges cannot)
nDirs = 5;
W = randn(nDirs, numHiddenUnits);
S = R.toStar;
wlb = zeros(nDirs, 1);
wub = zeros(nDirs, 1);
for k = 1:nDirs
    P = S.affineMap(W(k, :), []);
    [wlb(k), wub(k)] = P.getRanges;
end

% center, bound images, and random samples must all be enclosed
checkpoints = {x, x - ep, x + ep};
for k = 1:20
    checkpoints{end+1} = (x - ep) + 2*ep.*rand(inputSize, seqLength); %#ok<SAGROW>
end
for k = 1:length(checkpoints)
    y = L.evaluateSequence(checkpoints{k});
    assert(all(y(:) >= lb(:) - tol) && all(y(:) <= ub(:) + tol), ...
        sprintf('concrete trajectory %d must be enclosed by reachSequence output', k));
    wy = W*y(:);
    assert(all(wy >= wlb - tol) && all(wy <= wub + tol), ...
        sprintf('concrete trajectory %d must be enclosed in projections', k));
end

% 'exact-star' is accepted for backward compatibility (with a warning)
% and dispatches through the same over-approximating star path
warning('off', 'NNV:LstmLayer:approximateReach');
R2 = L.reachSequence(IM, 'exact-star');
warning('on', 'NNV:LstmLayer:approximateReach');
[lb2, ub2] = R2.getRanges;
y0 = L.evaluateSequence(x);
assert(all(y0(:) >= lb2(:) - tol) && all(y0(:) <= ub2(:) + tol), ...
    'center trajectory must be enclosed for exact-star as well');

%% Test 2: predicate-only ImageStar input (no im_lb/im_ub) is enclosed
rng(43);

inputSize = 3;
numHiddenUnits = 4;
seqLength = 4;
ep = 0.05;
tol = 1e-6;

inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
bias = 0.1*randn(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_seq_sound_pred', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

% one shared predicate variable: input = x + a*ep, a in [-1, 1]
x = randn(inputSize, seqLength);
V = zeros(inputSize, seqLength, 1, 2);
V(:,:,1,1) = x;
V(:,:,1,2) = ep*ones(inputSize, seqLength);
C = [1; -1];
d = [1; 1];
IM = ImageStar(V, C, d, -1, 1);

R = L.reachSequence(IM, 'approx-star');
assert(isa(R, 'ImageStar'), 'reachSequence output must be an ImageStar');
[lb, ub] = R.getRanges;

for k = 1:20
    a = -1 + 2*rand;
    y = L.evaluateSequence(x + a*ep);
    assert(all(y(:) >= lb(:) - tol) && all(y(:) <= ub(:) + tol), ...
        sprintf('trajectory for predicate sample %d must be enclosed', k));
end

%% Test 3: sequence output mode encloses every time step
rng(44);

inputSize = 3;
numHiddenUnits = 4;
seqLength = 5;
ep = 0.05;
tol = 1e-6;

inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
bias = 0.1*randn(4*numHiddenUnits, 1);
cellState = zeros(numHiddenUnits, 1);
hiddenState = zeros(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_seq_sound_seqmode', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'sequence');

x = randn(inputSize, seqLength);
IM = ImageStar(x - ep, x + ep);

R = L.reachSequence(IM, 'approx-star');
assert(isa(R, 'ImageStar'), 'reachSequence output must be an ImageStar');
assert(R.height == numHiddenUnits, 'output must have numHiddenUnits rows');
assert(R.width == seqLength, ...
    'sequence-mode output must have one column per time step');

[lb, ub] = R.getRanges;
lb = reshape(lb, numHiddenUnits, seqLength);
ub = reshape(ub, numHiddenUnits, seqLength);

% per-step hidden states of sampled trajectories must be enclosed
% column-by-column (h_t equals the final hidden state of the prefix
% x(:,1:t), since the initial states are fixed)
for k = 1:10
    xs = (x - ep) + 2*ep.*rand(inputSize, seqLength);
    for t = 1:seqLength
        ht = L.evaluateSequence(xs(:, 1:t));
        assert(all(ht(:) >= lb(:, t) - tol) && all(ht(:) <= ub(:, t) + tol), ...
            sprintf('sample %d, step %d must be enclosed by column %d', k, t, t));
    end
end

%% Test 4: nonzero initial hidden/cell states are handled soundly
rng(45);

inputSize = 3;
numHiddenUnits = 4;
seqLength = 4;
ep = 0.05;
tol = 1e-6;

inputWeights = 0.5*randn(4*numHiddenUnits, inputSize);
recurrentWeights = 0.5*randn(4*numHiddenUnits, numHiddenUnits);
bias = 0.1*randn(4*numHiddenUnits, 1);
cellState = 0.2*randn(numHiddenUnits, 1);
hiddenState = 0.2*randn(numHiddenUnits, 1);

L = LstmLayer('Name', 'lstm_seq_sound_state', ...
    'InputSize', inputSize, ...
    'NumHiddenUnits', numHiddenUnits, ...
    'InputWeights', inputWeights, ...
    'RecurrentWeights', recurrentWeights, ...
    'Bias', bias, ...
    'CellState', cellState, ...
    'HiddenState', hiddenState, ...
    'OutputMode', 'last');

x = randn(inputSize, seqLength);
IM = ImageStar(x - ep, x + ep);

R = L.reachSequence(IM, 'approx-star');
[lb, ub] = R.getRanges;

for k = 1:20
    xs = (x - ep) + 2*ep.*rand(inputSize, seqLength);
    y = L.evaluateSequence(xs);
    assert(all(y(:) >= lb(:) - tol) && all(y(:) <= ub(:) + tol), ...
        sprintf('trajectory %d must be enclosed with nonzero initial states', k));
end
