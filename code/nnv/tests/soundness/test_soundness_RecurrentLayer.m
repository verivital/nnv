% test_soundness_RecurrentLayer
% Tests for vanilla Recurrent layer
% To run: results = runtests('test_soundness_RecurrentLayer')

%% Test 1: RecurrentLayer constructor
rng(42);

nH = 4;  % hidden units
nI = 3;  % input size
nO = 2;  % output size

% Create struct with required fields
params = struct();
params.Wh = randn(nH, nH) * 0.1;  % hidden-to-hidden
params.bh = randn(nH, 1) * 0.1;   % hidden bias
params.fh = 'poslin';              % ReLU activation for hidden

params.Wo = randn(nO, nH) * 0.1;  % hidden-to-output
params.bo = randn(nO, 1) * 0.1;   % output bias
params.fo = 'poslin';              % ReLU activation for output

params.Wi = randn(nH, nI) * 0.1;  % input-to-hidden

L = RecurrentLayer(params);

assert(L.nH == nH, 'Hidden units should match');
assert(L.nI == nI, 'Input size should match');
assert(L.nO == nO, 'Output size should match');

%% Test 2: RecurrentLayer evaluate
rng(42);

nH = 3;
nI = 2;
nO = 2;

params = struct();
params.Wh = randn(nH, nH) * 0.1;
params.bh = zeros(nH, 1);
params.fh = 'poslin';

params.Wo = randn(nO, nH) * 0.1;
params.bo = zeros(nO, 1);
params.fo = 'poslin';

params.Wi = randn(nH, nI) * 0.1;

L2 = RecurrentLayer(params);

% Create input sequence: [nI x sequence_length]
seq_len = 5;
input_seq = randn(nI, seq_len);

output = L2.evaluate(input_seq);

assert(size(output, 1) == nO, 'Output dimension should match nO');
assert(size(output, 2) == seq_len, 'Output sequence length should match input');

%% Test 3: RecurrentLayer with LeakyReLU activation
rng(42);

nH = 3;
nI = 2;
nO = 2;

params = struct();
params.Wh = randn(nH, nH) * 0.1;
params.bh = zeros(nH, 1);
params.fh = 'leakyrelu';
params.gh = 0.01;  % leaky factor

params.Wo = randn(nO, nH) * 0.1;
params.bo = zeros(nO, 1);
params.fo = 'leakyrelu';
params.go = 0.01;

params.Wi = randn(nH, nI) * 0.1;

L3 = RecurrentLayer(params);

assert(L3.gh == 0.01, 'Hidden leaky factor should be set');
assert(L3.go == 0.01, 'Output leaky factor should be set');

%% Test 4: RecurrentLayer properties
rng(42);

nH = 4;
nI = 3;
nO = 2;

params = struct();
params.Wh = randn(nH, nH);
params.bh = randn(nH, 1);
params.fh = 'poslin';

params.Wo = randn(nO, nH);
params.bo = randn(nO, 1);
params.fo = 'poslin';

params.Wi = randn(nH, nI);

L4 = RecurrentLayer(params);

assert(size(L4.Wh, 1) == nH, 'Wh rows should equal nH');
assert(size(L4.Wh, 2) == nH, 'Wh cols should equal nH');
assert(size(L4.Wi, 2) == nI, 'Wi cols should equal nI');
assert(size(L4.Wo, 1) == nO, 'Wo rows should equal nO');

