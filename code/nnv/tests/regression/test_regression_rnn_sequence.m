% test_regression_rnn_sequence
% Regression tests for RNN sequence verification
% Tests RecurrentLayer and LstmLayer classes
% To run: results = runtests('test_regression_rnn_sequence')

%% Test 1: RecurrentLayer class exists
rng(42);
assert(exist('RecurrentLayer', 'class') == 8, 'RecurrentLayer class should exist');

%% Test 2: LstmLayer class exists
rng(42);
assert(exist('LstmLayer', 'class') == 8, 'LstmLayer class should exist');

%% Test 3: SequenceInputLayer class exists
rng(42);
assert(exist('SequenceInputLayer', 'class') == 8, 'SequenceInputLayer class should exist');

%% Test 4: RecurrentLayer methods via metaclass
rng(42);
mc = ?RecurrentLayer;
methods = {mc.MethodList.Name};
assert(ismember('evaluate', methods), 'RecurrentLayer should have evaluate method');
assert(ismember('reach', methods), 'RecurrentLayer should have reach method');

%% Test 5: LstmLayer methods via metaclass
rng(42);
mc = ?LstmLayer;
methods = {mc.MethodList.Name};
assert(ismember('evaluate', methods), 'LstmLayer should have evaluate method');
assert(ismember('reach', methods), 'LstmLayer should have reach method');

%% Test 6: Create simple RecurrentLayer
rng(42);
rnn.Wi = 0.1 * randn(4, 3);
rnn.Wh = 0.1 * randn(4, 4);
rnn.bh = zeros(4, 1);
rnn.fh = 'poslin';
rnn.Wo = 0.1 * randn(2, 4);
rnn.bo = zeros(2, 1);
rnn.fo = 'purelin';
layer = RecurrentLayer(rnn);
assert(~isempty(layer), 'RecurrentLayer should be created');

%% Test 7: RecurrentLayer evaluate works
rng(42);
rnn.Wi = 0.1 * randn(4, 3);
rnn.Wh = 0.1 * randn(4, 4);
rnn.bh = zeros(4, 1);
rnn.fh = 'poslin';
rnn.Wo = 0.1 * randn(2, 4);
rnn.bo = zeros(2, 1);
rnn.fo = 'purelin';
layer = RecurrentLayer(rnn);
x = randn(3, 5);
y = layer.evaluate(x);
assert(size(y, 1) == 2, 'Output should have 2 rows');
assert(size(y, 2) == 5, 'Output should have 5 columns');

%% Test 8: RecurrentLayer reachability
rng(42);
rnn.Wi = 0.1 * randn(4, 3);
rnn.Wh = 0.1 * randn(4, 4);
rnn.bh = zeros(4, 1);
rnn.fh = 'poslin';
rnn.Wo = 0.1 * randn(2, 4);
rnn.bo = zeros(2, 1);
rnn.fo = 'purelin';
layer = RecurrentLayer(rnn);
lb = -0.1 * ones(3, 1); ub = 0.1 * ones(3, 1);
I1 = Star(lb, ub);
I2 = Star(lb, ub);
R = layer.reach([I1 I2], 'approx-star');
assert(~isempty(R), 'Reachable set should be computed');
