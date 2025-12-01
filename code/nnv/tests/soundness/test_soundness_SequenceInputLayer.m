% test_soundness_SequenceInputLayer
% Tests for Sequence Input layer (for RNN/LSTM inputs)
% To run: results = runtests('test_soundness_SequenceInputLayer')

%% Test 1: SequenceInputLayer constructor with name
rng(42);

L = SequenceInputLayer('Name', 'seq_test', 'InputSize', 10);

assert(strcmp(L.Name, 'seq_test'), 'Name should be set');
assert(L.InputSize == 10, 'InputSize should be 10');

%% Test 2: SequenceInputLayer constructor with normalization
rng(42);

L2 = SequenceInputLayer('Name', 'seq_norm', 'InputSize', 5, 'Normalization', 'none');

assert(L2.InputSize == 5, 'InputSize should be 5');
assert(strcmp(L2.Normalization, 'none'), 'Normalization should be none');

%% Test 3: SequenceInputLayer evaluate with no normalization
rng(42);

L3 = SequenceInputLayer('Name', 'seq_eval', 'InputSize', 4, 'Normalization', 'none');

% Sequence input: features x timesteps
input = rand(4, 10);  % 4 features, 10 timesteps
output = L3.evaluate(input);

% With no normalization, output should equal input
assert(max(abs(output(:) - input(:))) < 1e-10, 'No normalization should preserve input');

%% Test 4: SequenceInputLayer evaluateSequence
rng(42);

L4 = SequenceInputLayer('Name', 'seq_evalseq', 'InputSize', 3, 'Normalization', 'none');

input = rand(3, 5);  % 3 features, 5 timesteps
output = L4.evaluateSequence(input);

assert(max(abs(output(:) - input(:))) < 1e-10, 'evaluateSequence should work');

%% Test 5: SequenceInputLayer with longer sequence
% NOTE: SequenceInputLayer uses reachSequence, not reach
rng(42);

L5 = SequenceInputLayer('Name', 'seq_long', 'InputSize', 4, 'Normalization', 'none');

% Longer sequence
input = rand(4, 20);  % 4 features, 20 timesteps
output = L5.evaluate(input);

assert(isequal(size(output), size(input)), 'Output shape should match input');
assert(max(abs(output(:) - input(:))) < 1e-10, 'Should preserve input');

