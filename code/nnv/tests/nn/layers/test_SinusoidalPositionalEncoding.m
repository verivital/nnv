%% Test Suite for SinusoidalPositionalEncoding
% Tests for fixed sinusoidal positional encoding layer
%
% PE(pos, 2i) = sin(pos / 10000^(2i/d_model))
% PE(pos, 2i+1) = cos(pos / 10000^(2i/d_model))
%
% Author: NNV Team
% Date: November 2025

%% Test 1: Basic Construction
fprintf('=== Test 1: Basic Construction ===\n');

embed_dim = 64;
layer = SinusoidalPositionalEncoding(embed_dim);

assert(layer.EmbedDim == embed_dim, 'EmbedDim not set correctly');
assert(layer.MaxSeqLen == 512, 'Default MaxSeqLen should be 512');
assert(size(layer.PE, 1) == 512, 'PE matrix rows should be MaxSeqLen');
assert(size(layer.PE, 2) == embed_dim, 'PE matrix cols should be EmbedDim');

fprintf('Test 1 PASSED: Basic construction\n');

%% Test 2: Construction with Custom MaxSeqLen
fprintf('\n=== Test 2: Construction with Custom MaxSeqLen ===\n');

layer2 = SinusoidalPositionalEncoding(32, 128);

assert(layer2.EmbedDim == 32, 'EmbedDim not set');
assert(layer2.MaxSeqLen == 128, 'MaxSeqLen not set');
assert(size(layer2.PE, 1) == 128, 'PE rows mismatch');

fprintf('Test 2 PASSED: Custom MaxSeqLen\n');

%% Test 3: Encoding Values at Position 0
fprintf('\n=== Test 3: Encoding Values at Position 0 ===\n');

layer3 = SinusoidalPositionalEncoding(4);
pe_0 = layer3.get_encoding(0);

% At position 0: sin(0) = 0 for even, cos(0) = 1 for odd
% PE(0, 0) = sin(0) = 0
% PE(0, 1) = cos(0) = 1
% PE(0, 2) = sin(0) = 0
% PE(0, 3) = cos(0) = 1

expected = [0; 1; 0; 1];
assert(all(abs(pe_0 - expected) < 1e-10), 'Position 0 encoding incorrect');

fprintf('Test 3 PASSED: Position 0 encoding\n');

%% Test 4: Encoding Uniqueness
fprintf('\n=== Test 4: Encoding Uniqueness ===\n');

layer4 = SinusoidalPositionalEncoding(64, 100);

% Each position should have a unique encoding
pe_all = layer4.PE;
for i = 1:99
    for j = i+1:100
        diff = norm(pe_all(i, :) - pe_all(j, :));
        assert(diff > 0.01, 'Positions %d and %d have same encoding', i, j);
    end
end

fprintf('Test 4 PASSED: All positions have unique encodings\n');

%% Test 5: Evaluation - Add to Input
fprintf('\n=== Test 5: Evaluation - Add to Input ===\n');

layer5 = SinusoidalPositionalEncoding(8);

x = ones(8, 1);  % Input embedding
y = layer5.evaluate(x, 0);  % Add PE for position 0

pe_0 = layer5.get_encoding(0);
expected = x + pe_0;

assert(all(abs(y - expected) < 1e-10), 'Evaluation failed');

fprintf('Test 5 PASSED: Evaluation adds PE correctly\n');

%% Test 6: Evaluation - Multiple Positions
fprintf('\n=== Test 6: Evaluation - Multiple Positions ===\n');

layer6 = SinusoidalPositionalEncoding(16);

x6 = randn(16, 4);  % [EmbedDim x SeqLen]
positions = [0, 1, 2, 3];
y6 = layer6.evaluate(x6, positions);

% Manual check
pe_expected = layer6.PE(positions + 1, :)';  % [EmbedDim x SeqLen]
expected6 = x6 + pe_expected;

assert(all(abs(y6(:) - expected6(:)) < 1e-10), 'Multi-position evaluation failed');

fprintf('Test 6 PASSED: Multiple position evaluation\n');

%% Test 7: Reachability - Star Set
fprintf('\n=== Test 7: Reachability - Star Set ===\n');

layer7 = SinusoidalPositionalEncoding(4);

lb7 = [-1; -1; -1; -1];
ub7 = [1; 1; 1; 1];
S_in = Star(lb7, ub7);

S_out = layer7.reach(S_in, 'approx-star', 0);

assert(~isempty(S_out), 'Reachability returned empty');
assert(isa(S_out, 'Star'), 'Output should be Star');

fprintf('Test 7 PASSED: Star reachability\n');

%% Test 8: Reachability - Soundness Check
fprintf('\n=== Test 8: Reachability - Soundness Check ===\n');

layer8 = SinusoidalPositionalEncoding(4);

lb8 = [-0.5; -0.5; -0.5; -0.5];
ub8 = [0.5; 0.5; 0.5; 0.5];
S_in8 = Star(lb8, ub8);

S_out8 = layer8.reach(S_in8, 'approx-star', 0);
[out_lb, out_ub] = S_out8.getRanges;

pe_0 = layer8.get_encoding(0);

% Sample and verify
n_samples = 50;
all_contained = true;

for i = 1:n_samples
    x_sample = lb8 + (ub8 - lb8) .* rand(4, 1);
    y_sample = x_sample + pe_0;

    tol = 1e-6;
    if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
        all_contained = false;
    end
end

assert(all_contained, 'Soundness violation');
fprintf('Test 8 PASSED: Soundness check (%d samples)\n', n_samples);

%% Test 9: Encoding Smoothness
fprintf('\n=== Test 9: Encoding Smoothness ===\n');

% Adjacent positions should have similar encodings
layer9 = SinusoidalPositionalEncoding(64, 100);

for pos = 0:98
    pe_curr = layer9.get_encoding(pos);
    pe_next = layer9.get_encoding(pos + 1);
    diff = norm(pe_curr - pe_next);

    % Difference should be bounded (smooth transition)
    assert(diff < 2.0, 'Encoding not smooth at position %d', pos);
end

fprintf('Test 9 PASSED: Encoding smoothness\n');

%% Test 10: Named Constructor
fprintf('\n=== Test 10: Named Constructor ===\n');

layer10 = SinusoidalPositionalEncoding(32, 64, 10000, 'my_pe');
assert(strcmp(layer10.Name, 'my_pe'), 'Name not set');

fprintf('Test 10 PASSED: Named constructor\n');

%% Test 11: Odd Embedding Dimension
fprintf('\n=== Test 11: Odd Embedding Dimension ===\n');

layer11 = SinusoidalPositionalEncoding(5);  % Odd dimension

pe_0_odd = layer11.get_encoding(0);
assert(length(pe_0_odd) == 5, 'Odd dimension PE wrong length');
assert(~any(isnan(pe_0_odd)), 'NaN in odd dimension PE');

fprintf('Test 11 PASSED: Odd embedding dimension handled\n');

%% Test 12: Large Sequence Length
fprintf('\n=== Test 12: Large Sequence Length ===\n');

layer12 = SinusoidalPositionalEncoding(64, 2048);

pe_last = layer12.get_encoding(2047);
assert(~any(isnan(pe_last)), 'NaN in large position encoding');
assert(~any(isinf(pe_last)), 'Inf in large position encoding');

fprintf('Test 12 PASSED: Large sequence length\n');

%% Summary
fprintf('\n=== All SinusoidalPositionalEncoding Tests PASSED ===\n');
fprintf('Total: 12 tests\n');
