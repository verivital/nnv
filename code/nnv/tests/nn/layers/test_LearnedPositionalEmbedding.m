%% Test Suite for LearnedPositionalEmbedding
% Tests for learned positional embedding layer
%
% Author: NNV Team
% Date: November 2025

%% Test 1: Basic Construction with Dimensions
fprintf('=== Test 1: Basic Construction with Dimensions ===\n');

embed_dim = 64;
max_seq_len = 128;
layer = LearnedPositionalEmbedding(embed_dim, max_seq_len);

assert(layer.EmbedDim == embed_dim, 'EmbedDim not set correctly');
assert(layer.MaxSeqLen == max_seq_len, 'MaxSeqLen not set correctly');
assert(size(layer.PE, 1) == max_seq_len, 'PE rows should be MaxSeqLen');
assert(size(layer.PE, 2) == embed_dim, 'PE cols should be EmbedDim');

fprintf('Test 1 PASSED: Basic construction with dimensions\n');

%% Test 2: Construction with PE Matrix
fprintf('\n=== Test 2: Construction with PE Matrix ===\n');

PE_matrix = randn(100, 32);
layer2 = LearnedPositionalEmbedding(PE_matrix);

assert(layer2.MaxSeqLen == 100, 'MaxSeqLen from matrix wrong');
assert(layer2.EmbedDim == 32, 'EmbedDim from matrix wrong');
assert(all(all(layer2.PE == PE_matrix)), 'PE matrix not set correctly');

fprintf('Test 2 PASSED: Construction with PE matrix\n');

%% Test 3: Named Constructor
fprintf('\n=== Test 3: Named Constructor ===\n');

layer3 = LearnedPositionalEmbedding(32, 64, 'my_pos_embed');
assert(strcmp(layer3.Name, 'my_pos_embed'), 'Name not set correctly');

PE_named = randn(50, 16);
layer3b = LearnedPositionalEmbedding(PE_named, 'pretrained_pe');
assert(strcmp(layer3b.Name, 'pretrained_pe'), 'Name from matrix constructor wrong');

fprintf('Test 3 PASSED: Named constructor\n');

%% Test 4: Get Encoding
fprintf('\n=== Test 4: Get Encoding ===\n');

PE4 = randn(10, 4);
layer4 = LearnedPositionalEmbedding(PE4);

pe_0 = layer4.get_encoding(0);
expected_0 = PE4(1, :)';  % Position 0 -> row 1
assert(all(abs(pe_0 - expected_0) < 1e-10), 'Position 0 encoding incorrect');

pe_5 = layer4.get_encoding(5);
expected_5 = PE4(6, :)';  % Position 5 -> row 6
assert(all(abs(pe_5 - expected_5) < 1e-10), 'Position 5 encoding incorrect');

fprintf('Test 4 PASSED: Get encoding\n');

%% Test 5: Evaluate - Single Position
fprintf('\n=== Test 5: Evaluate - Single Position ===\n');

PE5 = randn(20, 8);
layer5 = LearnedPositionalEmbedding(PE5);

x5 = randn(8, 1);
y5 = layer5.evaluate(x5, 3);

expected5 = x5 + PE5(4, :)';  % Position 3 -> row 4
assert(all(abs(y5 - expected5) < 1e-10), 'Single position evaluation failed');

fprintf('Test 5 PASSED: Single position evaluation\n');

%% Test 6: Evaluate - Multiple Positions
fprintf('\n=== Test 6: Evaluate - Multiple Positions ===\n');

PE6 = randn(50, 16);
layer6 = LearnedPositionalEmbedding(PE6);

seq_len = 4;
x6 = randn(16, seq_len);
positions = [0, 1, 2, 3];
y6 = layer6.evaluate(x6, positions);

pe_expected = PE6(positions + 1, :)';  % [EmbedDim x SeqLen]
expected6 = x6 + pe_expected;
assert(all(abs(y6(:) - expected6(:)) < 1e-10), 'Multi-position evaluation failed');

fprintf('Test 6 PASSED: Multiple position evaluation\n');

%% Test 7: Evaluate - Default Positions
fprintf('\n=== Test 7: Evaluate - Default Positions ===\n');

PE7 = randn(100, 8);
layer7 = LearnedPositionalEmbedding(PE7);

x7 = randn(8, 5);
y7 = layer7.evaluate(x7);  % Should use positions 0,1,2,3,4

expected_positions = 0:4;
pe_expected7 = PE7(expected_positions + 1, :)';
expected7 = x7 + pe_expected7;

assert(all(abs(y7(:) - expected7(:)) < 1e-10), 'Default position evaluation failed');

fprintf('Test 7 PASSED: Default positions\n');

%% Test 8: Reachability - Star Set
fprintf('\n=== Test 8: Reachability - Star Set ===\n');

PE8 = randn(20, 4);
layer8 = LearnedPositionalEmbedding(PE8);

lb8 = [-1; -1; -1; -1];
ub8 = [1; 1; 1; 1];
S_in = Star(lb8, ub8);

S_out = layer8.reach(S_in, 'approx-star', 0);

assert(~isempty(S_out), 'Reachability returned empty');
assert(isa(S_out, 'Star'), 'Output should be Star');

fprintf('Test 8 PASSED: Star reachability\n');

%% Test 9: Reachability - Soundness Check
fprintf('\n=== Test 9: Reachability - Soundness Check ===\n');

PE9 = randn(30, 4);
layer9 = LearnedPositionalEmbedding(PE9);

lb9 = [-0.5; -0.5; -0.5; -0.5];
ub9 = [0.5; 0.5; 0.5; 0.5];
S_in9 = Star(lb9, ub9);

position = 5;
S_out9 = layer9.reach(S_in9, 'approx-star', position);
[out_lb, out_ub] = S_out9.getRanges;

pe_5 = layer9.get_encoding(position);

% Sample and verify soundness
n_samples = 50;
all_contained = true;

for i = 1:n_samples
    x_sample = lb9 + (ub9 - lb9) .* rand(4, 1);
    y_sample = x_sample + pe_5;

    tol = 1e-6;
    if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
        all_contained = false;
    end
end

assert(all_contained, 'Soundness violation');
fprintf('Test 9 PASSED: Soundness check (%d samples)\n', n_samples);

%% Test 10: Position Similarity
fprintf('\n=== Test 10: Position Similarity ===\n');

PE10 = randn(20, 8);
% Normalize rows for predictable similarity
PE10 = PE10 ./ vecnorm(PE10, 2, 2);
layer10 = LearnedPositionalEmbedding(PE10);

% Self-similarity should be 1
sim_self = layer10.position_similarity(0, 0);
assert(abs(sim_self - 1) < 1e-6, 'Self-similarity should be 1');

% Other similarity should be in [-1, 1]
sim_other = layer10.position_similarity(0, 10);
assert(sim_other >= -1 - 1e-6 && sim_other <= 1 + 1e-6, 'Similarity out of range');

fprintf('Test 10 PASSED: Position similarity\n');

%% Test 11: Set Embeddings Method
fprintf('\n=== Test 11: Set Embeddings Method ===\n');

layer11 = LearnedPositionalEmbedding(10, 4);

new_PE = randn(50, 8);
layer11 = layer11.set_embeddings(new_PE);

assert(layer11.MaxSeqLen == 50, 'MaxSeqLen not updated');
assert(layer11.EmbedDim == 8, 'EmbedDim not updated');
assert(all(all(layer11.PE == new_PE)), 'PE not updated');

fprintf('Test 11 PASSED: Set embeddings method\n');

%% Test 12: Initialize Methods
fprintf('\n=== Test 12: Initialize Methods ===\n');

layer12 = LearnedPositionalEmbedding();
layer12 = layer12.initialize_zeros(16, 32);

assert(layer12.EmbedDim == 16, 'EmbedDim not set');
assert(layer12.MaxSeqLen == 32, 'MaxSeqLen not set');
assert(all(layer12.PE(:) == 0), 'PE should be all zeros');

layer12b = LearnedPositionalEmbedding();
layer12b = layer12b.initialize_random(16, 32);
assert(any(layer12b.PE(:) ~= 0), 'Random PE should not be all zeros');

fprintf('Test 12 PASSED: Initialize methods\n');

%% Test 13: Get Info
fprintf('\n=== Test 13: Get Info ===\n');

layer13 = LearnedPositionalEmbedding(64, 128, 'test_layer');
info = layer13.get_info();

assert(strcmp(info.Name, 'test_layer'), 'Name in info wrong');
assert(info.EmbedDim == 64, 'EmbedDim in info wrong');
assert(info.MaxSeqLen == 128, 'MaxSeqLen in info wrong');
assert(info.NumParameters == 64 * 128, 'NumParameters wrong');

fprintf('Test 13 PASSED: Get info\n');

%% Test 14: From Pretrained
fprintf('\n=== Test 14: From Pretrained ===\n');

pretrained_weights = randn(256, 768);
layer14 = LearnedPositionalEmbedding.from_pretrained(pretrained_weights, 'bert_pe');

assert(strcmp(layer14.Name, 'bert_pe'), 'Name not set from pretrained');
assert(layer14.MaxSeqLen == 256, 'MaxSeqLen from pretrained wrong');
assert(layer14.EmbedDim == 768, 'EmbedDim from pretrained wrong');

fprintf('Test 14 PASSED: From pretrained\n');

%% Test 15: Position Out of Range
fprintf('\n=== Test 15: Position Out of Range ===\n');

layer15 = LearnedPositionalEmbedding(10, 4);

try
    layer15.get_encoding(-1);  % Invalid: negative
    error('Should have thrown error for negative position');
catch ME
    assert(contains(ME.message, 'out of range'), 'Wrong error message for negative');
end

try
    layer15.get_encoding(10);  % Invalid: >= MaxSeqLen
    error('Should have thrown error for position >= MaxSeqLen');
catch ME
    assert(contains(ME.message, 'out of range'), 'Wrong error message for overflow');
end

fprintf('Test 15 PASSED: Position validation\n');

%% Summary
fprintf('\n=== All LearnedPositionalEmbedding Tests PASSED ===\n');
fprintf('Total: 15 tests\n');
