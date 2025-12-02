%% Test Suite for EmbeddingLayer
% Tests for word/token embedding lookup layer
%
% Author: NNV Team
% Date: November 2025

%% Test 1: Basic Construction with Matrix
fprintf('=== Test 1: Basic Construction with Matrix ===\n');

E = randn(100, 32);  % 100 vocab, 32 embed dim
layer = EmbeddingLayer(E);

assert(layer.VocabSize == 100, 'VocabSize not set');
assert(layer.EmbedDim == 32, 'EmbedDim not set');
assert(all(all(layer.E == E)), 'Embedding matrix not set');

fprintf('Test 1 PASSED: Basic construction\n');

%% Test 2: Construction with Dimensions
fprintf('\n=== Test 2: Construction with Dimensions ===\n');

layer2 = EmbeddingLayer(500, 64);

assert(layer2.VocabSize == 500, 'VocabSize not set');
assert(layer2.EmbedDim == 64, 'EmbedDim not set');
assert(size(layer2.E, 1) == 500, 'E rows wrong');
assert(size(layer2.E, 2) == 64, 'E cols wrong');

fprintf('Test 2 PASSED: Construction with dimensions\n');

%% Test 3: Single Token Lookup
fprintf('\n=== Test 3: Single Token Lookup ===\n');

E3 = randn(10, 4);
layer3 = EmbeddingLayer(E3);

token_id = 5;
y = layer3.evaluate(token_id);

expected = E3(token_id, :)';
assert(all(abs(y - expected) < 1e-10), 'Single token lookup failed');

fprintf('Test 3 PASSED: Single token lookup\n');

%% Test 4: Multiple Token Lookup
fprintf('\n=== Test 4: Multiple Token Lookup ===\n');

E4 = randn(20, 8);
layer4 = EmbeddingLayer(E4);

token_ids = [1, 5, 10, 15];
y4 = layer4.evaluate(token_ids);

expected4 = E4(token_ids, :)';  % [EmbedDim x SeqLen]
assert(all(abs(y4(:) - expected4(:)) < 1e-10), 'Multiple token lookup failed');

fprintf('Test 4 PASSED: Multiple token lookup\n');

%% Test 5: Output Shape
fprintf('\n=== Test 5: Output Shape ===\n');

layer5 = EmbeddingLayer(100, 32);

y5 = layer5.evaluate([1, 2, 3, 4, 5]);
assert(size(y5, 1) == 32, 'Output should have EmbedDim rows');
assert(size(y5, 2) == 5, 'Output should have SeqLen cols');

fprintf('Test 5 PASSED: Output shape\n');

%% Test 6: Reachability with Epsilon
fprintf('\n=== Test 6: Reachability with Epsilon ===\n');

E6 = randn(50, 16);
layer6 = EmbeddingLayer(E6);

token_ids6 = [10, 20, 30];
epsilon = 0.1;

R6 = layer6.reach_with_epsilon(token_ids6, epsilon);

assert(~isempty(R6), 'Reachability returned empty');
assert(isa(R6, 'Star'), 'Output should be Star');

% Check bounds
[lb, ub] = R6.getRanges;
base_emb = layer6.evaluate(token_ids6);
base_emb = base_emb(:);

assert(all(lb <= base_emb + 1e-6), 'Lower bound too high');
assert(all(ub >= base_emb - 1e-6), 'Upper bound too low');

fprintf('Test 6 PASSED: Reachability with epsilon\n');

%% Test 7: Reachability with Epsilon - Soundness
fprintf('\n=== Test 7: Reachability with Epsilon - Soundness ===\n');

E7 = randn(30, 8);
layer7 = EmbeddingLayer(E7);

token_ids7 = [5, 15];
epsilon7 = 0.05;

R7 = layer7.reach_with_epsilon(token_ids7, epsilon7);
[out_lb, out_ub] = R7.getRanges;

base_emb7 = layer7.evaluate(token_ids7);
base_emb7 = base_emb7(:);

n_samples = 50;
all_contained = true;

for i = 1:n_samples
    % Random perturbation within epsilon
    perturb = (rand(size(base_emb7)) * 2 - 1) * epsilon7;
    y_sample = base_emb7 + perturb;

    tol = 1e-6;
    if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
        all_contained = false;
    end
end

assert(all_contained, 'Soundness violation');
fprintf('Test 7 PASSED: Epsilon soundness (%d samples)\n', n_samples);

%% Test 8: Reachability with Synonyms
fprintf('\n=== Test 8: Reachability with Synonyms ===\n');

E8 = randn(20, 4);
layer8 = EmbeddingLayer(E8);

token_id8 = 5;
synonym_ids8 = [5, 6, 7];  % Token 5 and its synonyms

R8 = layer8.reach_with_synonyms(token_id8, synonym_ids8);

assert(~isempty(R8), 'Synonym reachability returned empty');
assert(isa(R8, 'Star'), 'Output should be Star');

% Check that all synonym embeddings are within bounds
[lb8, ub8] = R8.getRanges;
for id = synonym_ids8
    emb = layer8.get_embedding(id);
    tol = 1e-6;
    assert(all(emb >= lb8 - tol) && all(emb <= ub8 + tol), ...
        'Synonym embedding not contained');
end

fprintf('Test 8 PASSED: Reachability with synonyms\n');

%% Test 9: Cosine Similarity
fprintf('\n=== Test 9: Cosine Similarity ===\n');

E9 = randn(10, 8);
E9 = E9 ./ vecnorm(E9, 2, 2);  % Normalize rows
layer9 = EmbeddingLayer(E9);

% Self-similarity should be 1
sim_self = layer9.cosine_similarity(1, 1);
assert(abs(sim_self - 1) < 1e-6, 'Self-similarity should be 1');

% Similarity should be in [-1, 1]
sim_other = layer9.cosine_similarity(1, 5);
assert(sim_other >= -1 && sim_other <= 1, 'Similarity out of range');

fprintf('Test 9 PASSED: Cosine similarity\n');

%% Test 10: Find Nearest Neighbors
fprintf('\n=== Test 10: Find Nearest Neighbors ===\n');

E10 = randn(50, 16);
layer10 = EmbeddingLayer(E10);

[ids, sims] = layer10.find_nearest(1, 5);

assert(length(ids) == 5, 'Should return 5 neighbors');
assert(all(diff(sims) <= 0), 'Similarities should be sorted descending');
assert(~ismember(1, ids), 'Should not include query token');

fprintf('Test 10 PASSED: Find nearest neighbors\n');

%% Test 11: Named Constructor
fprintf('\n=== Test 11: Named Constructor ===\n');

E11 = randn(20, 8);
layer11 = EmbeddingLayer(E11, 'word_embedding');

assert(strcmp(layer11.Name, 'word_embedding'), 'Name not set');

fprintf('Test 11 PASSED: Named constructor\n');

%% Test 12: Set Embeddings Method
fprintf('\n=== Test 12: Set Embeddings Method ===\n');

layer12 = EmbeddingLayer(10, 4);

new_E = randn(25, 8);
layer12 = layer12.set_embeddings(new_E);

assert(layer12.VocabSize == 25, 'VocabSize not updated');
assert(layer12.EmbedDim == 8, 'EmbedDim not updated');

fprintf('Test 12 PASSED: Set embeddings method\n');

%% Test 13: Token Index Validation
fprintf('\n=== Test 13: Token Index Validation ===\n');

layer13 = EmbeddingLayer(10, 4);

try
    layer13.evaluate(0);  % Invalid: 0 index
    error('Should have thrown error for index 0');
catch ME
    assert(contains(ME.message, 'out of range'), 'Wrong error message');
end

try
    layer13.evaluate(11);  % Invalid: > VocabSize
    error('Should have thrown error for index > VocabSize');
catch ME
    assert(contains(ME.message, 'out of range'), 'Wrong error message');
end

fprintf('Test 13 PASSED: Token index validation\n');

%% Test 14: evaluate_single Method
fprintf('\n=== Test 14: evaluate_single Method ===\n');

E14 = randn(15, 6);
layer14 = EmbeddingLayer(E14);

y14 = layer14.evaluate_single(7);
expected14 = E14(7, :)';

assert(all(abs(y14 - expected14) < 1e-10), 'evaluate_single failed');
assert(isvector(y14), 'Output should be vector');
assert(length(y14) == 6, 'Output length wrong');

fprintf('Test 14 PASSED: evaluate_single method\n');

%% Summary
fprintf('\n=== All EmbeddingLayer Tests PASSED ===\n');
fprintf('Total: 14 tests\n');
