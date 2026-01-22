%% Test Suite for RotaryPositionalEmbedding (RoPE)
% Tests for Rotary Position Embedding layer
%
% RoPE applies rotation matrices to pairs of dimensions based on position.
% Key property: For fixed positions, RoPE is a LINEAR transformation.
%
% Author: NNV Team
% Date: November 2025

%% Test 1: Basic Construction
fprintf('=== Test 1: Basic Construction ===\n');

dim = 64;
layer = RotaryPositionalEmbedding(dim);

assert(layer.Dim == dim, 'Dim not set correctly');
assert(layer.Base == 10000, 'Default base should be 10000');
assert(layer.MaxSeqLen == 512, 'Default MaxSeqLen should be 512');
assert(length(layer.InvFreq) == dim/2, 'InvFreq length wrong');

fprintf('Test 1 PASSED: Basic construction\n');

%% Test 2: Construction with Custom Parameters
fprintf('\n=== Test 2: Construction with Custom Parameters ===\n');

layer2 = RotaryPositionalEmbedding(32, 8000, 1024, 'my_rope');

assert(layer2.Dim == 32, 'Dim not set');
assert(layer2.Base == 8000, 'Base not set');
assert(layer2.MaxSeqLen == 1024, 'MaxSeqLen not set');
assert(strcmp(layer2.Name, 'my_rope'), 'Name not set');

fprintf('Test 2 PASSED: Custom parameters\n');

%% Test 3: Odd Dimension Error
fprintf('\n=== Test 3: Odd Dimension Error ===\n');

try
    layer_bad = RotaryPositionalEmbedding(33);  % Odd dimension
    error('Should have thrown error for odd dimension');
catch ME
    assert(contains(ME.message, 'even'), 'Wrong error message');
end

fprintf('Test 3 PASSED: Odd dimension rejected\n');

%% Test 4: Inverse Frequency Computation
fprintf('\n=== Test 4: Inverse Frequency Computation ===\n');

layer4 = RotaryPositionalEmbedding(8, 10000);

% inv_freq[i] = 1 / 10000^(2i/8) for i = 0, 1, 2, 3
expected_inv_freq = [1; 1/10000^(2/8); 1/10000^(4/8); 1/10000^(6/8)];
assert(all(abs(layer4.InvFreq - expected_inv_freq) < 1e-10), 'InvFreq computation wrong');

fprintf('Test 4 PASSED: Inverse frequency computation\n');

%% Test 5: Get Cos/Sin Values
fprintf('\n=== Test 5: Get Cos/Sin Values ===\n');

layer5 = RotaryPositionalEmbedding(4);
[cos_vals, sin_vals] = layer5.get_cos_sin([0, 1, 2]);

% At position 0: cos(0) = 1, sin(0) = 0 for all frequencies
assert(all(abs(cos_vals(1, :) - 1) < 1e-10), 'cos(0) should be 1');
assert(all(abs(sin_vals(1, :)) < 1e-10), 'sin(0) should be 0');

% Shape check
assert(size(cos_vals, 1) == 3, 'Wrong number of positions');
assert(size(cos_vals, 2) == 2, 'Wrong number of frequencies (Dim/2)');

fprintf('Test 5 PASSED: Cos/Sin computation\n');

%% Test 6: Rotation Matrix Properties
fprintf('\n=== Test 6: Rotation Matrix Properties ===\n');

layer6 = RotaryPositionalEmbedding(8);
R = layer6.get_rotation_matrix(5);

% Rotation matrix should be orthogonal: R * R' = I
ortho_err = norm(R * R' - eye(8), 'fro');
assert(ortho_err < 1e-10, 'Rotation matrix not orthogonal');

% Determinant should be 1 (proper rotation)
assert(abs(det(R) - 1) < 1e-10, 'Determinant should be 1');

% Should be square with correct dimension
assert(all(size(R) == [8, 8]), 'Rotation matrix size wrong');

fprintf('Test 6 PASSED: Rotation matrix is orthogonal with det=1\n');

%% Test 7: Evaluation - Position 0
fprintf('\n=== Test 7: Evaluation - Position 0 ===\n');

layer7 = RotaryPositionalEmbedding(4);
x7 = [1; 2; 3; 4];
y7 = layer7.evaluate(x7, 0);

% At position 0: rotation angle is 0, so output = input
% cos(0) = 1, sin(0) = 0 -> identity rotation
assert(all(abs(y7 - x7) < 1e-10), 'Position 0 should be identity');

fprintf('Test 7 PASSED: Position 0 is identity\n');

%% Test 8: Evaluation - Non-zero Position
fprintf('\n=== Test 8: Evaluation - Non-zero Position ===\n');

layer8 = RotaryPositionalEmbedding(4);
x8 = [1; 0; 1; 0];  % Simple test vector
pos = 1;

y8 = layer8.evaluate(x8, pos);

% Manual computation
[cos_p, sin_p] = layer8.get_cos_sin(pos);
% x_even = [1; 1], x_odd = [0; 0]
% y_even = 1*cos - 0*sin = cos
% y_odd = 1*sin + 0*cos = sin
expected_even = cos_p(:);
expected_odd = sin_p(:);
expected8 = zeros(4, 1);
expected8(1:2:end) = expected_even;
expected8(2:2:end) = expected_odd;

assert(all(abs(y8 - expected8) < 1e-10), 'Non-zero position evaluation wrong');

fprintf('Test 8 PASSED: Non-zero position evaluation\n');

%% Test 9: Evaluation - Multiple Positions (Sequence)
fprintf('\n=== Test 9: Evaluation - Multiple Positions ===\n');

layer9 = RotaryPositionalEmbedding(8);
seq_len = 4;
x9 = randn(8, seq_len);
positions = [0, 1, 2, 3];

y9 = layer9.evaluate(x9, positions);

% Check shape
assert(all(size(y9) == size(x9)), 'Output shape should match input');

% First position (0) should be identity
assert(all(abs(y9(:, 1) - x9(:, 1)) < 1e-10), 'Position 0 should be identity');

fprintf('Test 9 PASSED: Multiple positions\n');

%% Test 10: Rotation Matrix Matches Evaluate
fprintf('\n=== Test 10: Rotation Matrix Matches Evaluate ===\n');

layer10 = RotaryPositionalEmbedding(8);
x10 = randn(8, 1);
pos = 7;

% Method 1: evaluate()
y_eval = layer10.evaluate(x10, pos);

% Method 2: rotation matrix
R10 = layer10.get_rotation_matrix(pos);
y_matrix = R10 * x10;

assert(all(abs(y_eval - y_matrix) < 1e-10), 'evaluate() and matrix multiplication should match');

fprintf('Test 10 PASSED: Rotation matrix matches evaluate\n');

%% Test 11: Reachability - Star Set (Single Position)
fprintf('\n=== Test 11: Reachability - Star Set ===\n');

layer11 = RotaryPositionalEmbedding(4);

lb11 = [-1; -1; -1; -1];
ub11 = [1; 1; 1; 1];
S_in = Star(lb11, ub11);

S_out = layer11.reach(S_in, 'approx-star', 3);

assert(~isempty(S_out), 'Reachability returned empty');
assert(isa(S_out, 'Star'), 'Output should be Star');
assert(S_out.dim == 4, 'Output dimension should match input');

fprintf('Test 11 PASSED: Star reachability\n');

%% Test 12: Reachability - Soundness Check
fprintf('\n=== Test 12: Reachability - Soundness Check ===\n');

layer12 = RotaryPositionalEmbedding(4);

lb12 = [-0.5; -0.5; -0.5; -0.5];
ub12 = [0.5; 0.5; 0.5; 0.5];
S_in12 = Star(lb12, ub12);

pos = 5;
S_out12 = layer12.reach(S_in12, 'approx-star', pos);
[out_lb, out_ub] = S_out12.getRanges;

% Sample and verify soundness
n_samples = 100;
all_contained = true;

for i = 1:n_samples
    x_sample = lb12 + (ub12 - lb12) .* rand(4, 1);
    y_sample = layer12.evaluate(x_sample, pos);

    tol = 1e-6;
    if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
        all_contained = false;
    end
end

assert(all_contained, 'Soundness violation');
fprintf('Test 12 PASSED: Soundness check (%d samples)\n', n_samples);

%% Test 13: Reachability - Exactness (RoPE is Linear!)
fprintf('\n=== Test 13: Reachability - Exactness ===\n');

% Since RoPE is linear, reachability should be EXACT (not over-approximation)
layer13 = RotaryPositionalEmbedding(4);

lb13 = [0; 0; 0; 0];
ub13 = [1; 1; 1; 1];
S_in13 = Star(lb13, ub13);

pos = 2;
S_out13 = layer13.reach(S_in13, 'approx-star', pos);
[out_lb, out_ub] = S_out13.getRanges;

% The output bounds should be tight (since rotation is linear)
% Check by evaluating at corners
corners = [0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1;
           0 1 0 0; 0 1 0 1; 0 1 1 0; 0 1 1 1;
           1 0 0 0; 1 0 0 1; 1 0 1 0; 1 0 1 1;
           1 1 0 0; 1 1 0 1; 1 1 1 0; 1 1 1 1]';

y_corners = zeros(4, 16);
for i = 1:16
    y_corners(:, i) = layer13.evaluate(corners(:, i), pos);
end

actual_lb = min(y_corners, [], 2);
actual_ub = max(y_corners, [], 2);

% For linear transformation, bounds should be exact (within numerical tolerance)
assert(all(abs(out_lb - actual_lb) < 1e-6), 'Lower bounds not exact');
assert(all(abs(out_ub - actual_ub) < 1e-6), 'Upper bounds not exact');

fprintf('Test 13 PASSED: Reachability is exact for linear RoPE\n');

%% Test 14: Multiple Positions Reachability
fprintf('\n=== Test 14: Multiple Positions Reachability ===\n');

layer14 = RotaryPositionalEmbedding(4);
seq_len = 3;

% Input Star for flattened sequence [Dim * SeqLen]
lb14 = -ones(4 * seq_len, 1);
ub14 = ones(4 * seq_len, 1);
S_in14 = Star(lb14, ub14);

positions = [0, 1, 2];
S_out14 = layer14.reach(S_in14, 'approx-star', positions);

assert(S_out14.dim == 4 * seq_len, 'Output dimension wrong for sequence');

fprintf('Test 14 PASSED: Multiple positions reachability\n');

%% Test 15: Get Info
fprintf('\n=== Test 15: Get Info ===\n');

layer15 = RotaryPositionalEmbedding(64, 10000, 512, 'test_rope');
info = layer15.get_info();

assert(strcmp(info.Name, 'test_rope'), 'Name wrong');
assert(info.Dim == 64, 'Dim wrong');
assert(info.Base == 10000, 'Base wrong');
assert(info.NumParameters == 0, 'RoPE should have 0 learnable parameters');

fprintf('Test 15 PASSED: Get info\n');

%% Test 16: Position Encoding Varies with Position
fprintf('\n=== Test 16: Position Encoding Varies ===\n');

layer16 = RotaryPositionalEmbedding(8);
x16 = randn(8, 1);

y_pos0 = layer16.evaluate(x16, 0);
y_pos1 = layer16.evaluate(x16, 1);
y_pos10 = layer16.evaluate(x16, 10);

% Different positions should give different outputs
assert(norm(y_pos0 - y_pos1) > 0.01, 'Different positions should give different outputs');
assert(norm(y_pos1 - y_pos10) > 0.01, 'Different positions should give different outputs');

fprintf('Test 16 PASSED: Position encoding varies\n');

%% Test 17: Large Position Values
fprintf('\n=== Test 17: Large Position Values ===\n');

layer17 = RotaryPositionalEmbedding(16, 10000, 4096);
x17 = randn(16, 1);

% Test large position values
y_large = layer17.evaluate(x17, 4000);

assert(~any(isnan(y_large)), 'NaN at large position');
assert(~any(isinf(y_large)), 'Inf at large position');

fprintf('Test 17 PASSED: Large position values\n');

%% Test 18: Integration with Attention-like Usage
fprintf('\n=== Test 18: Attention Integration ===\n');

% Simulate how RoPE would be used in attention
dim = 16;
seq_len = 4;
layer18 = RotaryPositionalEmbedding(dim);

% Create Q and K matrices
Q = randn(dim, seq_len);
K = randn(dim, seq_len);
positions = 0:(seq_len-1);

% Apply RoPE to Q and K
Q_rot = layer18.evaluate(Q, positions);
K_rot = layer18.evaluate(K, positions);

% Compute attention scores (Q_rot' * K_rot should be position-aware)
attn_scores = Q_rot' * K_rot;

assert(all(size(attn_scores) == [seq_len, seq_len]), 'Attention scores shape wrong');
assert(~any(isnan(attn_scores(:))), 'NaN in attention scores');

fprintf('Test 18 PASSED: Attention integration\n');

%% Summary
fprintf('\n=== All RotaryPositionalEmbedding Tests PASSED ===\n');
fprintf('Total: 18 tests\n');
