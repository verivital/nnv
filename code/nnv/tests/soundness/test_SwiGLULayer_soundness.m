%% Comprehensive Soundness Tests for SwiGLULayer
% Verifies that the SwiGLU reachability analysis is sound (overapproximates true outputs)
%
% SwiGLU(x) = SiLU(x @ W_gate + b_gate) .* (x @ W_up + b_up)
%
% Tests multiple input configurations:
% - Various dimensions
% - Various perturbation sizes
% - Different weight configurations
% - Edge cases
%
% Author: NNV Team
% Date: December 2025

%% Test 1: Small Network - Small Perturbation
fprintf('=== Test 1: Small Network - Small Perturbation ===\n');

rng(42);
input_dim = 2;
hidden_dim = 3;
W_gate = randn(input_dim, hidden_dim) * 0.5;
W_up = randn(input_dim, hidden_dim) * 0.5;
b_gate = zeros(hidden_dim, 1);
b_up = zeros(hidden_dim, 1);

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = [0.5; -0.3];
eps = 0.01;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: small network, small perturbation');
fprintf('Test 1 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 2: Small Network - Medium Perturbation
fprintf('\n=== Test 2: Small Network - Medium Perturbation ===\n');

rng(43);
input_dim = 2;
hidden_dim = 3;
W_gate = randn(input_dim, hidden_dim) * 0.5;
W_up = randn(input_dim, hidden_dim) * 0.5;
b_gate = zeros(hidden_dim, 1);
b_up = zeros(hidden_dim, 1);

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = [1.0; -1.0];
eps = 0.05;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: small network, medium perturbation');
fprintf('Test 2 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 3: Small Network - Large Perturbation
fprintf('\n=== Test 3: Small Network - Large Perturbation ===\n');

rng(44);
input_dim = 2;
hidden_dim = 3;
W_gate = randn(input_dim, hidden_dim) * 0.3;
W_up = randn(input_dim, hidden_dim) * 0.3;
b_gate = zeros(hidden_dim, 1);
b_up = zeros(hidden_dim, 1);

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = [0; 0];
eps = 0.2;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: small network, large perturbation');
fprintf('Test 3 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 4: Medium Network (4 -> 6)
fprintf('\n=== Test 4: Medium Network (4 -> 6) ===\n');

rng(45);
input_dim = 4;
hidden_dim = 6;
W_gate = randn(input_dim, hidden_dim) * 0.4;
W_up = randn(input_dim, hidden_dim) * 0.4;
b_gate = randn(hidden_dim, 1) * 0.1;
b_up = randn(hidden_dim, 1) * 0.1;

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = randn(input_dim, 1) * 0.5;
eps = 0.03;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 300);

assert(all_contained, 'Soundness violation: medium network');
fprintf('Test 4 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 5: With Non-zero Biases
fprintf('\n=== Test 5: With Non-zero Biases ===\n');

rng(46);
input_dim = 3;
hidden_dim = 4;
W_gate = randn(input_dim, hidden_dim) * 0.5;
W_up = randn(input_dim, hidden_dim) * 0.5;
b_gate = randn(hidden_dim, 1) * 0.5;
b_up = randn(hidden_dim, 1) * 0.5;

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = randn(input_dim, 1);
eps = 0.05;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: with biases');
fprintf('Test 5 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 6: Zero Center
fprintf('\n=== Test 6: Zero Center ===\n');

rng(47);
input_dim = 3;
hidden_dim = 4;
W_gate = randn(input_dim, hidden_dim) * 0.5;
W_up = randn(input_dim, hidden_dim) * 0.5;
b_gate = zeros(hidden_dim, 1);
b_up = zeros(hidden_dim, 1);

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = zeros(input_dim, 1);
eps = 0.1;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: zero center');
fprintf('Test 6 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 7: Large Input Values
fprintf('\n=== Test 7: Large Input Values ===\n');

rng(48);
input_dim = 2;
hidden_dim = 3;
W_gate = randn(input_dim, hidden_dim) * 0.3;
W_up = randn(input_dim, hidden_dim) * 0.3;
b_gate = zeros(hidden_dim, 1);
b_up = zeros(hidden_dim, 1);

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = [3; -3];  % Large values
eps = 0.1;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: large input values');
fprintf('Test 7 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 8: Identity-like Weights
fprintf('\n=== Test 8: Identity-like Weights ===\n');

input_dim = 3;
hidden_dim = 3;
W_gate = eye(input_dim, hidden_dim);
W_up = eye(input_dim, hidden_dim);
b_gate = zeros(hidden_dim, 1);
b_up = zeros(hidden_dim, 1);

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = [1; 0; -1];
eps = 0.1;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: identity weights');
fprintf('Test 8 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 9: Sparse Weights
fprintf('\n=== Test 9: Sparse Weights ===\n');

rng(49);
input_dim = 4;
hidden_dim = 4;
W_gate = randn(input_dim, hidden_dim) * 0.5;
W_up = randn(input_dim, hidden_dim) * 0.5;
% Make sparse
W_gate(abs(W_gate) < 0.3) = 0;
W_up(abs(W_up) < 0.3) = 0;
b_gate = zeros(hidden_dim, 1);
b_up = zeros(hidden_dim, 1);

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = randn(input_dim, 1) * 0.5;
eps = 0.05;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: sparse weights');
fprintf('Test 9 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 10: Stress Test (6 -> 8)
fprintf('\n=== Test 10: Stress Test (6 -> 8) ===\n');

rng(50);
input_dim = 6;
hidden_dim = 8;
W_gate = randn(input_dim, hidden_dim) * 0.3;
W_up = randn(input_dim, hidden_dim) * 0.3;
b_gate = randn(hidden_dim, 1) * 0.1;
b_up = randn(hidden_dim, 1) * 0.1;

layer = SwiGLULayer(W_gate, b_gate, W_up, b_up, 'swiglu_soundness');

center = randn(input_dim, 1) * 0.3;
eps = 0.02;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 500);

assert(all_contained, 'Soundness violation: stress test');
fprintf('Test 10 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Summary
fprintf('\n=== All SwiGLULayer Soundness Tests PASSED ===\n');
fprintf('Total: 10 tests\n');

%% Helper Function
function [all_contained, n_samples] = run_soundness_check(layer, center, eps, n_samples)
    input_dim = length(center);
    lb = center - eps;
    ub = center + eps;
    S_in = Star(lb, ub);

    S_out = layer.reach(S_in, 'approx-star');
    [out_lb, out_ub] = S_out.getRanges();
    out_lb = out_lb(:);
    out_ub = out_ub(:);

    all_contained = true;
    tol = 1e-5;

    for i = 1:n_samples
        x_sample = lb + (ub - lb) .* rand(input_dim, 1);
        y_sample = layer.evaluate(x_sample);
        if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
            all_contained = false;
            fprintf('  Sample %d OUTSIDE bounds\n', i);
            break;
        end
    end
end
