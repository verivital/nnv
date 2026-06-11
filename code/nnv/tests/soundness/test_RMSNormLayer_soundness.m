%% Comprehensive Soundness Tests for RMSNormLayer
% Verifies that the RMSNorm reachability analysis is sound (overapproximates true outputs)
%
% RMSNorm(x) = x / sqrt(mean(x^2) + eps) * gamma
%
% Tests multiple input configurations:
% - Various dimensions
% - Various perturbation sizes
% - Different gamma values
% - Edge cases (near zero, scale invariance)
%
% Author: NNV Team
% Date: December 2025

%% Test 1: 2D - Default Gamma - Small Perturbation
fprintf('=== Test 1: 2D Default Gamma - Small Perturbation ===\n');

layer = RMSNormLayer(2);  % Default gamma = ones
center = [0.5; 0.3];
eps = 0.01;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: 2D default gamma, small perturbation');
fprintf('Test 1 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 2: 2D - Default Gamma - Medium Perturbation
fprintf('\n=== Test 2: 2D Default Gamma - Medium Perturbation ===\n');

layer = RMSNormLayer(2);
center = [1.0; -1.0];
eps = 0.05;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: 2D default gamma, medium perturbation');
fprintf('Test 2 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 3: 4D - Custom Gamma
fprintf('\n=== Test 3: 4D Custom Gamma ===\n');

gamma = [1.5; 0.8; 1.2; 0.5];
layer = RMSNormLayer(gamma);
center = [0.5; -0.5; 0.3; -0.3];
eps = 0.05;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: 4D custom gamma');
fprintf('Test 3 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 4: 3D - Large Gamma Values
fprintf('\n=== Test 4: 3D Large Gamma Values ===\n');

gamma = [5; 10; 2];
layer = RMSNormLayer(gamma);
center = [0.2; 0.4; 0.6];
eps = 0.03;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: large gamma values');
fprintf('Test 4 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 5: 3D - Small Gamma Values
fprintf('\n=== Test 5: 3D Small Gamma Values ===\n');

gamma = [0.1; 0.2; 0.05];
layer = RMSNormLayer(gamma);
center = [1.0; 2.0; 0.5];
eps = 0.05;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: small gamma values');
fprintf('Test 5 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 6: Near-Zero Input (Edge Case)
fprintf('\n=== Test 6: Near-Zero Input ===\n');

layer = RMSNormLayer(3);
center = [0.001; 0.002; 0.001];  % Very small
eps = 0.0005;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: near-zero input');
fprintf('Test 6 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 7: Large Input Values
fprintf('\n=== Test 7: Large Input Values ===\n');

layer = RMSNormLayer(3);
center = [10; -5; 8];
eps = 0.5;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: large input values');
fprintf('Test 7 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 8: All Negative Values
fprintf('\n=== Test 8: All Negative Values ===\n');

layer = RMSNormLayer(4);
center = [-1; -2; -0.5; -1.5];
eps = 0.1;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: all negative values');
fprintf('Test 8 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 9: Mixed Signs
fprintf('\n=== Test 9: Mixed Signs ===\n');

gamma = [1.2; 0.8; 1.5; 0.6];
layer = RMSNormLayer(gamma);
center = [2; -2; 0.5; -0.5];
eps = 0.1;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: mixed signs');
fprintf('Test 9 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 10: Custom Epsilon
fprintf('\n=== Test 10: Custom Epsilon ===\n');

gamma = ones(3, 1);
custom_eps = 1e-8;
layer = RMSNormLayer(gamma, custom_eps);
center = [0.1; 0.2; 0.3];
eps = 0.02;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: custom epsilon');
fprintf('Test 10 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 11: 8D - High Dimension
fprintf('\n=== Test 11: 8D High Dimension ===\n');

rng(42);
gamma = 0.5 + rand(8, 1);
layer = RMSNormLayer(gamma);
center = randn(8, 1);
eps = 0.03;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 300);

assert(all_contained, 'Soundness violation: 8D high dimension');
fprintf('Test 11 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 12: Stress Test - 16D
fprintf('\n=== Test 12: Stress Test - 16D ===\n');

rng(43);
gamma = 0.5 + rand(16, 1);
layer = RMSNormLayer(gamma);
center = randn(16, 1) * 0.5;
eps = 0.02;

[all_contained, n_samples] = run_soundness_check(layer, center, eps, 500);

assert(all_contained, 'Soundness violation: 16D stress test');
fprintf('Test 12 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Summary
fprintf('\n=== All RMSNormLayer Soundness Tests PASSED ===\n');
fprintf('Total: 12 tests\n');

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
