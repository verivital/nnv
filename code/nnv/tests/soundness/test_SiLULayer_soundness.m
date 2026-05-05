%% Comprehensive Soundness Tests for SiLULayer
% Verifies that the SiLU reachability analysis is sound (overapproximates true outputs)
%
% Tests multiple input configurations:
% - Various dimensions (2D, 4D, 8D, 16D)
% - Various perturbation sizes (small, medium, large)
% - Star and ImageStar inputs
% - Edge cases (near zero, large values, negative values)
%
% Author: NNV Team
% Date: December 2025

%% Test 1: 2D Star - Small Perturbation
fprintf('=== Test 1: 2D Star - Small Perturbation ===\n');

layer = SiLULayer('silu_soundness');
input_dim = 2;
center = [0.5; -0.3];
eps = 0.01;

[all_contained, n_samples, out_lb, out_ub] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: 2D small perturbation');
fprintf('Test 1 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 2: 2D Star - Medium Perturbation
fprintf('\n=== Test 2: 2D Star - Medium Perturbation ===\n');

layer = SiLULayer('silu_soundness');
center = [1.0; -1.0];
eps = 0.1;

[all_contained, n_samples, ~, ~] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: 2D medium perturbation');
fprintf('Test 2 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 3: 2D Star - Large Perturbation
fprintf('\n=== Test 3: 2D Star - Large Perturbation ===\n');

layer = SiLULayer('silu_soundness');
center = [0; 0];
eps = 1.0;

[all_contained, n_samples, ~, ~] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: 2D large perturbation');
fprintf('Test 3 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 4: 4D Star - Random Center
fprintf('\n=== Test 4: 4D Star - Random Center ===\n');

layer = SiLULayer('silu_soundness');
rng(42);
center = randn(4, 1);
eps = 0.1;

[all_contained, n_samples, ~, ~] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: 4D random center');
fprintf('Test 4 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 5: 8D Star - High Dimension
fprintf('\n=== Test 5: 8D Star - High Dimension ===\n');

layer = SiLULayer('silu_soundness');
rng(43);
center = randn(8, 1);
eps = 0.05;

[all_contained, n_samples, ~, ~] = run_soundness_check(layer, center, eps, 300);

assert(all_contained, 'Soundness violation: 8D high dimension');
fprintf('Test 5 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 6: Near SiLU Minimum
fprintf('\n=== Test 6: Near SiLU Minimum ===\n');

layer = SiLULayer('silu_soundness');
[x_min, ~] = SiLU.get_minimum();
center = [x_min; x_min];  % Both dims at minimum
eps = 0.2;

[all_contained, n_samples, ~, ~] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: near minimum');
fprintf('Test 6 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 7: Large Positive Values
fprintf('\n=== Test 7: Large Positive Values ===\n');

layer = SiLULayer('silu_soundness');
center = [5; 10];  % Large positive
eps = 0.5;

[all_contained, n_samples, ~, ~] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: large positive');
fprintf('Test 7 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 8: Large Negative Values
fprintf('\n=== Test 8: Large Negative Values ===\n');

layer = SiLULayer('silu_soundness');
center = [-5; -10];  % Large negative
eps = 0.5;

[all_contained, n_samples, ~, ~] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: large negative');
fprintf('Test 8 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 9: Mixed Signs
fprintf('\n=== Test 9: Mixed Signs ===\n');

layer = SiLULayer('silu_soundness');
center = [2; -2; 0.5; -0.5];
eps = 0.3;

[all_contained, n_samples, ~, ~] = run_soundness_check(layer, center, eps, 200);

assert(all_contained, 'Soundness violation: mixed signs');
fprintf('Test 9 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 10: ImageStar 2x2
fprintf('\n=== Test 10: ImageStar 2x2 ===\n');

layer = SiLULayer('silu_soundness');

center_img = randn(2, 2, 1);
eps = 0.1;
lb_img = center_img - eps;
ub_img = center_img + eps;
IS_in = ImageStar(lb_img, ub_img);

IS_out = layer.reach_star_single_input(IS_in, 'approx-star', 0, [], 'linprog');
[out_lb, out_ub] = IS_out.getRanges();

n_samples = 200;
all_contained = true;
tol = 1e-5;

for i = 1:n_samples
    x_sample = lb_img + (ub_img - lb_img) .* rand(size(lb_img));
    y_sample = layer.evaluate(x_sample);
    if any(y_sample(:) < out_lb(:) - tol) || any(y_sample(:) > out_ub(:) + tol)
        all_contained = false;
        break;
    end
end

assert(all_contained, 'Soundness violation: ImageStar 2x2');
fprintf('Test 10 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 11: ImageStar 4x4x3 (RGB-like)
fprintf('\n=== Test 11: ImageStar 4x4x3 (RGB-like) ===\n');

layer = SiLULayer('silu_soundness');

center_img = randn(4, 4, 3) * 0.5;
eps = 0.05;
lb_img = center_img - eps;
ub_img = center_img + eps;
IS_in = ImageStar(lb_img, ub_img);

IS_out = layer.reach_star_single_input(IS_in, 'approx-star', 0, [], 'linprog');
[out_lb, out_ub] = IS_out.getRanges();

n_samples = 100;
all_contained = true;
tol = 1e-5;

for i = 1:n_samples
    x_sample = lb_img + (ub_img - lb_img) .* rand(size(lb_img));
    y_sample = layer.evaluate(x_sample);
    if any(y_sample(:) < out_lb(:) - tol) || any(y_sample(:) > out_ub(:) + tol)
        all_contained = false;
        break;
    end
end

assert(all_contained, 'Soundness violation: ImageStar 4x4x3');
fprintf('Test 11 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 12: 16D Star - Stress Test
fprintf('\n=== Test 12: 16D Star - Stress Test ===\n');

layer = SiLULayer('silu_soundness');
rng(44);
center = randn(16, 1);
eps = 0.02;

[all_contained, n_samples, ~, ~] = run_soundness_check(layer, center, eps, 500);

assert(all_contained, 'Soundness violation: 16D stress test');
fprintf('Test 12 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Summary
fprintf('\n=== All SiLULayer Soundness Tests PASSED ===\n');
fprintf('Total: 12 tests\n');

%% Helper Function
function [all_contained, n_samples, out_lb, out_ub] = run_soundness_check(layer, center, eps, n_samples)
    input_dim = length(center);
    lb = center - eps;
    ub = center + eps;
    S_in = Star(lb, ub);

    S_out = layer.reach_star_single_input(S_in, 'approx-star', 0, [], 'linprog');
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
