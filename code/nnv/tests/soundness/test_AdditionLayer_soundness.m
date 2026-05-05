%% Comprehensive Soundness Tests for AdditionLayer
% Verifies that the AdditionLayer reachability analysis is sound
%
% AdditionLayer computes: output = input1 + input2 + ... + inputN
%
% Tests multiple configurations:
% - Two inputs
% - Three inputs
% - Various perturbation sizes
% - Star and ImageStar inputs
%
% Author: NNV Team
% Date: December 2025

%% Test 1: Two Inputs - Small Perturbation
fprintf('=== Test 1: Two Inputs - Small Perturbation ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});
input_dim = 4;

rng(42);
center1 = randn(input_dim, 1);
center2 = randn(input_dim, 1);
eps = 0.01;

[all_contained, n_samples] = run_soundness_check_addition(layer, center1, center2, eps, 200);

assert(all_contained, 'Soundness violation: two inputs, small perturbation');
fprintf('Test 1 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 2: Two Inputs - Medium Perturbation
fprintf('\n=== Test 2: Two Inputs - Medium Perturbation ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});
input_dim = 4;

rng(43);
center1 = randn(input_dim, 1);
center2 = randn(input_dim, 1);
eps = 0.1;

[all_contained, n_samples] = run_soundness_check_addition(layer, center1, center2, eps, 200);

assert(all_contained, 'Soundness violation: two inputs, medium perturbation');
fprintf('Test 2 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 3: Two Inputs - Large Perturbation
fprintf('\n=== Test 3: Two Inputs - Large Perturbation ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});
input_dim = 4;

rng(44);
center1 = randn(input_dim, 1);
center2 = randn(input_dim, 1);
eps = 0.5;

[all_contained, n_samples] = run_soundness_check_addition(layer, center1, center2, eps, 200);

assert(all_contained, 'Soundness violation: two inputs, large perturbation');
fprintf('Test 3 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 4: Two Inputs - Zero Centers
fprintf('\n=== Test 4: Two Inputs - Zero Centers ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});
input_dim = 3;

center1 = zeros(input_dim, 1);
center2 = zeros(input_dim, 1);
eps = 0.1;

[all_contained, n_samples] = run_soundness_check_addition(layer, center1, center2, eps, 200);

assert(all_contained, 'Soundness violation: zero centers');
fprintf('Test 4 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 5: Two Inputs - Opposite Sign Centers
fprintf('\n=== Test 5: Two Inputs - Opposite Sign Centers ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});
input_dim = 4;

center1 = [1; 2; 3; 4];
center2 = [-1; -2; -3; -4];  % Should nearly cancel
eps = 0.1;

[all_contained, n_samples] = run_soundness_check_addition(layer, center1, center2, eps, 200);

assert(all_contained, 'Soundness violation: opposite sign centers');
fprintf('Test 5 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 6: Two Inputs - Large Values
fprintf('\n=== Test 6: Two Inputs - Large Values ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});
input_dim = 3;

center1 = [10; -5; 8];
center2 = [-3; 7; -2];
eps = 0.5;

[all_contained, n_samples] = run_soundness_check_addition(layer, center1, center2, eps, 200);

assert(all_contained, 'Soundness violation: large values');
fprintf('Test 6 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 7: Two Inputs - 8D High Dimension
fprintf('\n=== Test 7: Two Inputs - 8D High Dimension ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});
input_dim = 8;

rng(45);
center1 = randn(input_dim, 1);
center2 = randn(input_dim, 1);
eps = 0.05;

[all_contained, n_samples] = run_soundness_check_addition(layer, center1, center2, eps, 300);

assert(all_contained, 'Soundness violation: 8D high dimension');
fprintf('Test 7 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 8: Three Inputs
fprintf('\n=== Test 8: Three Inputs ===\n');

layer = AdditionLayer('add_soundness', 3, 1, {'in1', 'in2', 'in3'}, {'out'});
input_dim = 4;

rng(46);
center1 = randn(input_dim, 1);
center2 = randn(input_dim, 1);
center3 = randn(input_dim, 1);
eps = 0.05;

[all_contained, n_samples] = run_soundness_check_addition_3(layer, center1, center2, center3, eps, 200);

assert(all_contained, 'Soundness violation: three inputs');
fprintf('Test 8 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 9: ImageStar Inputs - 2D
fprintf('\n=== Test 9: ImageStar Inputs - 2D ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});

rng(47);
center_img1 = randn(2, 2, 1);
center_img2 = randn(2, 2, 1);
eps = 0.05;

lb1 = center_img1 - eps;
ub1 = center_img1 + eps;
lb2 = center_img2 - eps;
ub2 = center_img2 + eps;

IS_in1 = ImageStar(lb1, ub1);
IS_in2 = ImageStar(lb2, ub2);

inputs = {IS_in1, IS_in2};
IS_out = layer.reach(inputs, 'approx-star');
[out_lb, out_ub] = IS_out.getRanges();

n_samples = 200;
all_contained = true;
tol = 1e-5;

for i = 1:n_samples
    x1_sample = lb1 + (ub1 - lb1) .* rand(size(lb1));
    x2_sample = lb2 + (ub2 - lb2) .* rand(size(lb2));
    y_sample = layer.evaluate({x1_sample, x2_sample});
    if any(y_sample(:) < out_lb(:) - tol) || any(y_sample(:) > out_ub(:) + tol)
        all_contained = false;
        break;
    end
end

assert(all_contained, 'Soundness violation: ImageStar inputs');
fprintf('Test 9 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 10: ImageStar 4x4x2
fprintf('\n=== Test 10: ImageStar 4x4x2 ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});

rng(48);
center_img1 = randn(4, 4, 2) * 0.5;
center_img2 = randn(4, 4, 2) * 0.5;
eps = 0.03;

lb1 = center_img1 - eps;
ub1 = center_img1 + eps;
lb2 = center_img2 - eps;
ub2 = center_img2 + eps;

IS_in1 = ImageStar(lb1, ub1);
IS_in2 = ImageStar(lb2, ub2);

inputs = {IS_in1, IS_in2};
IS_out = layer.reach(inputs, 'approx-star');
[out_lb, out_ub] = IS_out.getRanges();

n_samples = 100;
all_contained = true;
tol = 1e-5;

for i = 1:n_samples
    x1_sample = lb1 + (ub1 - lb1) .* rand(size(lb1));
    x2_sample = lb2 + (ub2 - lb2) .* rand(size(lb2));
    y_sample = layer.evaluate({x1_sample, x2_sample});
    if any(y_sample(:) < out_lb(:) - tol) || any(y_sample(:) > out_ub(:) + tol)
        all_contained = false;
        break;
    end
end

assert(all_contained, 'Soundness violation: ImageStar 4x4x2');
fprintf('Test 10 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Test 11: Stress Test - 16D
fprintf('\n=== Test 11: Stress Test - 16D ===\n');

layer = AdditionLayer('add_soundness', 2, 1, {'in1', 'in2'}, {'out'});
input_dim = 16;

rng(49);
center1 = randn(input_dim, 1) * 0.5;
center2 = randn(input_dim, 1) * 0.5;
eps = 0.02;

[all_contained, n_samples] = run_soundness_check_addition(layer, center1, center2, eps, 500);

assert(all_contained, 'Soundness violation: 16D stress test');
fprintf('Test 11 PASSED: %d/%d samples contained\n', n_samples, n_samples);

%% Summary
fprintf('\n=== All AdditionLayer Soundness Tests PASSED ===\n');
fprintf('Total: 11 tests\n');

%% Helper Functions
function [all_contained, n_samples] = run_soundness_check_addition(layer, center1, center2, eps, n_samples)
    input_dim = length(center1);
    lb1 = center1 - eps;
    ub1 = center1 + eps;
    lb2 = center2 - eps;
    ub2 = center2 + eps;

    S_in1 = Star(lb1, ub1);
    S_in2 = Star(lb2, ub2);

    inputs = {S_in1, S_in2};
    S_out = layer.reach(inputs, 'approx-star');
    [out_lb, out_ub] = S_out.getRanges();
    out_lb = out_lb(:);
    out_ub = out_ub(:);

    all_contained = true;
    tol = 1e-5;

    for i = 1:n_samples
        x1_sample = lb1 + (ub1 - lb1) .* rand(input_dim, 1);
        x2_sample = lb2 + (ub2 - lb2) .* rand(input_dim, 1);
        y_sample = layer.evaluate({x1_sample, x2_sample});
        if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
            all_contained = false;
            fprintf('  Sample %d OUTSIDE bounds\n', i);
            break;
        end
    end
end

function [all_contained, n_samples] = run_soundness_check_addition_3(layer, center1, center2, center3, eps, n_samples)
    input_dim = length(center1);
    lb1 = center1 - eps;
    ub1 = center1 + eps;
    lb2 = center2 - eps;
    ub2 = center2 + eps;
    lb3 = center3 - eps;
    ub3 = center3 + eps;

    S_in1 = Star(lb1, ub1);
    S_in2 = Star(lb2, ub2);
    S_in3 = Star(lb3, ub3);

    inputs = {S_in1, S_in2, S_in3};
    S_out = layer.reach(inputs, 'approx-star');
    [out_lb, out_ub] = S_out.getRanges();
    out_lb = out_lb(:);
    out_ub = out_ub(:);

    all_contained = true;
    tol = 1e-5;

    for i = 1:n_samples
        x1_sample = lb1 + (ub1 - lb1) .* rand(input_dim, 1);
        x2_sample = lb2 + (ub2 - lb2) .* rand(input_dim, 1);
        x3_sample = lb3 + (ub3 - lb3) .* rand(input_dim, 1);
        y_sample = layer.evaluate({x1_sample, x2_sample, x3_sample});
        if any(y_sample < out_lb - tol) || any(y_sample > out_ub + tol)
            all_contained = false;
            fprintf('  Sample %d OUTSIDE bounds\n', i);
            break;
        end
    end
end
