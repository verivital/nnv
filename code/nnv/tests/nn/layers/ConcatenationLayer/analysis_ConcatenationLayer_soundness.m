%% Test ConcatenationLayer soundness
% This test verifies that ConcatenationLayer properly handles the
% concatenation of multiple ImageStar sets with independent predicate variables.
%
% Issue: When concatenating two ImageStar sets from different branches,
% each set has its own predicate variables and constraints. The concatenated
% result must preserve ALL constraints from both inputs, not just one.
%
% The correct approach is to use block-diagonal constraint matrices
% (like Star.concatenate does), not to discard constraints from some inputs.
%
% Author: Test for soundness issue
% Date: 2024

%% Test 1: Basic soundness test - verify constraints are preserved
% Create two simple ImageStars with different constraints

fprintf('\n=== Test 1: Basic Constraint Preservation ===\n');

% Create first ImageStar: 2x2x1 image with one predicate variable
% x1 in [-1, 1], center at 0.5
V1 = zeros(2, 2, 1, 2);
V1(:,:,1,1) = [0.5 0.5; 0.5 0.5];  % center
V1(:,:,1,2) = [1 0; 0 0];          % basis vector (only affects pixel (1,1))
C1 = [1; -1];   % constraints: alpha1 <= 1 and -alpha1 <= 1
d1 = [1; 1];    % so alpha1 in [-1, 1]
pred_lb1 = -1;
pred_ub1 = 1;
IS1 = ImageStar(V1, C1, d1, pred_lb1, pred_ub1);

% Create second ImageStar: 2x2x1 image with one predicate variable
% x2 in [-0.5, 0.5], center at 0
V2 = zeros(2, 2, 1, 2);
V2(:,:,1,1) = [0 0; 0 0];          % center
V2(:,:,1,2) = [0 0; 0 1];          % basis vector (only affects pixel (2,2))
C2 = [1; -1];   % constraints: alpha2 <= 0.5 and -alpha2 <= 0.5
d2 = [0.5; 0.5];  % so alpha2 in [-0.5, 0.5]
pred_lb2 = -0.5;
pred_ub2 = 0.5;
IS2 = ImageStar(V2, C2, d2, pred_lb2, pred_ub2);

% Print input set info
fprintf('Input Set 1 (IS1):\n');
fprintf('  Center pixel (1,1): %.2f, range: [%.2f, %.2f]\n', ...
    V1(1,1,1,1), V1(1,1,1,1) + pred_lb1*V1(1,1,1,2), V1(1,1,1,1) + pred_ub1*V1(1,1,1,2));
fprintf('  nVar: %d, Constraints: %d\n', IS1.numPred, size(C1, 1));

fprintf('Input Set 2 (IS2):\n');
fprintf('  Center pixel (2,2): %.2f, range: [%.2f, %.2f]\n', ...
    V2(2,2,1,1), V2(2,2,1,1) + pred_lb2*V2(2,2,1,2), V2(2,2,1,1) + pred_ub2*V2(2,2,1,2));
fprintf('  nVar: %d, Constraints: %d\n', IS2.numPred, size(C2, 1));

% Create ConcatenationLayer (concatenate along channel dimension = 3)
layer = ConcatenationLayer('test_concat', 3);

% Perform concatenation via layer
inputs = {IS1, IS2};
IS_concat_layer = layer.reach_single_input(inputs);

% Perform concatenation via ImageStar.concatenation (correct method)
IS_concat_correct = IS1.concatenation(IS2);

fprintf('\nConcatenation Results:\n');
fprintf('Via ConcatenationLayer:\n');
fprintf('  nVar: %d\n', IS_concat_layer.numPred);
fprintf('  Constraint rows: %d\n', size(IS_concat_layer.C, 1));

fprintf('Via ImageStar.concatenation:\n');
fprintf('  nVar: %d\n', IS_concat_correct.numPred);
fprintf('  Constraint rows: %d\n', size(IS_concat_correct.C, 1));

%% Test 2: Soundness verification through sampling
fprintf('\n=== Test 2: Soundness Through Sampling ===\n');

% The key test: sample points from the input sets and verify they are
% contained in the output set

% Sample corner cases from IS1: alpha1 = -1, 0, 1
alpha1_samples = [-1, 0, 1];
% Sample corner cases from IS2: alpha2 = -0.5, 0, 0.5
alpha2_samples = [-0.5, 0, 0.5];

num_violations_layer = 0;
num_violations_correct = 0;
num_total = 0;

fprintf('Testing all combinations of predicate values...\n');

for a1 = alpha1_samples
    for a2 = alpha2_samples
        num_total = num_total + 1;

        % Compute actual images for these predicate values
        img1 = V1(:,:,:,1) + a1 * V1(:,:,:,2);
        img2 = V2(:,:,:,1) + a2 * V2(:,:,:,2);

        % Concatenate the actual images
        img_concat = cat(3, img1, img2);

        % Check if this point should be valid
        % IS1 constraint: C1 * a1 <= d1
        constraint1_satisfied = all(C1 * a1 <= d1 + 1e-9);
        % IS2 constraint: C2 * a2 <= d2
        constraint2_satisfied = all(C2 * a2 <= d2 + 1e-9);

        should_be_valid = constraint1_satisfied && constraint2_satisfied;

        % Check containment in layer result
        % Convert to Star and check
        S_layer = IS_concat_layer.toStar;
        pt = reshape(img_concat, [], 1);

        % For the layer result, we need to check if the point can be
        % expressed as V*alpha where C*alpha <= d
        % This is a feasibility problem

        if should_be_valid
            % Try to find alpha such that V*[1; alpha] = pt and C*alpha <= d
            % V = [c v1 v2 ... vn], so pt = c + V_basis * alpha
            c_layer = S_layer.V(:,1);
            V_basis_layer = S_layer.V(:,2:end);

            if size(V_basis_layer, 2) > 0
                % Solve: V_basis * alpha = pt - c, subject to C*alpha <= d
                % Check if solvable
                residual = pt - c_layer;

                % Use least squares to find alpha
                if rank(V_basis_layer) == size(V_basis_layer, 2)
                    alpha_solve = V_basis_layer \ residual;

                    % Check if reconstruction is accurate
                    reconstruction_error = norm(V_basis_layer * alpha_solve - residual);

                    % Check if constraints are satisfied
                    if ~isempty(S_layer.C)
                        constraints_ok = all(S_layer.C * alpha_solve <= S_layer.d + 1e-6);
                    else
                        constraints_ok = true;
                    end

                    if reconstruction_error > 1e-6 || ~constraints_ok
                        num_violations_layer = num_violations_layer + 1;
                        fprintf('  VIOLATION (layer): a1=%.1f, a2=%.1f - ', a1, a2);
                        if reconstruction_error > 1e-6
                            fprintf('reconstruction error=%.6f ', reconstruction_error);
                        end
                        if ~constraints_ok
                            fprintf('constraints violated');
                        end
                        fprintf('\n');
                    end
                else
                    fprintf('  WARNING: V_basis not full rank for layer result\n');
                end
            end

            % Check correct method
            S_correct = IS_concat_correct.toStar;
            c_correct = S_correct.V(:,1);
            V_basis_correct = S_correct.V(:,2:end);

            if size(V_basis_correct, 2) > 0
                residual = pt - c_correct;

                if rank(V_basis_correct) == size(V_basis_correct, 2)
                    alpha_solve = V_basis_correct \ residual;
                    reconstruction_error = norm(V_basis_correct * alpha_solve - residual);

                    if ~isempty(S_correct.C)
                        constraints_ok = all(S_correct.C * alpha_solve <= S_correct.d + 1e-6);
                    else
                        constraints_ok = true;
                    end

                    if reconstruction_error > 1e-6 || ~constraints_ok
                        num_violations_correct = num_violations_correct + 1;
                        fprintf('  VIOLATION (correct): a1=%.1f, a2=%.1f\n', a1, a2);
                    end
                end
            end
        end
    end
end

fprintf('\nResults:\n');
fprintf('  Total valid input combinations: %d\n', num_total);
fprintf('  Violations in ConcatenationLayer result: %d\n', num_violations_layer);
fprintf('  Violations in ImageStar.concatenation result: %d\n', num_violations_correct);

%% Test 3: Constraint structure comparison
fprintf('\n=== Test 3: Constraint Structure Analysis ===\n');

% The correct concatenation should have:
% - nVar = nVar1 + nVar2 (2 predicate variables total)
% - Block diagonal constraint matrix
% - Combined predicate bounds

fprintf('Expected (correct) structure:\n');
fprintf('  nVar should be: %d + %d = %d\n', IS1.numPred, IS2.numPred, IS1.numPred + IS2.numPred);
fprintf('  Constraint rows should be: %d + %d = %d\n', size(C1,1), size(C2,1), size(C1,1)+size(C2,1));

fprintf('\nActual ConcatenationLayer structure:\n');
fprintf('  nVar: %d (expected: %d) - %s\n', IS_concat_layer.numPred, IS1.numPred + IS2.numPred, ...
    iif(IS_concat_layer.numPred == IS1.numPred + IS2.numPred, 'OK', 'WRONG'));
fprintf('  Constraint rows: %d (expected: %d) - %s\n', size(IS_concat_layer.C, 1), size(C1,1)+size(C2,1), ...
    iif(size(IS_concat_layer.C, 1) == size(C1,1)+size(C2,1), 'OK', 'WRONG'));

fprintf('\nActual ImageStar.concatenation structure:\n');
fprintf('  nVar: %d (expected: %d) - %s\n', IS_concat_correct.numPred, IS1.numPred + IS2.numPred, ...
    iif(IS_concat_correct.numPred == IS1.numPred + IS2.numPred, 'OK', 'WRONG'));
fprintf('  Constraint rows: %d (expected: %d) - %s\n', size(IS_concat_correct.C, 1), size(C1,1)+size(C2,1), ...
    iif(size(IS_concat_correct.C, 1) == size(C1,1)+size(C2,1), 'OK', 'WRONG'));

%% Test 4: Demonstrate the soundness problem with bounds
fprintf('\n=== Test 4: Output Bounds Comparison ===\n');

% Get ranges from both methods
[lb_layer, ub_layer] = IS_concat_layer.getRanges;
[lb_correct, ub_correct] = IS_concat_correct.getRanges;

fprintf('Channel 1 (from IS1), pixel (1,1):\n');
fprintf('  Expected range: [%.2f, %.2f]\n', -0.5, 1.5);
fprintf('  Layer result:   [%.4f, %.4f]\n', lb_layer(1,1,1), ub_layer(1,1,1));
fprintf('  Correct result: [%.4f, %.4f]\n', lb_correct(1,1,1), ub_correct(1,1,1));

fprintf('Channel 2 (from IS2), pixel (2,2):\n');
fprintf('  Expected range: [%.2f, %.2f]\n', -0.5, 0.5);
fprintf('  Layer result:   [%.4f, %.4f]\n', lb_layer(2,2,2), ub_layer(2,2,2));
fprintf('  Correct result: [%.4f, %.4f]\n', lb_correct(2,2,2), ub_correct(2,2,2));

% Check if bounds are correct
bounds_correct_layer = abs(lb_layer(1,1,1) - (-0.5)) < 0.01 && abs(ub_layer(1,1,1) - 1.5) < 0.01 && ...
                       abs(lb_layer(2,2,2) - (-0.5)) < 0.01 && abs(ub_layer(2,2,2) - 0.5) < 0.01;
bounds_correct_correct = abs(lb_correct(1,1,1) - (-0.5)) < 0.01 && abs(ub_correct(1,1,1) - 1.5) < 0.01 && ...
                         abs(lb_correct(2,2,2) - (-0.5)) < 0.01 && abs(ub_correct(2,2,2) - 0.5) < 0.01;

fprintf('\nBounds correctness:\n');
fprintf('  ConcatenationLayer: %s\n', iif(bounds_correct_layer, 'CORRECT', 'INCORRECT'));
fprintf('  ImageStar.concatenation: %s\n', iif(bounds_correct_correct, 'CORRECT', 'INCORRECT'));

%% Summary
fprintf('\n========================================\n');
fprintf('SOUNDNESS TEST SUMMARY\n');
fprintf('========================================\n');

% Determine if there's a soundness issue
has_soundness_issue = (IS_concat_layer.numPred ~= IS1.numPred + IS2.numPred) || ...
                      (size(IS_concat_layer.C, 1) ~= size(C1,1) + size(C2,1));

if has_soundness_issue
    fprintf('RESULT: SOUNDNESS ISSUE DETECTED\n\n');
    fprintf('The ConcatenationLayer.reach_single_input method does not properly\n');
    fprintf('combine constraints from multiple input ImageStar sets.\n\n');
    fprintf('Current behavior:\n');
    fprintf('  - Only uses constraints from input with largest V matrix\n');
    fprintf('  - Discards constraints from other inputs\n\n');
    fprintf('Correct behavior (as in ImageStar.concatenation):\n');
    fprintf('  - Use block-diagonal constraint matrix: blkdiag(C1, C2, ...)\n');
    fprintf('  - Concatenate constraint vectors: [d1; d2; ...]\n');
    fprintf('  - Concatenate predicate bounds\n');
    fprintf('  - Properly handle independent predicate variables\n\n');
    fprintf('Recommendation: Fix ConcatenationLayer.reach_single_input to use\n');
    fprintf('ImageStar.concatenation method (like DepthConcatenationLayer does).\n');
else
    fprintf('RESULT: No soundness issue detected\n');
end

%% Helper function
function result = iif(condition, true_val, false_val)
    if condition
        result = true_val;
    else
        result = false_val;
    end
end
