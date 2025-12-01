% test_soundness_ConcatenationLayer
% Soundness tests for ConcatenationLayer
% Verifies that constraints from ALL inputs are preserved (block-diagonal)
% To run: results = runtests('test_soundness_ConcatenationLayer')

%% Test 1: Two-input concatenation preserves all constraints
% Create first ImageStar: 2x2x1 with one predicate variable
V1 = zeros(2, 2, 1, 2);
V1(:,:,1,1) = [0.5 0.5; 0.5 0.5];  % center
V1(:,:,1,2) = [1 0; 0 0];          % basis (affects pixel (1,1))
C1 = [1; -1];
d1 = [1; 1];  % alpha1 in [-1, 1]
pred_lb1 = -1;
pred_ub1 = 1;
IS1 = ImageStar(V1, C1, d1, pred_lb1, pred_ub1);

% Create second ImageStar: 2x2x1 with one predicate variable
V2 = zeros(2, 2, 1, 2);
V2(:,:,1,1) = [0 0; 0 0];          % center
V2(:,:,1,2) = [0 0; 0 1];          % basis (affects pixel (2,2))
C2 = [1; -1];
d2 = [0.5; 0.5];  % alpha2 in [-0.5, 0.5]
pred_lb2 = -0.5;
pred_ub2 = 0.5;
IS2 = ImageStar(V2, C2, d2, pred_lb2, pred_ub2);

% Create ConcatenationLayer
layer = ConcatenationLayer('test_concat', 3);  % concat along channel dim

% Perform concatenation
inputs = {IS1, IS2};
IS_concat = layer.reach_single_input(inputs);

% Verify structural properties (block-diagonal constraints)
assert(IS_concat.numPred == 2, 'Should have 2 predicate variables (1+1)');
assert(size(IS_concat.C, 1) == 4, 'Should have 4 constraint rows (2+2)');
assert(size(IS_concat.C, 2) == 2, 'Constraint matrix should have 2 columns');
assert(IS_concat.numChannel == 2, 'Should have 2 channels after concat');

%% Test 2: Soundness via sampling - all corners must be contained
% Recreate all variables (test sections are independent)
V1 = zeros(2, 2, 1, 2);
V1(:,:,1,1) = [0.5 0.5; 0.5 0.5];
V1(:,:,1,2) = [1 0; 0 0];
C1 = [1; -1];
d1 = [1; 1];
IS1 = ImageStar(V1, C1, d1, -1, 1);

V2 = zeros(2, 2, 1, 2);
V2(:,:,1,1) = [0 0; 0 0];
V2(:,:,1,2) = [0 0; 0 1];
C2 = [1; -1];
d2 = [0.5; 0.5];
IS2 = ImageStar(V2, C2, d2, -0.5, 0.5);

layer = ConcatenationLayer('test_concat', 3);
inputs = {IS1, IS2};
IS_concat = layer.reach_single_input(inputs);

% Sample all corner combinations
alpha1_vals = [-1, 0, 1];
alpha2_vals = [-0.5, 0, 0.5];

for a1 = alpha1_vals
    for a2 = alpha2_vals
        % Check constraint satisfaction for original inputs
        if all(C1 * a1 <= d1) && all(C2 * a2 <= d2)
            % Compute concrete images
            img1 = V1(:,:,:,1) + a1 * V1(:,:,:,2);
            img2 = V2(:,:,:,1) + a2 * V2(:,:,:,2);
            img_concat = cat(3, img1, img2);

            % Verify containment in output
            contained = soundness_test_utils.verify_imagestar_containment(IS_concat, img_concat);
            assert(contained, ...
                'Soundness violation: corner (a1=%.1f, a2=%.1f) not in output', a1, a2);
        end
    end
end

%% Test 3: Three-input concatenation
% Recreate first two ImageStars
V1 = zeros(2, 2, 1, 2);
V1(:,:,1,1) = [0.5 0.5; 0.5 0.5];
V1(:,:,1,2) = [1 0; 0 0];
IS1 = ImageStar(V1, [1; -1], [1; 1], -1, 1);

V2 = zeros(2, 2, 1, 2);
V2(:,:,1,1) = [0 0; 0 0];
V2(:,:,1,2) = [0 0; 0 1];
IS2 = ImageStar(V2, [1; -1], [0.5; 0.5], -0.5, 0.5);

% Create third ImageStar
V3 = zeros(2, 2, 1, 2);
V3(:,:,1,1) = ones(2, 2);
V3(1,2,1,2) = 0.5;
C3 = [1; -1];
d3 = [0.3; 0.3];
pred_lb3 = -0.3;
pred_ub3 = 0.3;
IS3 = ImageStar(V3, C3, d3, pred_lb3, pred_ub3);

layer = ConcatenationLayer('test_concat', 3);
inputs3 = {IS1, IS2, IS3};
IS_concat3 = layer.reach_single_input(inputs3);

% Verify structure
assert(IS_concat3.numPred == 3, 'Should have 3 predicate variables');
assert(size(IS_concat3.C, 1) == 6, 'Should have 6 constraint rows');
assert(IS_concat3.numChannel == 3, 'Should have 3 channels');

%% Test 4: Random sampling soundness for 3-input case
% Recreate all variables
V1 = zeros(2, 2, 1, 2);
V1(:,:,1,1) = [0.5 0.5; 0.5 0.5];
V1(:,:,1,2) = [1 0; 0 0];
pred_lb1 = -1;
pred_ub1 = 1;
IS1 = ImageStar(V1, [1; -1], [1; 1], pred_lb1, pred_ub1);

V2 = zeros(2, 2, 1, 2);
V2(:,:,1,1) = [0 0; 0 0];
V2(:,:,1,2) = [0 0; 0 1];
pred_lb2 = -0.5;
pred_ub2 = 0.5;
IS2 = ImageStar(V2, [1; -1], [0.5; 0.5], pred_lb2, pred_ub2);

V3 = zeros(2, 2, 1, 2);
V3(:,:,1,1) = ones(2, 2);
V3(1,2,1,2) = 0.5;
pred_lb3 = -0.3;
pred_ub3 = 0.3;
IS3 = ImageStar(V3, [1; -1], [0.3; 0.3], pred_lb3, pred_ub3);

layer = ConcatenationLayer('test_concat', 3);
inputs3 = {IS1, IS2, IS3};
IS_concat3 = layer.reach_single_input(inputs3);

n_samples = 50;
for i = 1:n_samples
    % Generate valid random alpha values
    a1 = pred_lb1 + (pred_ub1 - pred_lb1) * rand();
    a2 = pred_lb2 + (pred_ub2 - pred_lb2) * rand();
    a3 = pred_lb3 + (pred_ub3 - pred_lb3) * rand();

    % Compute concrete images
    img1 = V1(:,:,:,1) + a1 * V1(:,:,:,2);
    img2 = V2(:,:,:,1) + a2 * V2(:,:,:,2);
    img3 = V3(:,:,:,1) + a3 * V3(:,:,:,2);
    img_concat = cat(3, img1, img2, img3);

    % Verify containment
    contained = soundness_test_utils.verify_imagestar_containment(IS_concat3, img_concat);
    assert(contained, 'Soundness violation at random sample %d', i);
end

%% Test 5: Block-diagonal structure verification
% The constraint matrix should be block-diagonal
layer5 = ConcatenationLayer('test_blkdiag', 3);

V5a = zeros(2, 2, 1, 3);  % 2 predicate vars
V5a(:,:,1,1) = ones(2,2);
V5a(1,1,1,2) = 1; V5a(1,2,1,3) = 1;
C5a = [1 0; 0 1; -1 0; 0 -1];
d5a = [0.5; 0.5; 0.5; 0.5];
IS5a = ImageStar(V5a, C5a, d5a, [-0.5; -0.5], [0.5; 0.5]);

V5b = zeros(2, 2, 1, 4);  % 3 predicate vars
V5b(:,:,1,1) = 2*ones(2,2);
V5b(2,1,1,2) = 1; V5b(2,2,1,3) = 1; V5b(1,1,1,4) = 1;
C5b = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
d5b = [0.3; 0.3; 0.3; 0.3; 0.3; 0.3];
IS5b = ImageStar(V5b, C5b, d5b, [-0.3; -0.3; -0.3], [0.3; 0.3; 0.3]);

inputs5 = {IS5a, IS5b};
IS_concat5 = layer5.reach_single_input(inputs5);

% Verify block-diagonal structure
assert(IS_concat5.numPred == 5, 'Expected 5 predicate vars');
assert(size(IS_concat5.C, 1) == 10, 'Expected 10 constraints');

% Off-diagonal blocks should be zero
offDiag1 = IS_concat5.C(1:4, 3:5);  % rows from IS5a, cols from IS5b
offDiag2 = IS_concat5.C(5:10, 1:2); % rows from IS5b, cols from IS5a
assert(all(offDiag1(:) == 0), 'Off-diagonal block 1 should be zero');
assert(all(offDiag2(:) == 0), 'Off-diagonal block 2 should be zero');

%% Test 6: Concatenation along dimension 1 (height)
layer6 = ConcatenationLayer('concat_dim1', 1);

V6a = zeros(2, 3, 1, 2);
V6a(:,:,1,1) = ones(2,3);
V6a(1,1,1,2) = 0.5;
IS6a = ImageStar(V6a, [1;-1], [1;1], -1, 1);

V6b = zeros(3, 3, 1, 2);
V6b(:,:,1,1) = 2*ones(3,3);
V6b(2,2,1,2) = 0.3;
IS6b = ImageStar(V6b, [1;-1], [0.5;0.5], -0.5, 0.5);

inputs6 = {IS6a, IS6b};
IS_concat6 = layer6.reach_single_input(inputs6);

% Verify dimensions (should be 5x3x1 after concatenating 2x3 and 3x3 along height)
assert(size(IS_concat6.V, 1) == 5, 'Output height should be 2+3=5');
assert(IS_concat6.numPred == 2, 'Should have 2 predicate variables');

% Verify soundness with sampling
for i = 1:20
    a1 = -1 + 2*rand();
    a2 = -0.5 + rand();
    img1 = V6a(:,:,:,1) + a1 * V6a(:,:,:,2);
    img2 = V6b(:,:,:,1) + a2 * V6b(:,:,:,2);
    img_concat = cat(1, img1, img2);
    contained = soundness_test_utils.verify_imagestar_containment(IS_concat6, img_concat);
    assert(contained, 'Dim-1 concatenation soundness violation at sample %d', i);
end
