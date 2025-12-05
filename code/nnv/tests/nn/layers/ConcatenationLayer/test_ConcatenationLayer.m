% Test ConcatenationLayer functionality
% To run: results = runtests('test_ConcatenationLayer')

%% Test 1: ConcatenationLayer constructor - with name and dim
L = ConcatenationLayer('concat1', 1);
assert(strcmp(L.Name, 'concat1'));
assert(L.Dim == 1);

%% Test 2: ConcatenationLayer constructor - with all parameters
L = ConcatenationLayer('concat2', 2, 1, {'in1', 'in2'}, {'out'}, 2);
assert(strcmp(L.Name, 'concat2'));
assert(L.Dim == 2);
assert(L.NumInputs == 2);

%% Test 3: ConcatenationLayer evaluate - concatenate along dim 1 (rows)
L = ConcatenationLayer('concat_dim1', 1);

% Create two matrices to concatenate
input1 = [1 2 3; 4 5 6];
input2 = [7 8 9; 10 11 12];

inputs = {input1, input2};
output = L.evaluate(inputs);

% Should concatenate vertically (along rows)
expected = [1 2 3; 4 5 6; 7 8 9; 10 11 12];
assert(isequal(output, expected), 'Concatenation along dim 1 failed');

%% Test 4: ConcatenationLayer evaluate - concatenate along dim 2 (columns)
L = ConcatenationLayer('concat_dim2', 2);

% Create two matrices to concatenate
input1 = [1 2; 3 4];
input2 = [5 6; 7 8];

inputs = {input1, input2};
output = L.evaluate(inputs);

% Should concatenate horizontally (along columns)
expected = [1 2 5 6; 3 4 7 8];
assert(isequal(output, expected), 'Concatenation along dim 2 failed');

%% Test 5: ConcatenationLayer evaluate - three inputs
L = ConcatenationLayer('concat_three', 1);

% Create three matrices
input1 = [1; 2];
input2 = [3; 4];
input3 = [5; 6];

inputs = {input1, input2, input3};
output = L.evaluate(inputs);

% Should concatenate all three vertically
expected = [1; 2; 3; 4; 5; 6];
assert(isequal(output, expected), 'Concatenation of three inputs failed');

%% Test 6: ConcatenationLayer evaluate - 3D arrays along dim 3
L = ConcatenationLayer('concat_3d', 3);

% Create two 3D arrays
input1(:,:,1) = [1 2; 3 4];
input1(:,:,2) = [5 6; 7 8];

input2(:,:,1) = [9 10; 11 12];

inputs = {input1, input2};
output = L.evaluate(inputs);

% Should concatenate along 3rd dimension
assert(size(output, 1) == 2, 'Dim 1 should be unchanged');
assert(size(output, 2) == 2, 'Dim 2 should be unchanged');
assert(size(output, 3) == 3, 'Dim 3 should be sum of inputs');

%% Test 8: ConcatenationLayer reach with ImageStar inputs
L = ConcatenationLayer('concat_imagestar', 3);

% Create two ImageStar inputs
IM1(:,:,1) = [1 1; 0 1];
LB1(:,:,1) = [-0.1 -0.1; 0 0];
UB1(:,:,1) = [0.1 0.1; 0 0];
image_star1 = ImageStar(IM1, LB1, UB1);

IM2(:,:,1) = [0 1; 1 0];
LB2(:,:,1) = [-0.1 -0.1; 0 0];
UB2(:,:,1) = [0.1 0.1; 0 0];
image_star2 = ImageStar(IM2, LB2, UB2);

inputs = {image_star1, image_star2};
output = L.reach(inputs, 'approx-star');

assert(isa(output, 'ImageStar'), 'reach should return ImageStar');

%% Test 9: ConcatenationLayer toGPU
L = ConcatenationLayer('concat', 1);
L_gpu = L.toGPU();
assert(isa(L_gpu, 'ConcatenationLayer'));

%% Test 10: ConcatenationLayer changeParamsPrecision
L = ConcatenationLayer('concat', 1);
L_single = L.changeParamsPrecision('single');
assert(isa(L_single, 'ConcatenationLayer'));

%% Test 11: REGRESSION - constraint preservation with independent ImageStars
% This test verifies that constraints from ALL input ImageStars are preserved
% when concatenating. This is critical for soundness.
% Issue: When concatenating two ImageStar sets from different branches,
% each has its own predicate variables. The result must preserve ALL constraints.

L = ConcatenationLayer('concat_soundness', 3);

% Create first ImageStar with specific constraints
V1 = zeros(2, 2, 1, 2);
V1(:,:,1,1) = [0.5 0.5; 0.5 0.5];  % center
V1(:,:,1,2) = [1 0; 0 0];          % basis vector (affects pixel (1,1))
C1 = [1; -1];
d1 = [1; 1];  % alpha1 in [-1, 1]
IS1 = ImageStar(V1, C1, d1, -1, 1);

% Create second ImageStar with DIFFERENT constraints
V2 = zeros(2, 2, 1, 2);
V2(:,:,1,1) = [0 0; 0 0];          % center
V2(:,:,1,2) = [0 0; 0 1];          % basis vector (affects pixel (2,2))
C2 = [1; -1];
d2 = [0.5; 0.5];  % alpha2 in [-0.5, 0.5] (tighter constraint)
IS2 = ImageStar(V2, C2, d2, -0.5, 0.5);

% Concatenate
inputs = {IS1, IS2};
IS_concat = L.reach_single_input(inputs);

% SOUNDNESS CHECK: Must have constraints from BOTH inputs
expected_nVar = IS1.numPred + IS2.numPred;  % 2
expected_constraints = size(C1,1) + size(C2,1);  % 4

assert(IS_concat.numPred == expected_nVar, ...
    sprintf('SOUNDNESS BUG: nVar=%d, expected %d. Constraints from some inputs lost!', ...
    IS_concat.numPred, expected_nVar));

assert(size(IS_concat.C, 1) == expected_constraints, ...
    sprintf('SOUNDNESS BUG: constraints=%d, expected %d. Constraints from some inputs lost!', ...
    size(IS_concat.C, 1), expected_constraints));

%% Test 12: REGRESSION - bounds correctness with independent constraints
% Verify that output bounds respect constraints from ALL inputs

L = ConcatenationLayer('concat_bounds_check', 3);

% Same setup as test 11
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

inputs = {IS1, IS2};
IS_concat = L.reach_single_input(inputs);

% Get bounds
[lb, ub] = IS_concat.getRanges;

% Channel 1, pixel (1,1): range should be [-0.5, 1.5] (from IS1: center 0.5 +/- 1)
assert(abs(lb(1,1,1) - (-0.5)) < 0.01, ...
    sprintf('SOUNDNESS BUG: lb(1,1,1)=%.2f, expected -0.5', lb(1,1,1)));
assert(abs(ub(1,1,1) - 1.5) < 0.01, ...
    sprintf('SOUNDNESS BUG: ub(1,1,1)=%.2f, expected 1.5', ub(1,1,1)));

% Channel 2, pixel (2,2): range should be [-0.5, 0.5] (from IS2: center 0 +/- 0.5)
% If constraints from IS2 are lost, this would incorrectly be [-1, 1]
assert(abs(lb(2,2,2) - (-0.5)) < 0.01, ...
    sprintf('SOUNDNESS BUG: lb(2,2,2)=%.2f, expected -0.5. IS2 constraints lost!', lb(2,2,2)));
assert(abs(ub(2,2,2) - 0.5) < 0.01, ...
    sprintf('SOUNDNESS BUG: ub(2,2,2)=%.2f, expected 0.5. IS2 constraints lost!', ub(2,2,2)));

%% Test 13: REGRESSION - three-way concatenation preserves all constraints
L = ConcatenationLayer('concat_three_soundness', 3);

V1 = zeros(2, 2, 1, 2);
V1(:,:,1,1) = ones(2,2);
V1(:,:,1,2) = [1 0; 0 0];
IS1 = ImageStar(V1, [1;-1], [0.3; 0.3], -0.3, 0.3);

V2 = zeros(2, 2, 1, 2);
V2(:,:,1,1) = 2*ones(2,2);
V2(:,:,1,2) = [0 1; 0 0];
IS2 = ImageStar(V2, [1;-1], [0.2; 0.2], -0.2, 0.2);

V3 = zeros(2, 2, 1, 2);
V3(:,:,1,1) = 3*ones(2,2);
V3(:,:,1,2) = [0 0; 1 0];
IS3 = ImageStar(V3, [1;-1], [0.1; 0.1], -0.1, 0.1);

inputs = {IS1, IS2, IS3};
IS_concat = L.reach_single_input(inputs);

% Must have 3 predicate variables and 6 constraints
assert(IS_concat.numPred == 3, ...
    sprintf('SOUNDNESS BUG: 3-way concat nVar=%d, expected 3', IS_concat.numPred));
assert(size(IS_concat.C, 1) == 6, ...
    sprintf('SOUNDNESS BUG: 3-way concat constraints=%d, expected 6', size(IS_concat.C,1)));
assert(IS_concat.numChannel == 3, '3-way concat should have 3 channels');

%% Test 14: HIGH-DIMENSIONAL - 5 inputs concatenation
% Verify soundness with 5 independent ImageStar inputs
L = ConcatenationLayer('concat_5way', 3);
numInputs = 5;
inputs = cell(1, numInputs);
expectedNVar = 0;
expectedConstraints = 0;

for i = 1:numInputs
    V = zeros(2, 2, 1, 2);
    V(:,:,1,1) = i * ones(2, 2);  % center varies by input
    V(:,:,1,2) = zeros(2, 2);
    V(mod(i-1,2)+1, mod(floor((i-1)/2),2)+1, 1, 2) = 1;  % different pixel affected
    bound = 0.1 * i;  % different bounds per input
    C = [1; -1];
    d = [bound; bound];
    inputs{i} = ImageStar(V, C, d, -bound, bound);
    expectedNVar = expectedNVar + 1;
    expectedConstraints = expectedConstraints + 2;
end

IS_concat = L.reach_single_input(inputs);

assert(IS_concat.numPred == expectedNVar, ...
    sprintf('5-way concat: nVar=%d, expected %d', IS_concat.numPred, expectedNVar));
assert(size(IS_concat.C, 1) == expectedConstraints, ...
    sprintf('5-way concat: constraints=%d, expected %d', size(IS_concat.C,1), expectedConstraints));
assert(IS_concat.numChannel == numInputs, ...
    sprintf('5-way concat: channels=%d, expected %d', IS_concat.numChannel, numInputs));

%% Test 15: HIGH-DIMENSIONAL - 9 inputs concatenation
L = ConcatenationLayer('concat_9way', 3);
numInputs = 9;
inputs = cell(1, numInputs);
expectedNVar = 0;
expectedConstraints = 0;

for i = 1:numInputs
    V = zeros(3, 3, 1, 2);
    V(:,:,1,1) = i * ones(3, 3);
    V(:,:,1,2) = zeros(3, 3);
    V(mod(i-1,3)+1, mod(floor((i-1)/3),3)+1, 1, 2) = 1;
    bound = 0.05 * i;
    C = [1; -1];
    d = [bound; bound];
    inputs{i} = ImageStar(V, C, d, -bound, bound);
    expectedNVar = expectedNVar + 1;
    expectedConstraints = expectedConstraints + 2;
end

IS_concat = L.reach_single_input(inputs);

assert(IS_concat.numPred == expectedNVar, ...
    sprintf('9-way concat: nVar=%d, expected %d', IS_concat.numPred, expectedNVar));
assert(size(IS_concat.C, 1) == expectedConstraints, ...
    sprintf('9-way concat: constraints=%d, expected %d', size(IS_concat.C,1), expectedConstraints));
assert(IS_concat.numChannel == numInputs, ...
    sprintf('9-way concat: channels=%d, expected %d', IS_concat.numChannel, numInputs));

%% Test 16: HIGH-DIMENSIONAL - 11 inputs concatenation
L = ConcatenationLayer('concat_11way', 3);
numInputs = 11;
inputs = cell(1, numInputs);
expectedNVar = 0;
expectedConstraints = 0;

for i = 1:numInputs
    V = zeros(4, 4, 1, 2);
    V(:,:,1,1) = i * ones(4, 4);
    V(:,:,1,2) = zeros(4, 4);
    V(mod(i-1,4)+1, mod(floor((i-1)/4),4)+1, 1, 2) = 1;
    bound = 0.02 * i;
    C = [1; -1];
    d = [bound; bound];
    inputs{i} = ImageStar(V, C, d, -bound, bound);
    expectedNVar = expectedNVar + 1;
    expectedConstraints = expectedConstraints + 2;
end

IS_concat = L.reach_single_input(inputs);

assert(IS_concat.numPred == expectedNVar, ...
    sprintf('11-way concat: nVar=%d, expected %d', IS_concat.numPred, expectedNVar));
assert(size(IS_concat.C, 1) == expectedConstraints, ...
    sprintf('11-way concat: constraints=%d, expected %d', size(IS_concat.C,1), expectedConstraints));
assert(IS_concat.numChannel == numInputs, ...
    sprintf('11-way concat: channels=%d, expected %d', IS_concat.numChannel, numInputs));

%% Test 17: HIGH-DIMENSIONAL - 16 inputs concatenation
L = ConcatenationLayer('concat_16way', 3);
numInputs = 16;
inputs = cell(1, numInputs);
expectedNVar = 0;
expectedConstraints = 0;

for i = 1:numInputs
    V = zeros(4, 4, 1, 2);
    V(:,:,1,1) = i * ones(4, 4);
    V(:,:,1,2) = zeros(4, 4);
    V(mod(i-1,4)+1, mod(floor((i-1)/4),4)+1, 1, 2) = 1;
    bound = 0.01 * i;
    C = [1; -1];
    d = [bound; bound];
    inputs{i} = ImageStar(V, C, d, -bound, bound);
    expectedNVar = expectedNVar + 1;
    expectedConstraints = expectedConstraints + 2;
end

IS_concat = L.reach_single_input(inputs);

assert(IS_concat.numPred == expectedNVar, ...
    sprintf('16-way concat: nVar=%d, expected %d', IS_concat.numPred, expectedNVar));
assert(size(IS_concat.C, 1) == expectedConstraints, ...
    sprintf('16-way concat: constraints=%d, expected %d', size(IS_concat.C,1), expectedConstraints));
assert(IS_concat.numChannel == numInputs, ...
    sprintf('16-way concat: channels=%d, expected %d', IS_concat.numChannel, numInputs));

%% Test 18: HIGH-DIMENSIONAL - 50 inputs concatenation
L = ConcatenationLayer('concat_50way', 3);
numInputs = 50;
inputs = cell(1, numInputs);
expectedNVar = 0;
expectedConstraints = 0;

for i = 1:numInputs
    V = zeros(8, 8, 1, 2);
    V(:,:,1,1) = i * ones(8, 8);
    V(:,:,1,2) = zeros(8, 8);
    V(mod(i-1,8)+1, mod(floor((i-1)/8),8)+1, 1, 2) = 1;
    bound = 0.005 * i;
    C = [1; -1];
    d = [bound; bound];
    inputs{i} = ImageStar(V, C, d, -bound, bound);
    expectedNVar = expectedNVar + 1;
    expectedConstraints = expectedConstraints + 2;
end

IS_concat = L.reach_single_input(inputs);

assert(IS_concat.numPred == expectedNVar, ...
    sprintf('50-way concat: nVar=%d, expected %d', IS_concat.numPred, expectedNVar));
assert(size(IS_concat.C, 1) == expectedConstraints, ...
    sprintf('50-way concat: constraints=%d, expected %d', size(IS_concat.C,1), expectedConstraints));
assert(IS_concat.numChannel == numInputs, ...
    sprintf('50-way concat: channels=%d, expected %d', IS_concat.numChannel, numInputs));

%% Test 19: HIGH-DIMENSIONAL - Multiple predicate variables per input
% Test with inputs that each have multiple predicate variables
L = ConcatenationLayer('concat_multi_pred', 3);
numInputs = 5;
inputs = cell(1, numInputs);
expectedNVar = 0;
expectedConstraints = 0;

for i = 1:numInputs
    nPred = i;  % input i has i predicate variables
    V = zeros(4, 4, 1, nPred + 1);
    V(:,:,1,1) = i * ones(4, 4);  % center
    for p = 1:nPred
        V(mod(p-1,4)+1, mod(floor((p-1)/4),4)+1, 1, p+1) = 0.1;
    end
    % Create constraints: each predicate in [-0.1*i, 0.1*i]
    bound = 0.1 * i;
    C = [eye(nPred); -eye(nPred)];
    d = bound * ones(2*nPred, 1);
    pred_lb = -bound * ones(nPred, 1);
    pred_ub = bound * ones(nPred, 1);
    inputs{i} = ImageStar(V, C, d, pred_lb, pred_ub);
    expectedNVar = expectedNVar + nPred;
    expectedConstraints = expectedConstraints + 2*nPred;
end

IS_concat = L.reach_single_input(inputs);

% Total: 1+2+3+4+5 = 15 predicate variables, 30 constraints
assert(IS_concat.numPred == expectedNVar, ...
    sprintf('Multi-pred concat: nVar=%d, expected %d', IS_concat.numPred, expectedNVar));
assert(size(IS_concat.C, 1) == expectedConstraints, ...
    sprintf('Multi-pred concat: constraints=%d, expected %d', size(IS_concat.C,1), expectedConstraints));
assert(IS_concat.numChannel == numInputs, ...
    sprintf('Multi-pred concat: channels=%d, expected %d', IS_concat.numChannel, numInputs));

%% Test 20: HIGH-DIMENSIONAL - Verify block-diagonal structure
% Ensure constraints remain properly isolated (block-diagonal)
L = ConcatenationLayer('concat_blkdiag_verify', 3);

% Create 3 inputs with distinctive constraint matrices
V1 = zeros(2, 2, 1, 3);  % 2 predicate vars
V1(:,:,1,1) = ones(2,2);
V1(1,1,1,2) = 1; V1(1,2,1,3) = 1;
C1 = [1 0; 0 1; -1 0; 0 -1];  % box constraints
d1 = [0.5; 0.5; 0.5; 0.5];
IS1 = ImageStar(V1, C1, d1, [-0.5; -0.5], [0.5; 0.5]);

V2 = zeros(2, 2, 1, 4);  % 3 predicate vars
V2(:,:,1,1) = 2*ones(2,2);
V2(2,1,1,2) = 1; V2(2,2,1,3) = 1; V2(1,1,1,4) = 1;
C2 = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
d2 = [0.3; 0.3; 0.3; 0.3; 0.3; 0.3];
IS2 = ImageStar(V2, C2, d2, [-0.3; -0.3; -0.3], [0.3; 0.3; 0.3]);

inputs = {IS1, IS2};
IS_concat = L.reach_single_input(inputs);

% Verify structure
assert(IS_concat.numPred == 5, 'Block-diag test: expected 5 predicate vars');
assert(size(IS_concat.C, 1) == 10, 'Block-diag test: expected 10 constraints');
assert(size(IS_concat.C, 2) == 5, 'Block-diag test: C should have 5 columns');

% Verify block-diagonal structure: C should have zeros in off-diagonal blocks
% Block 1: rows 1-4, cols 1-2 (from IS1)
% Block 2: rows 5-10, cols 3-5 (from IS2)
% Off-diagonal: rows 1-4, cols 3-5 should be zero; rows 5-10, cols 1-2 should be zero
offDiag1 = IS_concat.C(1:4, 3:5);
offDiag2 = IS_concat.C(5:10, 1:2);
assert(all(offDiag1(:) == 0), 'Block-diag test: off-diagonal block 1 should be zero');
assert(all(offDiag2(:) == 0), 'Block-diag test: off-diagonal block 2 should be zero');
