% test_soundness_AdditionLayer
% Soundness tests for AdditionLayer
% Verifies that computed reachable sets contain all actual outputs
% To run: results = runtests('test_soundness_AdditionLayer')

%% Test 1: Two-input addition preserves all constraints
L = AdditionLayer('add_test', 2, 1, {'in1', 'in2'}, {'out'});

% Create first ImageStar
V1 = zeros(3, 3, 1, 2);
V1(:,:,1,1) = ones(3,3);
V1(1,1,1,2) = 0.5;
C1 = [1; -1];
d1 = [1; 1];
IS1 = ImageStar(V1, C1, d1, -1, 1);

% Create second ImageStar
V2 = zeros(3, 3, 1, 2);
V2(:,:,1,1) = 2*ones(3,3);
V2(2,2,1,2) = 0.3;
C2 = [1; -1];
d2 = [0.5; 0.5];
IS2 = ImageStar(V2, C2, d2, -0.5, 0.5);

% Perform addition
inputs = {IS1, IS2};
IS_add = L.reach_single_input(inputs);

% Verify structure
assert(IS_add.numPred == 2, 'Should have 2 predicate variables');
assert(size(IS_add.C, 1) == 4, 'Should have 4 constraints');

%% Test 2: Soundness via sampling
% Recreate all variables (test sections are independent)
L = AdditionLayer('add_test', 2, 1, {'in1', 'in2'}, {'out'});

V1 = zeros(3, 3, 1, 2);
V1(:,:,1,1) = ones(3,3);
V1(1,1,1,2) = 0.5;
IS1 = ImageStar(V1, [1; -1], [1; 1], -1, 1);

V2 = zeros(3, 3, 1, 2);
V2(:,:,1,1) = 2*ones(3,3);
V2(2,2,1,2) = 0.3;
IS2 = ImageStar(V2, [1; -1], [0.5; 0.5], -0.5, 0.5);

inputs = {IS1, IS2};
IS_add = L.reach_single_input(inputs);

n_samples = 50;
for i = 1:n_samples
    a1 = -1 + 2*rand();
    a2 = -0.5 + rand();

    img1 = V1(:,:,:,1) + a1 * V1(:,:,:,2);
    img2 = V2(:,:,:,1) + a2 * V2(:,:,:,2);
    img_add = img1 + img2;

    contained = soundness_test_utils.verify_imagestar_containment(IS_add, img_add);
    assert(contained, 'Addition soundness violation at sample %d', i);
end

%% Test 3: Three-input addition
% Recreate all variables
L3 = AdditionLayer('add_3', 3, 1, {'in1', 'in2', 'in3'}, {'out'});

V1 = zeros(3, 3, 1, 2);
V1(:,:,1,1) = ones(3,3);
V1(1,1,1,2) = 0.5;
IS1 = ImageStar(V1, [1; -1], [1; 1], -1, 1);

V2 = zeros(3, 3, 1, 2);
V2(:,:,1,1) = 2*ones(3,3);
V2(2,2,1,2) = 0.3;
IS2 = ImageStar(V2, [1; -1], [0.5; 0.5], -0.5, 0.5);

V3 = zeros(3, 3, 1, 2);
V3(:,:,1,1) = -ones(3,3);
V3(3,3,1,2) = 0.2;
C3 = [1; -1];
d3 = [0.3; 0.3];
IS3 = ImageStar(V3, C3, d3, -0.3, 0.3);

inputs3 = {IS1, IS2, IS3};
IS_add3 = L3.reach_single_input(inputs3);

% Verify structure
assert(IS_add3.numPred == 3, 'Should have 3 predicate variables');

% Verify soundness
for i = 1:30
    a1 = -1 + 2*rand();
    a2 = -0.5 + rand();
    a3 = -0.3 + 0.6*rand();

    img1 = V1(:,:,:,1) + a1 * V1(:,:,:,2);
    img2 = V2(:,:,:,1) + a2 * V2(:,:,:,2);
    img3 = V3(:,:,:,1) + a3 * V3(:,:,:,2);
    img_add = img1 + img2 + img3;

    contained = soundness_test_utils.verify_imagestar_containment(IS_add3, img_add);
    assert(contained, '3-input addition soundness violation at sample %d', i);
end

%% Test 4: Multi-channel addition
L4 = AdditionLayer('add_mc', 2, 1, {'in1', 'in2'}, {'out'});

V4a = zeros(2, 2, 3, 2);
V4a(:,:,:,1) = rand(2,2,3);
V4a(1,1,1,2) = 0.1;
IS4a = ImageStar(V4a, [1;-1], [1;1], -1, 1);

V4b = zeros(2, 2, 3, 2);
V4b(:,:,:,1) = rand(2,2,3);
V4b(2,2,2,2) = 0.2;
IS4b = ImageStar(V4b, [1;-1], [0.5;0.5], -0.5, 0.5);

inputs4 = {IS4a, IS4b};
IS_add4 = L4.reach_single_input(inputs4);

% Verify soundness
for i = 1:30
    a1 = -1 + 2*rand();
    a2 = -0.5 + rand();

    img1 = V4a(:,:,:,1) + a1 * V4a(:,:,:,2);
    img2 = V4b(:,:,:,1) + a2 * V4b(:,:,:,2);
    img_add = img1 + img2;

    contained = soundness_test_utils.verify_imagestar_containment(IS_add4, img_add);
    assert(contained, 'Multi-channel addition soundness violation at sample %d', i);
end

%% Test 5: Corner case sampling
% Recreate all variables
L = AdditionLayer('add_test', 2, 1, {'in1', 'in2'}, {'out'});

V1 = zeros(3, 3, 1, 2);
V1(:,:,1,1) = ones(3,3);
V1(1,1,1,2) = 0.5;
IS1 = ImageStar(V1, [1; -1], [1; 1], -1, 1);

V2 = zeros(3, 3, 1, 2);
V2(:,:,1,1) = 2*ones(3,3);
V2(2,2,1,2) = 0.3;
IS2 = ImageStar(V2, [1; -1], [0.5; 0.5], -0.5, 0.5);

inputs = {IS1, IS2};
IS_add = L.reach_single_input(inputs);

% Explicitly test all 4 corners of the 2-predicate input space
corners = [-1 -0.5; -1 0.5; 1 -0.5; 1 0.5];

for i = 1:size(corners, 1)
    a1 = corners(i, 1);
    a2 = corners(i, 2);

    img1 = V1(:,:,:,1) + a1 * V1(:,:,:,2);
    img2 = V2(:,:,:,1) + a2 * V2(:,:,:,2);
    img_add = img1 + img2;

    contained = soundness_test_utils.verify_imagestar_containment(IS_add, img_add);
    assert(contained, 'Corner case %d not contained in output', i);
end
