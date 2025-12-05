% test_soundness_Box
% Tests for Box (hyper-rectangle) set representation
% To run: results = runtests('test_soundness_Box')

%% Test 1: Box constructor
rng(42);

lb = [-1; -2; -3];
ub = [1; 2; 3];

B = Box(lb, ub);

assert(isequal(B.lb, lb), 'Lower bound should be preserved');
assert(isequal(B.ub, ub), 'Upper bound should be preserved');
assert(B.dim == 3, 'Dimension should be 3');

%% Test 2: Box center computation
rng(42);

lb = [0; 0];
ub = [4; 6];

B2 = Box(lb, ub);

expected_center = [2; 3];
assert(max(abs(B2.center - expected_center)) < 1e-10, 'Center should be midpoint');

%% Test 3: Box generators
rng(42);

lb = [-1; -1];
ub = [1; 1];

B3 = Box(lb, ub);

% Generators should be diagonal matrix with half-widths
% For [-1,1] x [-1,1], generators should be [1,0; 0,1]
assert(~isempty(B3.generators), 'Generators should be defined');
assert(size(B3.generators, 1) == 2, 'Should have 2-dimensional generators');

%% Test 4: Box with tight bounds (lb == ub)
rng(42);

center = [1; 2; 3];
lb = center;
ub = center;

B4 = Box(lb, ub);

assert(B4.dim == 3, 'Dimension should be 3');
assert(isequal(B4.center, center), 'Center should equal the point');

%% Test 5: Box single partition
rng(42);

lb = [0; 0];
ub = [4; 4];

B5 = Box(lb, ub);

% Partition along first dimension into 2 parts
Bs = B5.singlePartition(1, 2);

assert(length(Bs) == 2, 'Should produce 2 boxes');
assert(Bs(1).ub(1) == 2, 'First box should have ub(1) = 2');
assert(Bs(2).lb(1) == 2, 'Second box should have lb(1) = 2');

%% Test 6: Box partition multiple dimensions
rng(42);

lb = [0; 0];
ub = [4; 4];

B6 = Box(lb, ub);

% Partition both dimensions into 2 parts each
Bs = B6.partition([1, 2], [2, 2]);

assert(length(Bs) == 4, 'Should produce 4 boxes (2x2 partition)');

%% Test 7: Box sample containment
rng(42);

lb = [-2; -3; -1];
ub = [2; 3; 1];

B7 = Box(lb, ub);

% Samples within bounds should satisfy containment
for i = 1:10
    sample = lb + rand(size(lb)) .* (ub - lb);
    assert(all(sample >= lb - 1e-10), 'Sample should be >= lb');
    assert(all(sample <= ub + 1e-10), 'Sample should be <= ub');
end

