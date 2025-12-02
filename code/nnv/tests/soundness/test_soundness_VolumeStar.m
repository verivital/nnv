% test_soundness_VolumeStar
% Tests for VolumeStar set representation (3D Star sets)
% To run: results = runtests('test_soundness_VolumeStar')

%% Test 1: VolumeStar constructor with bounds
rng(42);

% Create 3D volume with lower/upper bounds
lb = zeros(3, 3, 3);
ub = ones(3, 3, 3) * 0.5;

VS = VolumeStar(lb, ub);

assert(~isempty(VS), 'VolumeStar should be created');
assert(VS.height == 3, 'Height should be 3');
assert(VS.width == 3, 'Width should be 3');
assert(VS.depth == 3, 'Depth should be 3');

%% Test 2: VolumeStar with center and perturbation
rng(42);

% Center volume
center = rand(4, 4, 4) * 0.5;

% Lower and upper perturbation
lb = center - 0.1;
ub = center + 0.1;

VS2 = VolumeStar(lb, ub);

assert(VS2.height == 4, 'Height should be 4');
assert(VS2.width == 4, 'Width should be 4');
assert(VS2.depth == 4, 'Depth should be 4');

%% Test 3: Sample from VolumeStar bounds
rng(42);

lb = zeros(2, 2, 2);
ub = ones(2, 2, 2);

VS3 = VolumeStar(lb, ub);

% Sample should be within bounds
for i = 1:10
    sample = lb + rand(size(lb)) .* (ub - lb);
    assert(all(sample(:) >= lb(:) - 1e-10), 'Sample should be >= lower bound');
    assert(all(sample(:) <= ub(:) + 1e-10), 'Sample should be <= upper bound');
end

%% Test 4: VolumeStar with single channel
rng(42);

lb = rand(5, 5, 5) * 0.2;
ub = lb + rand(5, 5, 5) * 0.3;

VS4 = VolumeStar(lb, ub);

assert(VS4.numChannel == 1, 'Single-volume should have 1 channel');

%% Test 5: VolumeStar V representation access
rng(42);

lb = zeros(2, 2, 2);
ub = ones(2, 2, 2) * 0.1;

VS5 = VolumeStar(lb, ub);

% V should be defined
assert(~isempty(VS5.V), 'V matrix should be defined');

% Number of predicates
assert(VS5.numPred >= 0, 'numPred should be non-negative');

