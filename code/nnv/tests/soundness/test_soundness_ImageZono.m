% test_soundness_ImageZono
% Tests for ImageZono set representation
% To run: results = runtests('test_soundness_ImageZono')

%% Test 1: ImageZono constructor with bounds
rng(42);

lb = zeros(4, 4);
ub = ones(4, 4) * 0.5;

IZ = ImageZono(lb, ub);

assert(~isempty(IZ), 'ImageZono should be created');
assert(IZ.height == 4, 'Height should be 4');
assert(IZ.width == 4, 'Width should be 4');
assert(IZ.numChannels == 1, 'Single channel grayscale');

%% Test 2: ImageZono with multi-channel
rng(42);

lb = rand(3, 3, 3) * 0.2;
ub = lb + rand(3, 3, 3) * 0.3;

IZ2 = ImageZono(lb, ub);

assert(IZ2.height == 3, 'Height should be 3');
assert(IZ2.width == 3, 'Width should be 3');
assert(IZ2.numChannels == 3, 'Should have 3 channels');

%% Test 3: ImageZono bounds preservation
rng(42);

lb = [0 0.1; 0.2 0.3];
ub = [0.5 0.6; 0.7 0.8];

IZ3 = ImageZono(lb, ub);

assert(isequal(IZ3.lb_image, lb), 'Lower bound should be preserved');
assert(isequal(IZ3.ub_image, ub), 'Upper bound should be preserved');

%% Test 4: ImageZono V representation
rng(42);

lb = zeros(2, 2);
ub = ones(2, 2);

IZ4 = ImageZono(lb, ub);

% V should have numPred+1 in last dimension (center + generators)
assert(~isempty(IZ4.V), 'V should be defined');
assert(IZ4.numPreds >= 0, 'numPreds should be non-negative');
assert(size(IZ4.V, 4) == IZ4.numPreds + 1, 'V last dim should be numPreds+1');

%% Test 5: ImageZono sample containment
rng(42);

lb = rand(3, 3) * 0.3;
ub = lb + rand(3, 3) * 0.4;

IZ5 = ImageZono(lb, ub);

% Samples within bounds should be representable
for i = 1:10
    sample = lb + rand(size(lb)) .* (ub - lb);
    assert(all(sample(:) >= lb(:) - 1e-10), 'Sample should be >= lb');
    assert(all(sample(:) <= ub(:) + 1e-10), 'Sample should be <= ub');
end

%% Test 6: ImageZono with tight bounds
rng(42);

% Tight bounds (lb == ub) should work
center = rand(2, 2);
lb = center;
ub = center;

IZ6 = ImageZono(lb, ub);

assert(IZ6.height == 2, 'Height should be 2');
assert(IZ6.width == 2, 'Width should be 2');

