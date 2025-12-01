% test_soundness_Zono
% Tests for Zonotope set representation
% To run: results = runtests('test_soundness_Zono')

%% Test 1: Zono constructor with center and generators
rng(42);

c = [1; 2; 3];  % center
V = randn(3, 2);  % 2 generators

Z = Zono(c, V);

assert(~isempty(Z), 'Zono should be created');
assert(Z.dim == 3, 'Dimension should be 3');

%% Test 2: Zono getBounds
rng(42);

c = [0; 0];
V = [1 0; 0 1];  % unit square centered at origin

Z = Zono(c, V);
[lb, ub] = Z.getBounds();

assert(all(abs(lb - [-1; -1]) < 1e-10), 'Lower bounds should be [-1; -1]');
assert(all(abs(ub - [1; 1]) < 1e-10), 'Upper bounds should be [1; 1]');

%% Test 3: Zono affineMap
rng(42);

c = [1; 1];
V = [0.5 0; 0 0.5];

Z = Zono(c, V);

% Apply affine transformation
W = [2 0; 0 2];
b = [1; 1];
Z2 = Z.affineMap(W, b);

assert(Z2.dim == 2, 'Transformed Zono should have dim 2');

% Check new center
expected_center = W * c + b;
assert(max(abs(Z2.c - expected_center)) < 1e-10, 'Center should transform correctly');

%% Test 4: Zono contains center
rng(42);

c = [0.5; 0.5; 0.5];
V = randn(3, 3) * 0.1;

Z = Zono(c, V);

% Center should always be in the zonotope (alpha = 0)
[lb, ub] = Z.getBounds();
assert(all(c >= lb - 1e-10) && all(c <= ub + 1e-10), ...
    'Center should be within bounds');

%% Test 5: Zono sample containment
rng(42);

c = [0; 0; 0];
V = [1 0 0; 0 1 0; 0 0 1];  % unit cube

Z = Zono(c, V);

% Sample with alpha in [-1, 1]
for i = 1:20
    alpha = -1 + 2*rand(3, 1);
    sample = c + V * alpha;
    [lb, ub] = Z.getBounds();
    assert(all(sample >= lb - 1e-10) && all(sample <= ub + 1e-10), ...
        'Sample %d should be within bounds', i);
end

%% Test 6: Zono from interval bounds
rng(42);

lb = [0; 0];
ub = [2; 4];

% Create Zono from bounds (center + generator form)
c = (lb + ub) / 2;  % [1; 2]
V = diag((ub - lb) / 2);  % [[1 0]; [0 2]]

Z = Zono(c, V);

[lb_check, ub_check] = Z.getBounds();
assert(max(abs(lb_check - lb)) < 1e-10, 'Lower bounds should match');
assert(max(abs(ub_check - ub)) < 1e-10, 'Upper bounds should match');

