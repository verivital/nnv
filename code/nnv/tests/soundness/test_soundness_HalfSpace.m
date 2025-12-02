% test_soundness_HalfSpace
% Tests for HalfSpace set representation
% HalfSpace: {x | G*x <= g}
% To run: results = runtests('test_soundness_HalfSpace')

%% Test 1: HalfSpace constructor
rng(42);

G = [1 0; 0 1; -1 0; 0 -1];  % Box constraints
g = [1; 1; 0; 0];  % 0 <= x <= 1, 0 <= y <= 1

H = HalfSpace(G, g);

assert(isequal(H.G, G), 'G should match');
assert(isequal(H.g, g), 'g should match');

%% Test 2: HalfSpace dimension
rng(42);

G = [1 0 0; 0 1 0];  % 2 constraints, 3 dimensions
g = [1; 1];

H2 = HalfSpace(G, g);

assert(H2.dim == 3, 'Dimension should be 3');

%% Test 3: HalfSpace contains - point inside
rng(42);

G = [1 0; 0 1; -1 0; 0 -1];
g = [1; 1; 0; 0];  % 0 <= x <= 1, 0 <= y <= 1

H3 = HalfSpace(G, g);

% Point inside the box
x_inside = [0.5; 0.5];
assert(H3.contains(x_inside), 'Point (0.5, 0.5) should be inside');

%% Test 4: HalfSpace contains - point outside
rng(42);

G = [1 0; 0 1; -1 0; 0 -1];
g = [1; 1; 0; 0];  % 0 <= x <= 1, 0 <= y <= 1

H4 = HalfSpace(G, g);

% Point outside the box
x_outside = [2; 0.5];
assert(~H4.contains(x_outside), 'Point (2, 0.5) should be outside');

%% Test 5: HalfSpace contains - point on boundary
rng(42);

G = [1 0; 0 1; -1 0; 0 -1];
g = [1; 1; 0; 0];  % 0 <= x <= 1, 0 <= y <= 1

H5 = HalfSpace(G, g);

% Point on boundary
x_boundary = [1; 0.5];
assert(H5.contains(x_boundary), 'Point (1, 0.5) should be on boundary (inside)');

%% Test 6: HalfSpace single constraint
rng(42);

% Single hyperplane: x1 + x2 <= 1
G = [1 1];
g = [1];

H6 = HalfSpace(G, g);

assert(H6.dim == 2, 'Dimension should be 2');
assert(H6.contains([0; 0]), 'Origin should satisfy x1+x2 <= 1');
assert(~H6.contains([1; 1]), 'Point (1,1) should not satisfy x1+x2 <= 1');

