% Test SignLayer functionality
% To run: results = runtests('test_SignLayer')

%% Test 1: SignLayer constructor - default
L = SignLayer();
assert(L.gamma == 0);
assert(strcmp(L.mode, 'polar_zero_to_pos_one'));

%% Test 2: SignLayer constructor - with gamma
L = SignLayer(0.5);
assert(L.gamma == 0.5);
assert(strcmp(L.mode, 'polar_zero_to_pos_one'));

%% Test 3: SignLayer constructor - with gamma and mode
L = SignLayer(0, 'nonnegative_zero_to_pos_one');
assert(L.gamma == 0);
assert(strcmp(L.mode, 'nonnegative_zero_to_pos_one'));

%% Test 4: SignLayer evaluate - polar_zero_to_pos_one mode
L = SignLayer(0, 'polar_zero_to_pos_one');

% Test with simple input
input = [-2; -1; 0; 1; 2];
output = L.evaluate(input);

% Expected: [-1; -1; 1; 1; 1] (zero maps to +1)
expected = [-1; -1; 1; 1; 1];
assert(isequal(output, expected), 'SignLayer polar mode evaluate failed');

%% Test 5: SignLayer evaluate - nonnegative_zero_to_pos_one mode
L = SignLayer(0, 'nonnegative_zero_to_pos_one');

% Test with simple input
input = [-2; -1; 0; 1; 2];
output = L.evaluate(input);

% Expected: [0; 0; 1; 1; 1] (negative maps to 0, zero maps to +1)
expected = [0; 0; 1; 1; 1];
assert(isequal(output, expected), 'SignLayer nonnegative mode evaluate failed');

%% Test 6: SignLayer evaluate - 2D input
L = SignLayer();

% Create 2D input
input = [-3 -2 -1 0; 1 2 3 4];
output = L.evaluate(input);

% Sign should be applied element-wise
assert(all(output(1,:) == [-1 -1 -1 1]), 'First row should have correct signs');
assert(all(output(2,:) == [1 1 1 1]), 'Second row should have all +1');

%% Test 7: SignLayer sample
L = SignLayer();

% Create vertices
V = [-2 -1 0 1 2; -1 0 1 2 3];
Y = L.sample(V);

% Should apply sign to each vertex
assert(size(Y, 2) == size(V, 2), 'Should have same number of outputs as inputs');

%% Test 8: SignLayer reach with Star - exact-star
L = SignLayer();

% Create Star input
lb = [-1; -1];
ub = [1; 1];
B = Box(lb, ub);
I_star = B.toStar;

S = L.reach(I_star, 'exact-star');
assert(~isempty(S), 'SignLayer reach exact-star should return result');

%% Test 9: SignLayer reach with Star - approx-star
L = SignLayer();

% Create Star input
lb = [-1; -1];
ub = [1; 1];
B = Box(lb, ub);
I_star = B.toStar;

S = L.reach(I_star, 'approx-star');
assert(~isempty(S), 'SignLayer reach approx-star should return result');
