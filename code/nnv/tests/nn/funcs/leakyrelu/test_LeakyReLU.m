% Test LeakyReLU functionality
% To run: results = runtests('test_LeakyReLU')

%% Test 1: LeakyReLU evaluate
x = [-2; -1; 0; 1; 2];
gamma = 0.01;
y = LeakyReLU.evaluate(x, gamma);

expected = [-0.02; -0.01; 0; 1; 2];
assert(all(abs(y - expected) < 1e-6), 'LeakyReLU evaluate failed');

%% Test 2: LeakyReLU evaluate with different gamma
x = [-10; -5; 0; 5; 10];
gamma = 0.2;
y = LeakyReLU.evaluate(x, gamma);

expected = [-2; -1; 0; 5; 10];
assert(all(abs(y - expected) < 1e-6), 'LeakyReLU evaluate with gamma=0.2 failed');

%% Test 3: LeakyReLU reach exact star with crossing zero
% Create a star set that crosses zero
lb = [-1; -1];
ub = [1; 1];

I = Star(lb, ub);
gamma = 0.01;

S = LeakyReLU.reach(I, gamma, 'exact-star');
assert(~isempty(S), 'LeakyReLU reach exact-star should return result');

%% Test 4: LeakyReLU reach star approx
% Create input star
lb = [-1; -1];
ub = [1; 1];
I_star = Star(lb,ub);

gamma = 0.1;
S = LeakyReLU.reach(I_star, gamma, 'approx-star');
assert(~isempty(S), 'LeakyReLU reach approx-star should return result');

%% Test 5: LeakyReLU reach zono approx
lb = [-1; -1];
ub = [1; 1];
B = Box(lb, ub);
I_zono = B.toZono;

gamma = 0.01;
Z = LeakyReLU.reach(I_zono, gamma, 'approx-zono');
assert(isa(Z, 'Zono'), 'LeakyReLU reach approx-zono should return Zono');

%% Test 6: LeakyReLU stepReach with positive range
lb = [0.5;1];
ub = [2.5;2];
I = Star(lb, ub);

gamma = 0.01;
index = 1;
S = LeakyReLU.stepReach(I, index, gamma);

% Should return unchanged for positive values
assert(isa(S, 'Star'), 'stepReach should return Star');
assert(isequal(S.V, I.V), 'stepReach with positive range should not change V');

%% Test 7: LeakyReLU stepReach with negative range
lb = [-1;-10];
ub = [-0.1; -8];
I = Star(lb, ub);

gamma = 0.1;
index = 1;
S = LeakyReLU.stepReach(I, index, gamma);

% Should scale by gamma for negative values
assert(isa(S, 'Star'), 'stepReach should return Star');
assert(S.V(index, 1) == gamma * I.V(index, 1), 'stepReach should scale negative values');

%% Test 8: LeakyReLU stepReach with range crossing zero
V = [0 1; 0 1];
C = [1; -1];
d = [1; 1];
lb = [-1];
ub = [1];
I = Star(V, C, d, lb, ub);

gamma = 0.01;
index = 1;
S = LeakyReLU.stepReach(I, index, gamma);

% Should return multiple stars for crossing zero
assert(length(S) == 2, 'stepReach with crossing zero should return 2 stars');
