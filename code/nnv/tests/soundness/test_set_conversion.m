% test_set_conversion
% Tests for set representation conversions (Star, Zono, Box)
% To run: results = runtests('test_set_conversion')

%% Test 1: Star to Box conversion
rng(42);
lb = [-1; -2; -3];
ub = [1; 2; 3];
S = Star(lb, ub);
B = S.getBox();
assert(~isempty(B), 'Box should be created from Star');

%% Test 2: Star getRanges
rng(42);
lb = [-1; -2];
ub = [1; 2];
S = Star(lb, ub);
[lb_out, ub_out] = S.getRanges();
assert(all(lb_out == lb), 'Lower bounds should match');
assert(all(ub_out == ub), 'Upper bounds should match');

%% Test 3: Zono to Star conversion
rng(42);
c = [0; 0];
V = [1 0; 0 1];
Z = Zono(c, V);
S = Z.toStar();
assert(isa(S, 'Star'), 'Should convert to Star');

%% Test 4: Star affine map preserves containment
rng(42);
lb = [-1; -1];
ub = [1; 1];
S = Star(lb, ub);
W = [2 1; 1 2];
b = [0.5; -0.5];
S2 = S.affineMap(W, b);
x = [0.5; -0.3];
y = W * x + b;
assert(S2.contains(y), 'Transformed point should be in transformed set');

%% Test 5: Box contains original bounds
rng(42);
lb = [-2; -3; -1];
ub = [2; 3; 1];
B = Box(lb, ub);
center = (lb + ub) / 2;
assert(all(B.lb == lb), 'Box lower bound should match');
assert(all(B.ub == ub), 'Box upper bound should match');

%% Test 6: Zono generators define extent
rng(42);
c = [1; 2];
V = [1 0; 0 2];
Z = Zono(c, V);
[lb, ub] = Z.getBounds();
assert(lb(1) == 0, 'Lower bound X should be center - gen');
assert(ub(1) == 2, 'Upper bound X should be center + gen');

%% Test 7: Star intersection with HalfSpace
rng(42);
lb = [-1; -1];
ub = [1; 1];
S = Star(lb, ub);
G = [1 0];
g = 0.5;
S2 = S.intersectHalfSpace(G, g);
assert(~isempty(S2), 'Intersection should produce a set');

%% Test 8: Multiple set conversions roundtrip
rng(42);
lb = [-1; -1];
ub = [1; 1];
S1 = Star(lb, ub);
[lb1, ub1] = S1.getRanges();
B = Box(lb1, ub1);
S2 = B.toStar();
[lb2, ub2] = S2.getRanges();
tol = 1e-6;
assert(max(abs(lb1 - lb2)) < tol, 'Lower bounds should match after roundtrip');
assert(max(abs(ub1 - ub2)) < tol, 'Upper bounds should match after roundtrip');
