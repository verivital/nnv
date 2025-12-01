% test_lpsolver
% Unit tests for NNV LP solver interface
% Tests lpsolver function with various LP problems
% To run: results = runtests('test_lpsolver')

%% Test 1: lpsolver function exists
rng(42);
assert(exist('lpsolver', 'file') == 2, 'lpsolver function should exist');

%% Test 2: Simple bounded LP with linprog
rng(42);
% min x1 + x2 subject to x1 >= 0, x2 >= 0, x1 + x2 <= 1
f = [1; 1];  % objective: minimize x1 + x2
A = [1, 1];  % constraint: x1 + x2 <= 1
b = 1;
lb = [0; 0]; % lower bounds
ub = [1; 1]; % upper bounds
[fval, exitflag] = lpsolver(f, A, b, [], [], lb, ub, 'linprog');
assert(abs(fval) < 1e-6, 'Optimal value should be 0');
assert(startsWith(string(exitflag), 'l') || startsWith(string(exitflag), 'g'), 'Exit flag should indicate success');

%% Test 3: LP with unique solution
rng(42);
% min -x1 - x2 subject to x1 + x2 <= 1, x1 >= 0, x2 >= 0
f = [-1; -1];  % maximize x1 + x2
A = [1, 1];    % x1 + x2 <= 1
b = 1;
lb = [0; 0];
ub = [inf; inf];
[fval, exitflag] = lpsolver(f, A, b, [], [], lb, ub, 'linprog');
assert(abs(fval - (-1)) < 1e-6, 'Optimal value should be -1');

%% Test 4: LP with equality constraints
rng(42);
% min x1 subject to x1 + x2 = 2, x1 >= 0, x2 >= 0
f = [1; 0];    % minimize x1
A = [];        % no inequality constraints
b = [];
Aeq = [1, 1];  % x1 + x2 = 2
Beq = 2;
lb = [0; 0];
ub = [inf; inf];
[fval, exitflag] = lpsolver(f, A, b, Aeq, Beq, lb, ub, 'linprog');
assert(abs(fval) < 1e-6, 'Optimal value should be 0 (x1=0, x2=2)');

%% Test 5: Handles single precision input
rng(42);
f = single([1; 1]);
A = single([1, 1]);
b = single(1);
lb = single([0; 0]);
ub = single([1; 1]);
[fval, exitflag] = lpsolver(f, A, b, [], [], lb, ub, 'linprog');
assert(abs(fval) < 1e-6, 'Should handle single precision');

%% Test 6: emptySet option for feasibility check
rng(42);
% Check if region x1 >= 0, x2 >= 0, x1 + x2 <= 1 is non-empty
f = [0; 0];    % constant objective
A = [1, 1];
b = 1;
lb = [0; 0];
ub = [inf; inf];
[~, exitflag] = lpsolver(f, A, b, [], [], lb, ub, 'linprog', 'emptySet');
assert(startsWith(string(exitflag), 'l') || startsWith(string(exitflag), 'g'), 'Should find feasible region');

%% Test 7: GLPK solver (if available)
rng(42);
f = [1; 1];
A = [1, 1];
b = 1;
lb = [0; 0];
ub = [1; 1];
try
    [fval, exitflag] = lpsolver(f, A, b, [], [], lb, ub, 'glpk');
    assert(abs(fval) < 1e-6, 'GLPK should find optimal value');
catch ME
    % GLPK might not be available, which is acceptable
    assert(contains(ME.message, 'glpk') || contains(ME.message, 'Undefined'), ...
        'Error should be about GLPK availability');
end

%% Test 8: Return values are valid types
rng(42);
f = [1; 1];
A = [1, 1];
b = 1;
lb = [0; 0];
ub = [1; 1];
[fval, exitflag] = lpsolver(f, A, b, [], [], lb, ub, 'linprog');
assert(isnumeric(fval), 'fval should be numeric');
assert(isscalar(fval), 'fval should be scalar');
assert(isstring(exitflag) || ischar(exitflag), 'exitflag should be string');
