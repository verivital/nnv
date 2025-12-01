% test_soundness_NonlinearNNCS
% Soundness tests for NNCS-related classes
% To run: results = runtests('test_soundness_NonlinearNNCS')

%% Test 1: NonlinearNNCS class exists
rng(42);
assert(exist('NonlinearNNCS', 'class') == 8, 'NonlinearNNCS class should exist');

%% Test 2: LinearODE class exists
rng(42);
assert(exist('LinearODE', 'class') == 8, 'LinearODE class should exist');

%% Test 3: LinearODE plant creation
rng(42);
A = [-1 0; 0 -2];
B = [1; 0];
C = [1 0];
D = 0;
plant = LinearODE(A, B, C, D);
assert(~isempty(plant), 'LinearODE should be created');
assert(plant.dim == 2, 'Plant should have 2 states');

%% Test 4: HybridA class exists
rng(42);
assert(exist('HybridA', 'class') == 8, 'HybridA class should exist');

%% Test 5: DLinearODE class exists
rng(42);
assert(exist('DLinearODE', 'class') == 8, 'DLinearODE class should exist');

%% Test 6: NonLinearODE class exists
rng(42);
assert(exist('NonLinearODE', 'class') == 8, 'NonLinearODE class should exist');

%% Test 7: DNonLinearODE class exists
rng(42);
assert(exist('DNonLinearODE', 'class') == 8, 'DNonLinearODE class should exist');

%% Test 8: NNCS class exists
rng(42);
assert(exist('NNCS', 'class') == 8, 'NNCS class should exist');
