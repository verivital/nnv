% test_soundness_DNonlinearNNCS
% Soundness tests for Discrete Nonlinear NNCS class
% Tests DNonlinearNNCS and related plant classes
% To run: results = runtests('test_soundness_DNonlinearNNCS')

%% Test 1: DNonlinearNNCS class exists
rng(42);
assert(exist('DNonlinearNNCS', 'class') == 8, 'DNonlinearNNCS class should exist');

%% Test 2: DLinearODE class exists
rng(42);
assert(exist('DLinearODE', 'class') == 8, 'DLinearODE class should exist');

%% Test 3: DNonLinearODE class exists
rng(42);
assert(exist('DNonLinearODE', 'class') == 8, 'DNonLinearODE class should exist');

%% Test 4: DLinearNNCS class exists
rng(42);
assert(exist('DLinearNNCS', 'class') == 8, 'DLinearNNCS class should exist');

%% Test 5: Create simple DLinearODE plant
rng(42);
A = [0.9 0.1; 0 0.95];
B = [0; 0.1];
C = eye(2);
D = zeros(2, 1);
Ts = 0.1;
plant = DLinearODE(A, B, C, D, Ts);
assert(~isempty(plant), 'DLinearODE plant should be created');

%% Test 6: DLinearODE has reach method
rng(42);
mc = ?DLinearODE;
methods = {mc.MethodList.Name};
assert(ismember('stepReachStar', methods) || ismember('reach', methods), ...
    'DLinearODE should have reachability method');

%% Test 7: DNonlinearNNCS has reach method
rng(42);
mc = ?DNonlinearNNCS;
methods = {mc.MethodList.Name};
assert(ismember('reach', methods), 'DNonlinearNNCS should have reach method');

%% Test 8: DNonlinearNNCS has verify method
rng(42);
mc = ?DNonlinearNNCS;
methods = {mc.MethodList.Name};
assert(ismember('verify', methods), 'DNonlinearNNCS should have verify method');
