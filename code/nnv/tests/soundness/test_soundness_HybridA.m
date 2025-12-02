% test_soundness_HybridA
% Soundness tests for Hybrid Automata class
% Tests HybridA and HybridANNCS classes (CORA wrapper)
% To run: results = runtests('test_soundness_HybridA')

%% Test 1: HybridA class exists
rng(42);
assert(exist('HybridA', 'class') == 8, 'HybridA class should exist');

%% Test 2: HybridANNCS class exists
rng(42);
assert(exist('HybridANNCS', 'class') == 8, 'HybridANNCS class should exist');

%% Test 3: HybridA has stepReachStar method
rng(42);
mc = ?HybridA;
methods = {mc.MethodList.Name};
assert(ismember('stepReachStar', methods), 'HybridA should have stepReachStar');

%% Test 4: HybridA has reach_zono method
rng(42);
mc = ?HybridA;
methods = {mc.MethodList.Name};
assert(ismember('reach_zono', methods), 'HybridA should have reach_zono');

%% Test 5: HybridANNCS has reach method
rng(42);
mc = ?HybridANNCS;
methods = {mc.MethodList.Name};
assert(ismember('reach', methods), 'HybridANNCS should have reach method');

%% Test 6: HybridANNCS has verify method
rng(42);
mc = ?HybridANNCS;
methods = {mc.MethodList.Name};
assert(ismember('verify', methods), 'HybridANNCS should have verify method');

%% Test 7: CORA installed check
rng(42);
cora_path = [nnvroot(), filesep, 'code', filesep, 'nnv', filesep, 'engine', filesep, 'cora'];
assert(exist(cora_path, 'dir') == 7, 'CORA directory should exist');

%% Test 8: HybridA properties
rng(42);
mc = ?HybridA;
props = {mc.PropertyList.Name};
assert(ismember('dim', props), 'HybridA should have dim property');
assert(ismember('modes', props), 'HybridA should have modes property');
assert(ismember('controlPeriod', props), 'HybridA should have controlPeriod');
