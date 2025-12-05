% test_soundness_SetTree
% Soundness tests for SetTree class
% SetTree is used for tracking reachable sets in NNCS
% To run: results = runtests('test_soundness_SetTree')

%% Test 1: SetTree constructor
rng(42);
height = 5;
tree = SetTree(height);
assert(tree.height == height, 'Height should be set correctly');
assert(length(tree.S) == height, 'S should have correct size');

%% Test 2: SetTree requires positive height
rng(42);
try
    tree = SetTree(0);
    assert(false, 'Should error on height = 0');
catch ME
    assert(contains(ME.message, '> 0'), 'Should mention positive height');
end

%% Test 3: Add reach set to tree
rng(42);
tree = SetTree(3);
lb = [-1; -1]; ub = [1; 1];
R = Star(lb, ub);
tree.addReachSet(R, 1);
extracted = tree.extractReachSet(1);
assert(~isempty(extracted), 'Should be able to extract added set');

%% Test 4: Add multiple reach sets
rng(42);
tree = SetTree(3);
R1 = Star([-1; -1], [1; 1]);
R2 = Star([-2; -2], [2; 2]);
tree.addReachSet(R1, 1);
tree.addReachSet(R2, 2);
assert(~isempty(tree.extractReachSet(1)), 'Position 1 should have set');
assert(~isempty(tree.extractReachSet(2)), 'Position 2 should have set');

%% Test 5: Position validation
rng(42);
tree = SetTree(3);
try
    tree.addReachSet(Star([-1;-1], [1;1]), 5);
    assert(false, 'Should error on position > height');
catch ME
    assert(true, 'Correctly rejects invalid position');
end

%% Test 6: Extract feedback reach set
rng(42);
tree = SetTree(3);
R1 = Star([-1; -1], [1; 1]);
tree.addReachSet(R1, 1);
fb = tree.extract_fb_ReachSet(1);
assert(~isempty(fb), 'Feedback set should be extractable');

%% Test 7: SetTree with array of sets
rng(42);
tree = SetTree(3);
R = [Star([-1; -1], [0; 0]), Star([0; 0], [1; 1])];
tree.addReachSet(R, 1);
extracted = tree.extractReachSet(1);
assert(length(extracted) == 2, 'Should store array of sets');

%% Test 8: Sequential additions build feedback tree
rng(42);
tree = SetTree(3);
R1 = Star([-1; -1], [1; 1]);
R2 = Star([-0.5; -0.5], [0.5; 0.5]);
tree.addReachSet(R1, 1);
tree.addReachSet(R2, 2);
fb1 = tree.extract_fb_ReachSet(1);
fb2 = tree.extract_fb_ReachSet(2);
assert(~isempty(fb1), 'Feedback at pos 1 should exist');
assert(~isempty(fb2), 'Feedback at pos 2 should exist');
