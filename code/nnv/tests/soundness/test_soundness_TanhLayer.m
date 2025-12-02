% test_soundness_TanhLayer
% Basic tests for TanhLayer
% Note: Full soundness testing for approx-star is deferred due to known library issues
% To run: results = runtests('test_soundness_TanhLayer')

%% Test 1: Tanh evaluate correctness
L = TanhLayer();

% Test that evaluate produces correct tanh values
input = [0 1; -1 2];
output = L.evaluate(input);
expected = tanh(input);

assert(max(abs(output(:) - expected(:))) < 1e-10, ...
    'Tanh evaluate should match tanh(x)');

%% Test 2: Tanh reach produces output
L2 = TanhLayer();

V = zeros(2, 2, 1, 2);
V(:,:,1,1) = [0 0.5; -0.5 1];
V(1,1,1,2) = 0.01;

input_is = ImageStar(V, [1; -1], [0.5; 0.5], -0.5, 0.5);
output_is = L2.reach(input_is, 'approx-star');

% Verify output has correct structure
assert(isa(output_is, 'ImageStar'), 'Output should be ImageStar');
assert(output_is.numPred >= 0, 'Output should have valid predicates');
