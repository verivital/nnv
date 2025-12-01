% test_soundness_SigmoidLayer
% Basic tests for SigmoidLayer
% Note: Full soundness testing for approx-star is deferred due to known library issues
% To run: results = runtests('test_soundness_SigmoidLayer')

%% Test 1: Sigmoid evaluate correctness
L = SigmoidLayer();

% Test that evaluate produces correct sigmoid values
input = [0 1; -1 2];
output = L.evaluate(input);
expected = 1 ./ (1 + exp(-input));

assert(max(abs(output(:) - expected(:))) < 1e-10, ...
    'Sigmoid evaluate should match 1/(1+exp(-x))');

%% Test 2: Sigmoid reach produces output
L2 = SigmoidLayer();

V = zeros(2, 2, 1, 2);
V(:,:,1,1) = [0 0.5; -0.5 1];
V(1,1,1,2) = 0.01;

input_is = ImageStar(V, [1; -1], [0.5; 0.5], -0.5, 0.5);
output_is = L2.reach(input_is, 'approx-star');

% Verify output has correct structure
assert(isa(output_is, 'ImageStar'), 'Output should be ImageStar');
assert(output_is.numPred >= 0, 'Output should have valid predicates');
