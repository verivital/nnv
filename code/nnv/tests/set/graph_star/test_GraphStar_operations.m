% test operations for GraphStar class

numNodes = 5;
numFeatures = 4;
NF = rand(numNodes, numFeatures);
LB = -0.1 * ones(numNodes, numFeatures);
UB = 0.1 * ones(numNodes, numFeatures);

GS = GraphStar(NF, LB, UB);

%% 1) toStar conversion
S = GS.toStar();
assert(isa(S, 'Star'), 'toStar should return a Star object');
assert(size(S.V, 1) == numNodes * numFeatures, 'Star dimension mismatch');
assert(S.nVar == GS.numPred, 'Number of predicates should match');

%% 2) affineMap
scale = 2;
offset = 1;
GS2 = GS.affineMap(scale, offset);
assert(isa(GS2, 'GraphStar'), 'affineMap should return a GraphStar');
assert(GS2.numNodes == numNodes, 'numNodes should be preserved');
assert(GS2.numFeatures == numFeatures, 'numFeatures should be preserved');

% Verify center is scaled and offset
expected_center = scale * GS.V(:,:,1) + offset;
assert(max(abs(GS2.V(:,:,1) - expected_center), [], 'all') < 1e-10, 'Center should be scaled and offset');

%% 3) MinkowskiSum
NF2 = rand(numNodes, numFeatures);
GS3 = GraphStar(NF2, LB, UB);
GS4 = GS.MinkowskiSum(GS3);
assert(isa(GS4, 'GraphStar'), 'MinkowskiSum should return a GraphStar');
assert(GS4.numNodes == numNodes, 'numNodes should be preserved');
assert(GS4.numFeatures == numFeatures, 'numFeatures should be preserved');
% Note: Star.MinkowskiSum may combine generators when constraints match,
% so we just verify numPred >= original (could be same or doubled)
assert(GS4.numPred >= GS.numPred, 'Predicates should be preserved or combined');

%% 4) estimateRanges
[lb_est, ub_est] = GS.estimateRanges();
assert(size(lb_est, 1) == numNodes, 'lb_est rows should match numNodes');
assert(size(lb_est, 2) == numFeatures, 'lb_est cols should match numFeatures');
assert(size(ub_est, 1) == numNodes, 'ub_est rows should match numNodes');
assert(size(ub_est, 2) == numFeatures, 'ub_est cols should match numFeatures');
assert(all(lb_est <= ub_est, 'all'), 'Lower bound should be <= upper bound');

%% 5) estimateRange (single element)
[xmin, xmax] = GS.estimateRange(1, 1);
assert(xmin <= xmax, 'xmin should be <= xmax');
assert(xmin >= lb_est(1,1) - 1e-10, 'xmin should match estimateRanges');
assert(xmax <= ub_est(1,1) + 1e-10, 'xmax should match estimateRanges');

%% 6) sample
samples = GS.sample(5);
assert(length(samples) == 5, 'Should return 5 samples');
assert(size(samples{1}, 1) == numNodes, 'Sample should have correct numNodes');
assert(size(samples{1}, 2) == numFeatures, 'Sample should have correct numFeatures');

%% 7) evaluate
pred_val = zeros(GS.numPred, 1);  % evaluate at center
nf_center = GS.evaluate(pred_val);
assert(size(nf_center, 1) == numNodes, 'Evaluated result should have correct numNodes');
assert(size(nf_center, 2) == numFeatures, 'Evaluated result should have correct numFeatures');
assert(max(abs(nf_center - GS.V(:,:,1)), [], 'all') < 1e-10, 'Zero predicates should give center');

%% 8) isEmptySet
res = GS.isEmptySet();
assert(res == 0, 'GraphStar should not be empty');

%% 9) contains
% The center should be contained in the set
center_nf = (GS.nf_lb + GS.nf_ub) / 2;
res_contains = GS.contains(center_nf);
assert(res_contains == 1, 'Center should be contained in the GraphStar');

% A point outside should not be contained
outside_nf = GS.nf_ub + 10;
res_outside = GS.contains(outside_nf);
assert(res_outside == 0, 'Point outside should not be contained');

disp('All GraphStar operations tests passed!');
