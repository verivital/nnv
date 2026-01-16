% test_GraphStar_constructor.m - Unit tests for GraphStar constructors
%
% Tests various constructor forms: empty, bounds, center+perturbations, Star representation

% Shared setup (before any %% sections)
numNodes = 5;
numFeatures = 4;
NF = rand(numNodes, numFeatures);
LB = -0.1 * ones(numNodes, numFeatures);
UB = 0.1 * ones(numNodes, numFeatures);

%% 1) Empty constructor
GS0 = GraphStar();
assert(GS0.numNodes == 0, 'numNodes should be 0 for empty');
assert(GS0.numFeatures == 0, 'numFeatures should be 0 for empty');
assert(isempty(GS0.V), 'V should be empty');

%% 2) Constructor from bounds (lb, ub)
lb = NF + LB;
ub = NF + UB;
GS2 = GraphStar(lb, ub);

assert(GS2.numNodes == numNodes, 'numNodes mismatch (from bounds)');
assert(GS2.numFeatures == numFeatures, 'numFeatures mismatch (from bounds)');
assert(GS2.numPred > 0, 'Should have predicates (from bounds)');
assert(isequal(GS2.nf_lb, lb), 'nf_lb should match input lb');
assert(isequal(GS2.nf_ub, ub), 'nf_ub should match input ub');

%% 3) Constructor from center + perturbations (NF, LB, UB)
GS3 = GraphStar(NF, LB, UB);

assert(GS3.numNodes == numNodes, 'numNodes mismatch (center+perturbations)');
assert(GS3.numFeatures == numFeatures, 'numFeatures mismatch (center+perturbations)');
assert(GS3.numPred > 0, 'Should have predicates for perturbations');
assert(size(GS3.V, 1) == numNodes, 'V dimension 1 mismatch');
assert(size(GS3.V, 2) == numFeatures, 'V dimension 2 mismatch');
assert(size(GS3.V, 3) == GS3.numPred + 1, 'V dimension 3 mismatch');

%% 4) Constructor from Star representation (V, C, d, pred_lb, pred_ub)
V = rand(numNodes, numFeatures, 3);  % 2 predicates
C = [1 0; 0 1; -1 0; 0 -1];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

GS4 = GraphStar(V, C, d, pred_lb, pred_ub);

assert(GS4.numNodes == numNodes, 'numNodes mismatch (Star representation)');
assert(GS4.numFeatures == numFeatures, 'numFeatures mismatch (Star representation)');
assert(GS4.numPred == 2, 'numPred should be 2 (Star representation)');

%% 5) Full constructor with cached bounds (V, C, d, pred_lb, pred_ub, nf_lb, nf_ub)
V5 = rand(numNodes, numFeatures, 3);  % 2 predicates
C5 = [1 0; 0 1; -1 0; 0 -1];
d5 = [1; 1; 1; 1];
pred_lb5 = [-1; -1];
pred_ub5 = [1; 1];
nf_lb = rand(numNodes, numFeatures);
nf_ub = nf_lb + 0.1;

GS5 = GraphStar(V5, C5, d5, pred_lb5, pred_ub5, nf_lb, nf_ub);

assert(GS5.numNodes == numNodes, 'numNodes mismatch (full specification)');
assert(isequal(GS5.nf_lb, nf_lb), 'nf_lb should match (full specification)');
assert(isequal(GS5.nf_ub, nf_ub), 'nf_ub should match (full specification)');

disp('All GraphStar constructor tests passed!');
