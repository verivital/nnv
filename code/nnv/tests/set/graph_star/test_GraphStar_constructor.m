% test constructor for GraphStar class

%% Test case 3: center + perturbations
numNodes = 5;
numFeatures = 4;

NF = rand(numNodes, numFeatures);  % center node features
LB = -0.1 * ones(numNodes, numFeatures);  % lower perturbation
UB = 0.1 * ones(numNodes, numFeatures);   % upper perturbation

GS = GraphStar(NF, LB, UB);

assert(GS.numNodes == numNodes, 'numNodes mismatch');
assert(GS.numFeatures == numFeatures, 'numFeatures mismatch');
assert(GS.numPred > 0, 'Should have predicates for perturbations');
assert(size(GS.V, 1) == numNodes, 'V dimension 1 mismatch');
assert(size(GS.V, 2) == numFeatures, 'V dimension 2 mismatch');
assert(size(GS.V, 3) == GS.numPred + 1, 'V dimension 3 mismatch');

%% Test case 2: from bounds
lb = NF + LB;
ub = NF + UB;
GS2 = GraphStar(lb, ub);

assert(GS2.numNodes == numNodes, 'numNodes mismatch (case 2)');
assert(GS2.numFeatures == numFeatures, 'numFeatures mismatch (case 2)');
assert(GS2.numPred > 0, 'Should have predicates (case 2)');
assert(isequal(GS2.nf_lb, lb), 'nf_lb should match input lb');
assert(isequal(GS2.nf_ub, ub), 'nf_ub should match input ub');

%% Test case 5: Star representation
V = rand(numNodes, numFeatures, 3);  % 2 predicates
C = [1 0; 0 1; -1 0; 0 -1];
d = [1; 1; 1; 1];
pred_lb = [-1; -1];
pred_ub = [1; 1];

GS5 = GraphStar(V, C, d, pred_lb, pred_ub);

assert(GS5.numNodes == numNodes, 'numNodes mismatch (case 5)');
assert(GS5.numFeatures == numFeatures, 'numFeatures mismatch (case 5)');
assert(GS5.numPred == 2, 'numPred should be 2 (case 5)');

%% Test case 7: Full specification with cached bounds
nf_lb = rand(numNodes, numFeatures);
nf_ub = nf_lb + 0.1;

GS7 = GraphStar(V, C, d, pred_lb, pred_ub, nf_lb, nf_ub);

assert(GS7.numNodes == numNodes, 'numNodes mismatch (case 7)');
assert(isequal(GS7.nf_lb, nf_lb), 'nf_lb should match (case 7)');
assert(isequal(GS7.nf_ub, nf_ub), 'nf_ub should match (case 7)');

%% Test case 0: Empty GraphStar
GS0 = GraphStar();

assert(GS0.numNodes == 0, 'numNodes should be 0 for empty');
assert(GS0.numFeatures == 0, 'numFeatures should be 0 for empty');
assert(isempty(GS0.V), 'V should be empty');

disp('All GraphStar constructor tests passed!');
