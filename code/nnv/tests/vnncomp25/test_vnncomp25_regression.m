%% test_vnncomp25_regression
% Regression tests for the V04 fixes that landed in NNV's master tree
% during the VNN-COMP 2025 importer work. Each test reproduces the bug
% the fix addresses, so future changes can't silently regress.
%
% Bugs covered:
%  (1) NN.evaluate_withConns / reach_withConns DAG `y`/`outSet` not
%      refreshed when source layer was already evaluated by an earlier
%      iteration. Symptom: residual networks gave wrong outputs.
%  (2) FlattenLayer.reach_multipleInputs rejected Star inputs (only
%      handled ImageStar/ImageZono/VolumeStar). Symptom: FC-only nets
%      with a leading Flatten errored on reach.
%  (3) AdditionLayer.reach_single_input crashed when handed a single
%      Star (non-cell). Symptom: DAG residual reach error.
%  (4) load_nnv_from_mat ConcatenationLayer constructor signature
%      (5 vs 6 args).
%
% Run:  results = runtests('test_vnncomp25_regression')

%% Test 1: NN.evaluate_withConns refreshes y when source already evaluated
% Build a tiny DAG: input -> fc1 -> relu -> fc2 -> add(fc2, fc_short).
% input is also fed directly to fc_short (residual). With the bug, the
% second iteration reading from "input" would inherit `y` from the
% previous iteration (relu's output) and feed it to fc_short instead of
% the actual input.

rng(0);
in_size = 4;

W1 = randn(8, in_size, 'single');  b1 = zeros(8,1,'single');
W2 = randn(4, 8,       'single');  b2 = zeros(4,1,'single');
Ws = randn(4, in_size, 'single');  bs = zeros(4,1,'single');

L1 = FeatureInputLayer('input', in_size, 'none', 'auto', [], [], [], []);
L2 = FullyConnectedLayer('fc1',   double(W1),  double(b1));
L3 = ReluLayer('relu1');
L4 = FullyConnectedLayer('fc2',   double(W2),  double(b2));
L5 = FullyConnectedLayer('fc_short', double(Ws), double(bs));
L6 = AdditionLayer('add', 2, 1, {'in1','in2'}, {'out'});

Layers = {L1, L2, L3, L4, L5, L6};

Sources      = ["input"; "fc1";  "relu1"; "input";    "fc2";    "fc_short"];
Destinations = ["fc1";   "relu1";"fc2";   "fc_short"; "add/in1";"add/in2"];
Connections  = table(Sources, Destinations, 'VariableNames', {'Source','Destination'});

net = NN(Layers, Connections);
names = cellfun(@(L) L.Name, Layers, 'UniformOutput', false);
net.name2indx = containers.Map(names, 1:numel(Layers));

x = single(randn(in_size, 1));

% Reference: hand-compute
ref = double(W2) * max(double(W1)*double(x) + double(b1), 0) + double(b2);
ref = ref + double(Ws)*double(x) + double(bs);

y = double(net.evaluate(x));
assert(max(abs(y(:) - ref(:))) < 1e-5, ...
    'Test 1 failed (DAG y-refresh): NN.evaluate output does not match hand-computed residual');
fprintf('Test 1 PASSED (DAG y-refresh): max diff %.4e\n', max(abs(y(:) - ref(:))));


%% Test 2: FlattenLayer.reach passes Star through (V04 fix)
L = FlattenLayer('flat');
L.Type = 'nnet.onnx.layer.FlattenInto2dLayer';
S = Star(zeros(8,1), ones(8,1));
out = L.reach(S, 'approx-star', 'single');
assert(isa(out, 'Star'), 'Test 2 failed: expected Star out');
assert(out.dim == S.dim, 'Test 2 failed: dim changed');
fprintf('Test 2 PASSED (FlattenLayer Star passthrough)\n');


%% Test 3: AdditionLayer.reach_single_input handles non-cell Star
L = AdditionLayer('add', 2, 1, {'in1','in2'}, {'out'});
S = Star(zeros(4,1), ones(4,1));
out = L.reach_single_input(S);  % non-cell input — should not crash
assert(isa(out, 'Star'), 'Test 3 failed: expected Star out');
fprintf('Test 3 PASSED (AdditionLayer single Star non-cell)\n');


%% Test 4: ConcatenationLayer 6-arg constructor (axis pass-through)
% Just construct it; no network run needed.
L = ConcatenationLayer('cc', 2, 1, {'in1','in2'}, {'out'}, 1);
assert(strcmp(L.Name, 'cc'));
fprintf('Test 4 PASSED (ConcatenationLayer 6-arg)\n');


%% Test 5: TransposedConv2DLayer accepts [1,1,F] bias
W = randn(3,3,8,4,'single');   % [kH, kW, NumFilters=8, NumChannels=4]
b = reshape(zeros(8,1,'single'), 1, 1, 8);
pads = [0 0 0 0]; strides = [1 1];
L = TransposedConv2DLayer('ct', double(W), double(b), pads, strides);
assert(strcmp(L.Name, 'ct'));
assert(L.NumFilters == 8 && L.NumChannels == 4);
fprintf('Test 5 PASSED (TransposedConv2DLayer 3D bias)\n');


%% Test 6: TransposedConv2DLayer accepts NumFilters=1 with [1,1] bias
% MATLAB strips trailing singletons so reshape([0],1,1,1) returns [1,1].
% The constructor must not crash on this case (case-2 / 4 / 6 paths).
W = randn(3,3,1,4,'single');   % NumF=1, NumC=4
b1 = reshape(zeros(1,1,'single'), 1, 1, 1);   % becomes [1,1]
% case 2 (no pad/stride)
L = TransposedConv2DLayer(double(W), double(b1));
assert(L.NumFilters == 1 && L.NumChannels == 4);
% case 4 (no name)
L = TransposedConv2DLayer(double(W), double(b1), [0 0 0 0], [1 1]);
assert(L.NumFilters == 1 && L.NumChannels == 4);
% case 6 (with name)
L = TransposedConv2DLayer('ct1', double(W), double(b1), [0 0 0 0], [1 1]);
assert(L.NumFilters == 1 && L.NumChannels == 4);
fprintf('Test 6 PASSED (TransposedConv2DLayer single-filter trailing-singleton bias)\n');


%% Test 7: ReshapeLayer ONNX-BCHW semantics (flat -> [H,W,C])
% Fresh-from-FC flat vector, ONNX would reshape to [B,C,H,W] in C-order.
% NNV's HWC tensor must equal that BCHW tensor for matching downstream.
H = 2; W = 3; C = 4;
n = H*W*C;
flat = reshape(single(0:n-1), n, 1);
% Expected HWC built from ONNX BCHW interpretation:
% BCHW[c,h,w] = flat[c*H*W + h*W + w]; HWC[h,w,c] = BCHW[c,h,w].
expected = zeros(H, W, C, 'single');
for c = 1:C
    for h = 1:H
        for w = 1:W
            expected(h,w,c) = flat((c-1)*H*W + (h-1)*W + (w-1) + 1);
        end
    end
end
L = ReshapeLayer('rs', [H W C]);
L.OnnxBCHW = true;
got = L.evaluate(flat);
assert(isequal(size(got), [H W C]), 'Test 7 failed: shape');
assert(max(abs(got(:) - expected(:))) < 1e-6, ...
    'Test 7 failed: BCHW->HWC reshape did not match ONNX semantics');
fprintf('Test 7 PASSED (ReshapeLayer OnnxBCHW BCHW->HWC)\n');


%% Test 8: BatchNormalizationLayer accepts [1,1,C] params on HWC input
% load_nnv_from_mat now reshapes BN params to [1,1,C] so the rank-3 HWC
% branch of BN.evaluate works. Build a layer and run a forward pass.
nC = 4;
mean_ = reshape(single(zeros(nC,1)), 1, 1, nC);
var_  = reshape(single(ones(nC,1)),  1, 1, nC);
off   = reshape(single(zeros(nC,1)), 1, 1, nC);
sca   = reshape(single(ones(nC,1)),  1, 1, nC);
L = BatchNormalizationLayer('Name','bn','NumChannels',nC, ...
    'TrainedMean', mean_, 'TrainedVariance', var_, ...
    'Offset', off, 'Scale', sca, 'Epsilon', 1e-5);
x = single(randn(2, 3, nC));
y = L.evaluate(x);
assert(isequal(size(y), [2 3 nC]), 'Test 8 failed: BN HWC shape mismatch');
% With mean=0, var=1, scale=1, off=0, BN is approximately identity.
assert(max(abs(y(:) - x(:))) < 1e-3, 'Test 8 failed: BN identity within eps');
fprintf('Test 8 PASSED (BatchNormalizationLayer [1,1,C] params + HWC input)\n');


%% Test 9: LeakyReluLayer constructor accepts 6 args (loader call shape)
% load_nnv_from_mat constructs LeakyReluLayer via the 6-arg path; previously
% was 2-arg which crashed at varargin{6}.
L = LeakyReluLayer('lrelu1', 1, {'in'}, 1, {'out'}, 0.1);
assert(strcmp(L.Name, 'lrelu1') && abs(L.gamma - 0.1) < 1e-9);
fprintf('Test 9 PASSED (LeakyReluLayer 6-arg ctor)\n');


%% Test 10: ConcatenationLayer variadic NumInputs (was hard-coded to 2)
% nn4sys uses Concat with 6 inputs.
L = ConcatenationLayer('cc6', 6, 1, {'in1','in2','in3','in4','in5','in6'}, {'out'}, 1);
assert(L.NumInputs == 6);
ins = cell(6,1);
for k=1:6, ins{k} = single([k; k+0.5]); end
y = L.evaluate(ins);
assert(isequal(size(y), [12 1]), 'Test 10 failed: variadic Concat shape');
expected = [1; 1.5; 2; 2.5; 3; 3.5; 4; 4.5; 5; 5.5; 6; 6.5];
assert(max(abs(double(y) - expected)) < 1e-6, 'Test 10 failed: cat values');
fprintf('Test 10 PASSED (ConcatenationLayer variadic NumInputs)\n');


%% Test 11: PlaceholderLayer with Perm applies a real permutation
L = PlaceholderLayer('tr', 'Transpose');
L.Perm = [2 1];   % swap dims 1,2
x = single([1 2 3; 4 5 6]);   % 2x3
y = L.evaluate(x);
assert(isequal(size(y), [3 2]));
assert(isequal(double(y), double(x.')));
fprintf('Test 11 PASSED (PlaceholderLayer Perm permute)\n');


%% Test 12: PlaceholderLayer Sign and Abs apply real element-wise ops
Ls = PlaceholderLayer('s', 'Sign');
xs = single([-2 -0.5 0 0.5 2]);
ys = Ls.evaluate(xs);
assert(isequal(double(ys), [-1 -1 0 1 1]), 'Test 12 failed: Sign');
La = PlaceholderLayer('a', 'Abs');
ya = La.evaluate(xs);
assert(isequal(double(ya), [2 0.5 0 0.5 2]), 'Test 12 failed: Abs');
fprintf('Test 12 PASSED (PlaceholderLayer Sign + Abs)\n');


%% Test 13: SignLayer (NNV's polar_zero_to_pos_one mode)
Ls13 = SignLayer(0, 'polar_zero_to_pos_one');
ys13 = Ls13.evaluate(single([-2 -0.5 0 0.5 2]));
assert(isequal(double(ys13), [-1 -1 1 1 1]), 'Test 13 failed: SignLayer eval');
fprintf('Test 13 PASSED (SignLayer wired with polar_zero_to_pos_one)\n');


%% Test 14: ElementwiseAffineLayer with multi-dim bias preserves shape
% VGG-style per-pixel mean Sub: bias is [H,W,C] not flat.
H = 4; W = 4; C = 3;
b3d = ones(H,W,C, 'double') * 0.5;
s3d = ones(1,1,1, 'double');  % no scale
L14 = ElementwiseAffineLayer('eaff_3d', s3d, -b3d, false, true);  % subtract
x14 = ones(H, W, C, 'double');
y14 = L14.evaluate(x14);
assert(isequal(size(y14), [H W C]), 'Test 14 failed: shape');
assert(max(abs(y14(:) - 0.5)) < 1e-9, 'Test 14 failed: subtraction values');
fprintf('Test 14 PASSED (ElementwiseAffineLayer 3D HWC bias)\n');


%% Test 15: SoftmaxLayer reach soundness (intermediate softmax must produce
%% bounds in [0, 1], identity-only behavior is unsound for non-final softmax).
L15 = SoftmaxLayer('sm15');
L15.IsFinalLayer = false;
S15 = Star([0; 0; 0], [1; 1; 1]);
Sout15 = L15.reach(S15, 'approx-star');
if ~isempty(Sout15) && isa(Sout15, 'Star')
    [lb15, ub15] = Sout15.getRanges();
    assert(all(lb15 >= -1e-3), 'Test 15 failed: lb < 0');
    assert(all(ub15 <= 1+1e-3), 'Test 15 failed: ub > 1');
    fprintf('Test 15 PASSED (SoftmaxLayer non-final reach in [0,1])\n');
else
    error('Test 15 failed: SoftmaxLayer non-final returned empty/wrong type');
end


%% Test 16: SoftmaxLayer final-layer reach is identity (legacy preservation)
L16 = SoftmaxLayer('sm16');   % default IsFinalLayer = true
S16 = Star([1; 2; 3], [4; 5; 6]);
Sout16 = L16.reach(S16, 'approx-star');
[lb16a, ub16a] = S16.getRanges();
[lb16b, ub16b] = Sout16.getRanges();
assert(max(abs(lb16a - lb16b)) < 1e-9 && max(abs(ub16a - ub16b)) < 1e-9, ...
    'Test 16 failed: final-layer softmax should be identity');
fprintf('Test 16 PASSED (SoftmaxLayer final-layer identity preserved)\n');


fprintf('\n=== All V04/V05/V06/V07 NNV regression tests PASSED ===\n');
