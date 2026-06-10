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


%% Test 17: ElementwiseProductLayer (V08, dynamic Mul)
L17 = ElementwiseProductLayer('mul17', 2, 1, {'in1','in2'}, {'out'});
ins17 = {single([1; 2; 3]), single([4; 5; 6])};
y17 = L17.evaluate(ins17);
assert(isequal(double(y17), [4; 10; 18]), 'Test 17 failed: Hadamard product');
% Reach (sound interval over-approx)
Sa = Star([0; 0], [1; 1]);   % [0, 1]^2
Sb = Star([1; 1], [2; 2]);   % [1, 2]^2
Sout17 = L17.reach({Sa, Sb}, 'approx-star');
[lb17, ub17] = Sout17.getRanges();
% [0,1] * [1,2] elementwise product should be [0, 2] each
assert(all(abs(lb17 - 0) < 1e-9) && all(abs(ub17 - 2) < 1e-9), ...
    'Test 17 failed: bilinear product bounds wrong');
fprintf('Test 17 PASSED (ElementwiseProductLayer Hadamard + interval reach)\n');


%% Test 18: ElementwiseDivisionLayer (V08, dynamic Div)
L18 = ElementwiseDivisionLayer('div18', 2, 1, {'in1','in2'}, {'out'});
ins18 = {single([4; 10; 18]), single([2; 5; 6])};
y18 = L18.evaluate(ins18);
assert(isequal(double(y18), [2; 2; 3]), 'Test 18 failed: Hadamard division');
% Reach (divisor strictly positive — well-defined)
Sn = Star([4; 6], [8; 12]);   % numerator in [4, 8] x [6, 12]
Sd = Star([2; 2], [4; 4]);    % divisor in [2, 4] x [2, 4]
Sout18 = L18.reach({Sn, Sd}, 'approx-star');
[lb18, ub18] = Sout18.getRanges();
% [4,8]/[2,4] = [1, 4]; [6,12]/[2,4] = [1.5, 6]
assert(all(abs(lb18 - [1; 1.5]) < 1e-9), 'Test 18 failed: division lb');
assert(all(abs(ub18 - [4; 6]) < 1e-9), 'Test 18 failed: division ub');
fprintf('Test 18 PASSED (ElementwiseDivisionLayer + interval reach)\n');


%% Test 19: DynamicMatmulLayer (V09, dynamic Q*K^T attention)
L19 = DynamicMatmulLayer('mm19', 2, 1, {'in1','in2'}, {'out'});
A = rand(3, 4, 'single');
B = rand(4, 5, 'single');
y19 = L19.evaluate({A, B});
expected = A * B;
assert(isequal(size(y19), [3 5]), 'Test 19 failed: 2D matmul shape');
assert(max(abs(double(y19(:)) - double(expected(:)))) < 1e-5, 'Test 19 failed: matmul values');
fprintf('Test 19 PASSED (DynamicMatmulLayer 2D matmul)\n');


%% Test 20: ReshapeLayer -1 sentinel resolves correctly (V09)
L20 = ReshapeLayer('rs20', [1 48 -1]);
x20 = ones(2, 2, 48, 'single');   % numel = 192
y20 = L20.evaluate(x20);
assert(isequal(size(y20), [1 48 4]), 'Test 20 failed: -1 should resolve to 4');
fprintf('Test 20 PASSED (ReshapeLayer -1 sentinel)\n');


%% Test 21: ElementwiseAffineLayer align_to_input (V09)
% Bias [1,1,1,5] aligned to input [5,1] should reshape to [5,1] not broadcast 4D.
L21 = ElementwiseAffineLayer('ea21', single(1), reshape(single(1:5), 1,1,1,5), false, true);
x21 = single([10; 20; 30; 40; 50]);   % column [5,1]
y21 = L21.evaluate(x21);
assert(isequal(size(y21), [5 1]), 'Test 21 failed: shape mismatch (expected [5,1])');
expected21 = single([11; 22; 33; 44; 55]);
assert(max(abs(double(y21) - double(expected21))) < 1e-6, 'Test 21 failed: values');
fprintf('Test 21 PASSED (ElementwiseAffineLayer align_to_input squeeze [1,1,1,5] -> [5,1])\n');


%% Test 22: PlaceholderLayer Constant (V09 — initializer-as-data input)
L22 = PlaceholderLayer('cls', 'Constant');
L22.Constant = single(reshape(1:48, 1, 1, 48));
y22 = L22.evaluate([]);  % constant doesn't depend on input
assert(isequal(size(y22), [1 1 48]), 'Test 22 failed: shape');
assert(double(y22(1,1,5)) == 5, 'Test 22 failed: value');
fprintf('Test 22 PASSED (PlaceholderLayer Constant value)\n');

%% Test 23: BHWC-to-Flatten transpose pre-permute order (V10 — traffic_signs)
% In ONNX, when a Conv (BCHW) is followed by Transpose [0,2,3,1] -> Flatten,
% the flatten happens on BHWC data (c-fastest C-order). NNV's HWC convention
% with its ONNX-style Flatten produces a w-fastest order by default, so the
% importer now inserts a TransposeLayer with MATLAB perm [2,3,1] BEFORE the
% Flatten so the final flat order matches ONNX.
%
% This test reproduces the data path: HWC -> Transpose [2,3,1] -> Flatten,
% checks that the resulting flat layout reads channel-fastest.

% Build a tiny tensor with distinguishable values
H_t = 2; W_t = 3; C_t = 4;
x_hwc = single(reshape(1:H_t*W_t*C_t, H_t, W_t, C_t));   % MATLAB column-major fill

% Manual reference: ONNX BHWC C-order flat = (h, w, c) iteration with c fastest
ref = zeros(1, H_t*W_t*C_t, 'single');
k = 0;
for h_t = 1:H_t
    for w_t = 1:W_t
        for c_t = 1:C_t
            k = k + 1;
            ref(k) = x_hwc(h_t, w_t, c_t);
        end
    end
end

% NNV path: TransposeLayer [3,4,2] (1-indexed MATLAB perm) then Flatten
ph23 = PlaceholderLayer('t23', 'Transpose');
ph23.Perm = [2, 3, 1];
y_perm = ph23.evaluate(x_hwc);
assert(isequal(size(y_perm), [W_t, C_t, H_t]), 'Test 23 failed: pre-permute shape');

fl23 = FlattenLayer('f23');
fl23.Type = 'nnet.onnx.layer.FlattenInto2dLayer';
y_flat = fl23.evaluate(y_perm);
y_flat_v = squeeze(y_flat);
assert(numel(y_flat_v) == H_t*W_t*C_t, 'Test 23 failed: flat numel');
assert(max(abs(double(y_flat_v(:)) - double(ref(:)))) < 1e-6, ...
    'Test 23 failed: BHWC-Flatten order does not match ONNX C-order');
fprintf('Test 23 PASSED (BHWC-to-Flatten transpose pre-permute)\n');


%% Test 24: ml4acopf-style PWL Sigmoid expansion (V10b)
% Verify the FC-replicate + ElementwiseAffine + ReLU + FC-sum chain that
% the importer emits for the Unsqueeze->Sub->Relu->MatMul->Add Sigmoid PWL
% pattern computes the right values.
%
% Pattern math (per element x_i of input):
%   y_i = b + sum_j (delta_m[j] * max(0, x_i - thresholds[j]))

in_size = 4; N = 3;
thresholds = single([-1; 0; 1]);
delta_m = single([0.5; 0.3; 0.2]);
b_val = single(0.1);
x_in = single([-2; -0.5; 0.5; 2]);

% Reference output
ref = zeros(in_size, 1, 'single');
for ii = 1:in_size
    s = 0;
    for jj = 1:N
        s = s + delta_m(jj) * max(0, x_in(ii) - thresholds(jj));
    end
    ref(ii) = b_val + s;
end

% NNV chain: FC expand, ElementwiseAffine sub, ReLU, FC sum-with-bias
W_exp = zeros(in_size*N, in_size, 'single');
for ii = 1:in_size
    for jj = 1:N
        W_exp((ii-1)*N + jj, ii) = 1;
    end
end
b_exp = zeros(in_size*N, 1, 'single');
fc_exp = FullyConnectedLayer('exp', W_exp, b_exp);
y1 = fc_exp.evaluate(x_in);

tiled_th = repmat(thresholds, in_size, 1);
ea = ElementwiseAffineLayer('sub', single(1), -tiled_th, false, true);
y2 = ea.evaluate(y1);

rl = ReluLayer('relu24');
y3 = rl.evaluate(y2);

W_sum = zeros(in_size, in_size*N, 'single');
for ii = 1:in_size
    for jj = 1:N
        W_sum(ii, (ii-1)*N + jj) = delta_m(jj);
    end
end
b_sum = b_val * ones(in_size, 1, 'single');
fc_sum = FullyConnectedLayer('sum', W_sum, b_sum);
y_out = fc_sum.evaluate(y3);

assert(isequal(size(y_out), size(ref)), 'Test 24 failed: shape');
assert(max(abs(double(y_out(:)) - double(ref(:)))) < 1e-5, ...
    sprintf('Test 24 failed: max diff %g', max(abs(double(y_out(:)) - double(ref(:))))));
fprintf('Test 24 PASSED (PWL Sigmoid expansion: FC + EA + ReLU + FC)\n');


%% Test 25: PlaceholderLayer trig/floor element-wise op handlers (V10c)
% PlaceholderLayer with Type set to Floor/Sin/Cos/Tan/Exp/Log/Sqrt/Ceil/Round
% should compute the actual element-wise op rather than identity-passthrough.

x_t25 = single([-2.7; -0.5; 0; 0.5; 2.7]);

L_floor = PlaceholderLayer('ph_floor', 'Floor');
assert(max(abs(double(L_floor.evaluate(x_t25)) - double(floor(x_t25)))) < 1e-5, 'Test 25 floor');
L_ceil = PlaceholderLayer('ph_ceil', 'Ceil');
assert(max(abs(double(L_ceil.evaluate(x_t25)) - double(ceil(x_t25)))) < 1e-5, 'Test 25 ceil');
L_round = PlaceholderLayer('ph_round', 'Round');
assert(max(abs(double(L_round.evaluate(x_t25)) - double(round(x_t25)))) < 1e-5, 'Test 25 round');
L_sin = PlaceholderLayer('ph_sin', 'Sin');
assert(max(abs(double(L_sin.evaluate(x_t25)) - double(sin(x_t25)))) < 1e-5, 'Test 25 sin');
L_cos = PlaceholderLayer('ph_cos', 'Cos');
assert(max(abs(double(L_cos.evaluate(x_t25)) - double(cos(x_t25)))) < 1e-5, 'Test 25 cos');
L_tan = PlaceholderLayer('ph_tan', 'Tan');
assert(max(abs(double(L_tan.evaluate(x_t25)) - double(tan(x_t25)))) < 1e-5, 'Test 25 tan');
L_exp = PlaceholderLayer('ph_exp', 'Exp');
assert(max(abs(double(L_exp.evaluate(x_t25)) - double(exp(x_t25)))) < 1e-5, 'Test 25 exp');
L_sqrt = PlaceholderLayer('ph_sqrt', 'Sqrt');
x_pos = abs(x_t25);
assert(max(abs(double(L_sqrt.evaluate(x_pos)) - double(sqrt(x_pos)))) < 1e-5, 'Test 25 sqrt');
fprintf('Test 25 PASSED (PlaceholderLayer Floor/Ceil/Round/Sin/Cos/Tan/Exp/Sqrt)\n');

%% Test 26: ElementwiseAffineLayer reach with ND-shaped Offset (V10c)
% ACAS Xu leading norm has Offset shape [1,1,1,5]; prior reach() failed
% because Star.affineMap expects a flat [dim,1] bias.

off26 = single(zeros(1,1,1,5));
off26(1,1,1,1) = 1; off26(1,1,1,2) = 2; off26(1,1,1,3) = 3;
off26(1,1,1,4) = 4; off26(1,1,1,5) = 5;
L26 = ElementwiseAffineLayer('ea26', single(1), off26, false, true);
lb26 = double([0;0;0;0;0]);
ub26 = double([1;1;1;1;1]);
S26 = Star(lb26, ub26);
out26 = L26.reach(S26, 'approx-star', 'single');
[olb, oub] = out26.getRanges();
exp_lb = lb26 + (1:5)';
exp_ub = ub26 + (1:5)';
assert(max(abs(olb - exp_lb)) < 1e-6, 'Test 26 failed: lb mismatch');
assert(max(abs(oub - exp_ub)) < 1e-6, 'Test 26 failed: ub mismatch');
fprintf('Test 26 PASSED (ElementwiseAffineLayer reach with ND offset)\n');


%% Test 27: ElementwiseProductLayer reach is SOUND (MC-containment, signed ranges)
rng(27);
L27 = ElementwiseProductLayer('mul27', 2, 1, {'in1','in2'}, {'out'});
d = 4;
a_lo = randn(d,1) - 0.5; a_hi = a_lo + rand(d,1) + 0.2;   % straddles 0
b_lo = randn(d,1) - 0.5; b_hi = b_lo + rand(d,1) + 0.2;
S27 = L27.reach({Star(a_lo, a_hi), Star(b_lo, b_hi)}, 'approx-star');
[lb27, ub27] = S27.getRanges(); lb27 = lb27(:); ub27 = ub27(:);
viol = 0;
for k = 1:1000
    xa = a_lo + (a_hi - a_lo).*rand(d,1);
    xb = b_lo + (b_hi - b_lo).*rand(d,1);
    y = xa .* xb;
    if any(y < lb27 - 1e-9) || any(y > ub27 + 1e-9), viol = viol + 1; end
end
assert(viol == 0, 'Test 27 failed: product reach UNSOUND (%d/1000 MC violations)', viol);
% wrong arity must error, not silently pass an operand through
errored = false;
try, L27.reach(Star(a_lo, a_hi), 'approx-star'); catch, errored = true; end
assert(errored, 'Test 27 failed: non-cell input must error');
fprintf('Test 27 PASSED (ElementwiseProductLayer MC-sound + fail-loud arity)\n');

%% Test 28: ElementwiseDivisionLayer reach SOUND (MC) + divisor-straddling-zero ERRORS
rng(28);
L28 = ElementwiseDivisionLayer('div28', 2, 1, {'in1','in2'}, {'out'});
d = 4;
a_lo = randn(d,1); a_hi = a_lo + rand(d,1) + 0.2;          % numerator: any sign
b_lo = -2 - rand(d,1); b_hi = b_lo + rand(d,1);            % divisor strictly NEGATIVE
S28 = L28.reach({Star(a_lo, a_hi), Star(b_lo, b_hi)}, 'approx-star');
[lb28, ub28] = S28.getRanges(); lb28 = lb28(:); ub28 = ub28(:);
viol = 0;
for k = 1:1000
    xa = a_lo + (a_hi - a_lo).*rand(d,1);
    xb = b_lo + (b_hi - b_lo).*rand(d,1);
    y = xa ./ xb;
    if any(y < lb28 - 1e-9) || any(y > ub28 + 1e-9), viol = viol + 1; end
end
assert(viol == 0, 'Test 28 failed: division reach UNSOUND (%d/1000 MC violations)', viol);
% divisor interval straddling zero: quotient unbounded -> MUST error (the old
% behavior warned "reach unsound" and passed the numerator through!)
errored = false;
try
    L28.reach({Star(a_lo, a_hi), Star(-ones(d,1), ones(d,1))}, 'approx-star');
catch ME28
    errored = strcmp(ME28.identifier, 'ElementwiseDivisionLayer:divisorStraddlesZero');
end
assert(errored, 'Test 28 failed: straddle-zero divisor must raise divisorStraddlesZero');
fprintf('Test 28 PASSED (ElementwiseDivisionLayer MC-sound + straddle-zero fail-loud)\n');

%% Test 29: DynamicMatmulLayer reach ERRORS (no sound bound implemented); evaluate exact
% The previous reach took ELEMENT-WISE interval products of the flattened
% operands -- not a sound bound for a matrix product (out(i,j) sums k bilinear
% terms and the output shape differs). It must refuse rather than mislead.
L29 = DynamicMatmulLayer('mm29', 2, 1, {'in1','in2'}, {'out'});
errored = false;
try
    L29.reach({Star(zeros(4,1), ones(4,1)), Star(zeros(4,1), ones(4,1))}, 'approx-star');
catch ME29
    errored = strcmp(ME29.identifier, 'DynamicMatmulLayer:reachNotImplemented');
end
assert(errored, 'Test 29 failed: DynamicMatmul reach must raise reachNotImplemented');
A29 = rand(3,4,'single'); B29 = rand(4,2,'single');
assert(max(abs(double(L29.evaluate({A29,B29}) - A29*B29)), [], 'all') < 1e-5, ...
    'Test 29 failed: evaluate must stay exact');
fprintf('Test 29 PASSED (DynamicMatmulLayer reach fail-loud; evaluate exact)\n');

%% Test 30: PlaceholderLayer ACTIVE ops reach soundly or refuse (identity was unsound)
% evaluate() computes Sign/Abs/Exp/... and Constant/Transpose, but reach() was
% an unconditional identity -> evaluate and reach disagreed (silent unsoundness).
rng(30);
d = 5; lo = randn(d,1) - 0.5; hi = lo + rand(d,1) + 0.3;   % straddles 0
S30 = Star(lo, hi);
% (a) Exp: monotone box, MC-containment
Lx = PlaceholderLayer('ph_exp30', 'Exp');
Sx = Lx.reach(S30, 'approx-star');
[xl, xu] = Sx.getRanges(); xl = xl(:); xu = xu(:);
viol = 0;
for k = 1:500
    x = lo + (hi - lo).*rand(d,1); y = exp(x);
    if any(y < xl - 1e-9) || any(y > xu + 1e-9), viol = viol + 1; end
end
assert(viol == 0, 'Test 30a failed: Exp reach unsound (%d/500)', viol);
% (b) Abs: interval abs, MC-containment
La = PlaceholderLayer('ph_abs30', 'Abs');
Sa30 = La.reach(S30, 'approx-star');
[al, au] = Sa30.getRanges(); al = al(:); au = au(:);
viol = 0;
for k = 1:500
    x = lo + (hi - lo).*rand(d,1); y = abs(x);
    if any(y < al - 1e-9) || any(y > au + 1e-9), viol = viol + 1; end
end
assert(viol == 0, 'Test 30b failed: Abs reach unsound (%d/500)', viol);
% (c) Sin: sound codomain box [-1,1]
Ls = PlaceholderLayer('ph_sin30', 'Sin');
Ss30 = Ls.reach(S30, 'approx-star');
[sl, su] = Ss30.getRanges();
assert(all(sl(:) <= -1 + 1e-9) && all(su(:) >= 1 - 1e-9), 'Test 30c failed: Sin box');
% (d) Constant: exact point set
Lc = PlaceholderLayer('ph_const30', 'Constant');
Lc.Constant = single([1; -2; 3]);
Sc30 = Lc.reach(S30, 'approx-star');
[cl, cu] = Sc30.getRanges();
assert(max(abs(cl(:) - [1; -2; 3])) < 1e-6 && max(abs(cu(:) - [1; -2; 3])) < 1e-6, ...
    'Test 30d failed: Constant point set');
% (e) Tan and non-identity Transpose must REFUSE
Lt = PlaceholderLayer('ph_tan30', 'Tan');
errored = false; try, Lt.reach(S30, 'approx-star'); catch, errored = true; end
assert(errored, 'Test 30e failed: Tan reach must error');
Lp = PlaceholderLayer('ph_perm30', 'Transpose'); Lp.Perm = [2 1];
errored = false; try, Lp.reach(S30, 'approx-star'); catch, errored = true; end
assert(errored, 'Test 30e failed: non-identity Transpose reach must error');
% (f) inert placeholder stays identity
Ld = PlaceholderLayer('ph_drop30', 'nnet.cnn.layer.DropoutLayer');
Sout = Ld.reach(S30, 'approx-star');
[il, iu] = Sout.getRanges();
assert(max(abs(il(:) - lo)) < 1e-9 && max(abs(iu(:) - hi)) < 1e-9, ...
    'Test 30f failed: inert placeholder must stay identity');
fprintf('Test 30 PASSED (PlaceholderLayer active ops: sound reach or refuse; inert identity)\n');

%% Test 31: PlaceholderLayer 'UnsupportedOp:<op>' REFUSES in evaluate AND reach
% The loader tags active ops it can't represent (Where/ScatterND/ArgMax/Expand/
% dynamic Reshape/...) as 'UnsupportedOp:<op>'. Treating them as identity would
% silently evaluate/verify a DIFFERENT network (the cctsdb_yolo / collins class).
% Both paths must error -- regression so the tag can never silently pass through.
Lu = PlaceholderLayer('ph_unsup', 'UnsupportedOp:Where');
errEval = false;
try, Lu.evaluate(single([1; 2; 3])); catch ME31a, errEval = strcmp(ME31a.identifier,'PlaceholderLayer:unsupportedOp'); end
assert(errEval, 'Test 31 failed: UnsupportedOp evaluate must error (identity would be a wrong network)');
errReach = false;
try, Lu.reach(Star(zeros(3,1), ones(3,1)), 'approx-star'); catch ME31b, errReach = strcmp(ME31b.identifier,'PlaceholderLayer:unsupportedOp'); end
assert(errReach, 'Test 31 failed: UnsupportedOp reach must error');
% reachSequence path must refuse too
errSeq = false;
try, Lu.reachSequence(Star(zeros(3,1), ones(3,1)), 'approx-star'); catch, errSeq = true; end
assert(errSeq, 'Test 31 failed: UnsupportedOp reachSequence must error');
fprintf('Test 31 PASSED (PlaceholderLayer UnsupportedOp refuses in evaluate + reach + reachSequence)\n');

fprintf('\n=== All V04/V05/V06/V07/V08/V09/V10 NNV regression tests PASSED ===\n');
