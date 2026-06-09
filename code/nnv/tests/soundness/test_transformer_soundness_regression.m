%% Transformer soundness regression tests
% Locks in the soundness fixes for Softmax (mid-network), attention (multi-token
% fail-loud), and confirms BatchNorm soundness. Runs in matrix CI (fast, <5s).
%
% Each %% section is SELF-CONTAINED (runtests runs each in a fresh workspace), so
% do not rely on variables defined in other sections.

%% Test 1: intermediate softmax reach is sound (MC-containment) and in [0,1]
rng(0); d = 4; c = randn(d,1); ev = 0.5;
I = Star(c - ev, c + ev);
L = SoftmaxLayer('sm'); L.IsFinalLayer = false;   % intermediate (e.g. attention) softmax
S = L.reach(I);
[lb, ub] = S.getRanges; lb = lb(:); ub = ub(:);
assert(all(lb >= -1e-9) && all(ub <= 1 + 1e-9), ...
    'intermediate softmax reach not within codomain [0,1]');
viol = 0;
for k = 1:1000
    x = c + ev*(2*rand(d,1) - 1);
    y = L.evaluate(x); y = y(:);
    if any(y < lb - 1e-6) || any(y > ub + 1e-6), viol = viol + 1; end
end
assert(viol == 0, 'intermediate softmax UNSOUND: %d/1000 MC violations', viol);
fprintf('Test 1 PASSED (intermediate softmax sound + in [0,1])\n');

%% Test 2: final-layer softmax reach is identity (legacy, sound for argmax-on-logits)
rng(0); d = 4; c = randn(d,1); ev = 0.5;
I = Star(c - ev, c + ev);
L = SoftmaxLayer('sm');   % default IsFinalLayer = true
S = L.reach(I);
[lbi, ubi] = I.getRanges; [lbs, ubs] = S.getRanges;
assert(max(abs(lbi(:) - lbs(:))) < 1e-9 && max(abs(ubi(:) - ubs(:))) < 1e-9, ...
    'final-layer softmax reach should be identity');
fprintf('Test 2 PASSED (final-layer softmax identity)\n');

%% Test 3: intermediate softmax errors on unsupported set type (no silent unsound)
L = SoftmaxLayer('sm'); L.IsFinalLayer = false;
errored = false;
try, L.reach(42); catch, errored = true; end
assert(errored, 'intermediate softmax must error on unsupported set type, not silently passthrough');
fprintf('Test 3 PASSED (softmax fails loud on unsupported type)\n');

%% Test 4: SDPA multi-token reach errors (would be unsound); single-token still works
L = ScaledDotProductAttentionLayer('a');
L.ValueDim = 4; L.QueryDim = 4; L.KeyDim = 4; L.Scale = 0.5;
Q = Star(zeros(4,1), ones(4,1)); K = Q;
errMulti = false;
try, L.reach(Q, K, Star(zeros(8,1), ones(8,1)), 'approx-star'); catch, errMulti = true; end
assert(errMulti, 'multi-token SDPA reach must error (single-token V-bounds would be unsound)');
S = L.reach(Q, K, Star(zeros(4,1), ones(4,1)), 'approx-star');
assert(S.dim == 4, 'single-token SDPA reach broken');
fprintf('Test 4 PASSED (SDPA multi-token fail-loud, single-token ok)\n');

%% Test 5: MHA multi-token reach errors (would be unsound)
L = MultiHeadAttentionLayer();
L.NumHeads = 1; L.EmbedDim = 4; L.HeadDim = 4;
L.W_Q = eye(4); L.W_K = eye(4); L.W_V = eye(4); L.W_O = eye(4);
errMulti = false;
try, L.reach(Star(zeros(8,1), ones(8,1)), 'approx-star'); catch, errMulti = true; end
assert(errMulti, 'multi-token MHA reach must error (single-token V-bounds would be unsound)');
fprintf('Test 5 PASSED (MHA multi-token fail-loud)\n');

%% Test 6: BatchNorm reach is sound (MC-containment), ImageStar path
rng(0); H = 2; W = 2; C = 3;
mean_  = reshape([0.5 -0.2 1.0], 1, 1, C);
var_   = reshape([0.3  0.5 0.2], 1, 1, C);
scale  = reshape([1.2  0.8 1.5], 1, 1, C);
offset = reshape([0.1 -0.1 0.2], 1, 1, C);
L = BatchNormalizationLayer('Name','bn','NumChannels',C,'TrainedMean',mean_, ...
    'TrainedVariance',var_,'Scale',scale,'Offset',offset,'Epsilon',1e-5);
lb = rand(H,W,C)*0.5; ub = lb + 0.4;
IS = ImageStar(lb, ub);
OS = L.reach(IS, 'approx-star');
[olb, oub] = OS.getRanges; olb = olb(:); oub = oub(:);
viol = 0;
for k = 1:1000
    x = lb + (ub - lb).*rand(H,W,C);
    yv = L.evaluate(x); yv = yv(:);
    if any(yv < olb - 1e-6) || any(yv > oub + 1e-6), viol = viol + 1; end
end
assert(viol == 0, 'BatchNorm reach UNSOUND: %d/1000 MC violations', viol);
fprintf('Test 6 PASSED (BatchNorm ImageStar sound)\n');
