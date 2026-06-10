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

%% Test 7: MHA with weights set but UNRESOLVED dims (EmbedDim=0) must still
%% fail loud on multi-token input -- EmbedDim=0 previously DISARMED the guard
%% (the real-attention ViT bug: R2026a selfAttentionLayer property names don't
%% match parse's reads, so EmbedDim stayed 0, HeadDim=0 gave an empty per-head
%% index range, and the output set came back EMPTY -> garbage downstream).
L = MultiHeadAttentionLayer();
L.NumHeads = 2; L.EmbedDim = 0; L.HeadDim = 0;   % simulate failed property read
L.W_Q = eye(4); L.W_K = eye(4); L.W_V = eye(4); L.W_O = eye(4);
errored = false;
try, L.reach(Star(zeros(8,1), ones(8,1)), 'approx-star'); catch, errored = true; end
assert(errored, 'dims-unset MHA on multi-token input must error, not return an empty/garbage set');
fprintf('Test 7 PASSED (MHA derives dims from weights; multi-token still fail-loud)\n');

%% Test 8: parse derives EmbedDim/HeadDim from weights on R2026a selfAttentionLayer
%% (whose property names are InputSize/NumQueryChannels, not NumChannels).
dln = dlnetwork([sequenceInputLayer(8) selfAttentionLayer(2, 8, 'Name', 'sa')]);
L = MultiHeadAttentionLayer.parse(dln.Layers(2));
assert(~isempty(L.W_Q) && isequal(size(L.W_Q), [8 8]), 'parse failed to extract weights');
assert(L.EmbedDim == 8 && L.HeadDim == 4, ...
    'parse failed to derive dims from weights (EmbedDim=%g, HeadDim=%g)', L.EmbedDim, L.HeadDim);
fprintf('Test 8 PASSED (parse derives EmbedDim/HeadDim from extracted weights)\n');

%% Test 9: matlab2nnv mid-network softmax -> SOUND SoftmaxLayer (not identity
%% placeholder); final softmax stays a placeholder (legacy, argmax-on-logits).
dln = dlnetwork([featureInputLayer(4, 'Name', 'in') ...
                 fullyConnectedLayer(4, 'Name', 'fc1') ...
                 softmaxLayer('Name', 'sm_mid') ...
                 fullyConnectedLayer(3, 'Name', 'fc2') ...
                 softmaxLayer('Name', 'sm_out')]);
nnvnet = matlab2nnv(dln);
Ls = nnvnet.Layers; if ~iscell(Ls), Ls = num2cell(Ls); end
idx_mid = 0; idx_out = 0;
for k = 1:numel(Ls)
    nm = ''; if isprop(Ls{k}, 'Name'), nm = Ls{k}.Name; end
    if strcmp(nm, 'sm_mid'), idx_mid = k; end
    if strcmp(nm, 'sm_out'), idx_out = k; end
end
assert(idx_mid > 0 && idx_out > 0, 'softmax layers not found after conversion');
assert(isa(Ls{idx_mid}, 'SoftmaxLayer') && ~Ls{idx_mid}.IsFinalLayer, ...
    'mid-network MATLAB softmax must convert to a non-final (sound) SoftmaxLayer, got %s', class(Ls{idx_mid}));
assert(isa(Ls{idx_out}, 'PlaceholderLayer'), ...
    'final MATLAB softmax should remain an identity placeholder, got %s', class(Ls{idx_out}));
fprintf('Test 9 PASSED (matlab2nnv: mid softmax sound, final softmax identity)\n');

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

%% Test 10: LayerNorm reach is SOUND -- regression for the unsound var-of-center bug
% The old code estimated variance from the CENTER point only (var_center*0.5 /
% max(.)^2/4), giving std_lb=sqrt(eps) and a lower bound that EXCLUDED reachable
% outputs (2532/5000 MC violations). Concrete witness below + randomized MC.
% (a) the exact Codex counterexample: x in [0,2]x[-3,5], witness x=[0;5]->[-1,1].
L = LayerNormalizationLayer('Name','ln','NumChannels',2,'Epsilon',1e-5,'Scale',[1 1],'Offset',[0 0]);
S = Star([0; -3], [2; 5]);
OS = L.reach(S, 'approx-star');
[olb, oub] = OS.getRanges; olb = olb(:); oub = oub(:);
yw = L.evaluate([0; 5]); yw = yw(:);
assert(all(yw >= olb - 1e-6) && all(yw <= oub + 1e-6), ...
    'Test 10a failed: LayerNorm reach excludes the witness x=[0;5] (lb(1)=%.3f, y(1)=%.3f)', olb(1), yw(1));
% (b) randomized MC over varied n / centers / radii / scale-signs / eps
rng(10); totalViol = 0; cases = 0;
for t = 1:120
    n = randi([1 6]); c = randn(n,1)*3; r = rand(n,1)*4 + 0.05; lb = c - r; ub = c + r;
    sc = randn(n,1); of = randn(n,1); ep = 10^(-randi([1 6]));
    Lr = LayerNormalizationLayer('Name','l','NumChannels',n,'Epsilon',ep,'Scale',sc,'Offset',of);
    Or = Lr.reach(Star(lb, ub), 'approx-star'); [a, b] = Or.getRanges; a = a(:); b = b(:);
    for k = 1:30
        x = lb + (ub - lb).*rand(n,1); y = Lr.evaluate(x); y = y(:);
        if any(y < a - 1e-6) || any(y > b + 1e-6), totalViol = totalViol + 1; end
        cases = cases + 1;
    end
end
assert(totalViol == 0, 'Test 10b failed: LayerNorm reach UNSOUND (%d/%d MC violations)', totalViol, cases);
fprintf('Test 10 PASSED (LayerNorm reach sound: witness + %d randomized MC)\n', cases);
