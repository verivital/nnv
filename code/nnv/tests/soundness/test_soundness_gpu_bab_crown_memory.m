% test_soundness_gpu_bab_crown_memory
% Soundness + equivalence tests for the GPU-BaB conv MEMORY handling in gpu_bab_crown_tight:
%   (1) chunked eye(nk) intermediate-bound seed (NNV_CROWN_CHUNK) -- bit-identical (B>=2) + sound (B=1);
%   (2) first-layer IBP short-circuit (box-exact inputs, NNV_CROWN_NO_SHORTCUT kill-switch) -- transparent
%       (== full CROWN, since IBP of an affine map of a box is exact) + sound + no looser than IBP.
% Net is small so the dense eye(nk) pass and the chunked/short-circuit paths all fit and can be compared.

rng(7);
layers = [imageInputLayer([8 8 3],'Normalization','none','Name','in')
          convolution2dLayer(3,5,'Padding',1,'Name','c1'); reluLayer('Name','r1')
          convolution2dLayer(3,4,'Padding',1,'Name','c2'); reluLayer('Name','r2')
          fullyConnectedLayer(6,'Name','fc')];
nnvnet = matlab2nnv(dlnetwork(layers));
ops = nn_to_ops(nnvnet);
ish = [8 8 3];
C   = [eye(5) -ones(5,1)];
xc  = rand(ish); lb = xc(:) - 0.05; ub = xc(:) + 0.05;
Xs  = lb + (ub - lb) .* rand(numel(lb), 4000);
tmp = inf(5,1);
for s = 1:size(Xs,2), y = nnvnet.evaluate(reshape(Xs(:,s), ish)); tmp = min(tmp, C * y(:)); end
trueMin = tmp;

%% Test 1: chunked seed (B>=2) is BIT-IDENTICAL to the dense single eye(nk) pass
setenv('NNV_CROWN_CHUNK','1000000'); mFull = gpu_bab_crown_tight(ops, lb, ub, C, 'double');
for B = [2 3 7]
    setenv('NNV_CROWN_CHUNK', num2str(B));
    mB = gpu_bab_crown_tight(ops, lb, ub, C, 'double');
    assert(max(abs(mB(:) - mFull(:))) == 0, sprintf('chunked B=%d must be bit-identical to the dense pass', B));
end
setenv('NNV_CROWN_CHUNK','');

%% Test 2: chunking is SOUND under maximal chunking (B=1)
setenv('NNV_CROWN_CHUNK','1'); m1 = gpu_bab_crown_tight(ops, lb, ub, C, 'double'); setenv('NNV_CROWN_CHUNK','');
assert(all(m1(:) <= trueMin + 1e-6), 'chunked B=1 margin must be sound (<= true min)');

%% Test 3: first-layer IBP short-circuit is TRANSPARENT (== full CROWN, to FP rounding)
setenv('NNV_CROWN_NO_SHORTCUT','');  mON  = gpu_bab_crown_tight(ops, lb, ub, C, 'double');
setenv('NNV_CROWN_NO_SHORTCUT','1'); mOFF = gpu_bab_crown_tight(ops, lb, ub, C, 'double');
setenv('NNV_CROWN_NO_SHORTCUT','');
assert(max(abs(mON(:) - mOFF(:))) < 1e-9, 'IBP short-circuit must match full CROWN (box-exact => exact)');

%% Test 4: short-circuit margins are SOUND and no looser than IBP
mON = gpu_bab_crown_tight(ops, lb, ub, C, 'double');   % short-circuit is on by default
assert(all(mON(:) <= trueMin + 1e-6), 'short-circuit margin must be sound (<= true min)');
[ilb, iub] = gpu_bab_ibp(ops, lb, ub, 'double');
assert(all(mON(:) >= max(C,0)*ilb + min(C,0)*iub - 1e-6), 'short-circuit CROWN must be no looser than IBP');

%% Test 5: gpu_bab_ibp box-exactness flags the input (true) and a coordinate-mixing op (false)
[~, ~, ~, ~, ~, boxExact] = gpu_bab_ibp(ops, lb, ub, 'double');
assert(~isempty(boxExact) && boxExact{1} == true, 'the input set must be box-exact');
assert(any(cellfun(@(v) ~isempty(v) && ~all(logical(v)), boxExact)), ...
    'a coordinate-mixing op (conv) must break box-exactness downstream');

%% Test 6: chunking + short-circuit COMPOSED is still sound (the production default config)
setenv('NNV_CROWN_CHUNK','2');  mc = gpu_bab_crown_tight(ops, lb, ub, C, 'double'); setenv('NNV_CROWN_CHUNK','');
assert(all(mc(:) <= trueMin + 1e-6), 'chunk+short-circuit composed margin must be sound');

%% Test 7: AUTO-CHUNK halve-on-OOM yields the SAME bounds as a clean run (retry is exact + B-invariant)
% NNV_CROWN_TEST_OOM forces a synthetic device OOM while B > thr, so the catch/halve/retry path runs:
% B steps 64 -> 32 -> 16 -> 8, then proceeds. The result must be bit-identical to a clean B=64 pass.
setenv('NNV_CROWN_CHUNK','64'); mClean = gpu_bab_crown_tight(ops, lb, ub, C, 'double');
setenv('NNV_CROWN_TEST_OOM','8'); mHalved = gpu_bab_crown_tight(ops, lb, ub, C, 'double');
setenv('NNV_CROWN_TEST_OOM',''); setenv('NNV_CROWN_CHUNK','');
assert(max(abs(mHalved(:) - mClean(:))) == 0, 'auto-chunk halved result must be bit-identical to a clean run');
assert(all(mHalved(:) <= trueMin + 1e-6), 'auto-chunk halved margin must be sound');

%% Test 8: a NON-memory error inside the chunk loop PROPAGATES (i_is_oom must not swallow real bugs)
setenv('NNV_CROWN_CHUNK','64'); setenv('NNV_CROWN_TEST_OOM','bug');
threw = false;
try, gpu_bab_crown_tight(ops, lb, ub, C, 'double'); catch ME, threw = strcmp(ME.identifier, 'NNV:test:notMemory'); end
setenv('NNV_CROWN_TEST_OOM',''); setenv('NNV_CROWN_CHUNK','');
assert(threw, 'a non-memory error must propagate out of the auto-chunk retry (not be retried as OOM)');

%% Summary
disp('test_soundness_gpu_bab_crown_memory: all sections passed');
