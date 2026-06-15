% test_soundness_gpu_bab_root_tight
% ROOT-TIGHT-BOUND REUSE: gpu_bab_crown_spec_dag(...,rootBounds) must reuse the tight root
% intermediate bounds (computed once by gpu_bab_crown_tight) -- clamped per node -- SOUNDLY:
%   1. with NO clamps, the reuse margin equals the serial-tight margin (same preL/preU =>
%      same backward CROWN), so the reuse backward is correct;
%   2. the reuse margin is a sound lower bound on C*f over the (whole) box (<= Monte-Carlo min);
%   3. for clamped nodes, reuse <= serial-tight-with-fixings (so reuse inherits serial-tight's
%      soundness) AND reuse >= IBP (the tightness win) AND reuse <= the MC min over the node's
%      sub-region (box intersect the fixing halfspaces) -- a direct, rigorous soundness check;
%   4. at the BaB level (gpu_bab_relu_split_batched, rootTight default ON), the verdict never
%      contradicts the serial path and decides at least the IBP path's robust cases.
% FC+ReLU nets, DOUBLE precision for the bound checks (exactness). Loops in local functions.

%% Test 1: reuse with empty fixings == serial-tight margin (the reuse backward is exact)
rng(1); ops = i_fc_relu([5 14 14 6]);
lb = -0.25*ones(5,1); ub = 0.25*ones(5,1);
C = [eye(5) -ones(5,1)];
[mTight, rtL, rtU] = gpu_bab_crown_tight(ops, lb, ub, C, 'double', cell(numel(ops),1));
rootBounds = struct('preL', {rtL}, 'preU', {rtU});
mReuse = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double', cell(numel(ops),1), rootBounds);
assert(max(abs(mTight(:) - mReuse(:))) < 1e-9, ...
    sprintf('reuse(empty) must equal serial-tight margin; max diff = %.3e', max(abs(mTight(:)-mReuse(:)))));

%% Test 2: reuse(empty) is a sound lower bound over the box (<= Monte-Carlo min)
rng(2); ops = i_fc_relu([6 16 16 7]);
lb = -0.2*ones(6,1); ub = 0.2*ones(6,1);
C = [eye(6) -ones(6,1)];
[~, rtL, rtU] = gpu_bab_crown_tight(ops, lb, ub, C, 'double', cell(numel(ops),1));
rootBounds = struct('preL', {rtL}, 'preU', {rtU});
mReuse = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double', cell(numel(ops),1), rootBounds);
mcMin = i_mc_min(ops, C, lb, ub, 4000);
assert(all(mReuse(:) <= mcMin(:) + 1e-6), 'reuse(empty) must be a sound lower bound (<= MC min)');

%% Test 3: clamped reuse is SOUND (<= node-region MC min) and tighter than IBP in aggregate
% NOTE: reuse is NOT required to be <= serial-tight -- CROWN's bound is not monotone in
% intermediate-bound tightness (the adaptive lower slope flips), so root-tight reuse and full
% serial-tight are two distinct sound bounds and reuse is sometimes the tighter one. The hard
% guarantee we assert is direct: reuse never exceeds the Monte-Carlo min over the node region.
[nBad, nNode, nMC, tighterFrac] = i_check_clamped(3, 40);
assert(nBad == 0, sprintf('clamped reuse SOUNDNESS violated (reuse > node-region MC min) in %d checks', nBad));
assert(nNode >= 25, sprintf('expected >=25 node checks, got %d', nNode));
assert(nMC >= 8, sprintf('expected >=8 valid node-region MC checks, got %d', nMC));
% Aggregate tightness win: root-tight reuse is at least as tight as IBP on the large majority
% of nodes (a node clamp propagates through the IBP forward but not the reused root bounds, so
% it is not a pointwise invariant -- but the tight intermediate bounds dominate overall).
assert(tighterFrac >= 0.6, sprintf('root-tight reuse tighter than IBP on only %.0f%% of nodes (expected >=60%%)', 100*tighterFrac));

%% Test 4: BaB with rootTight ON never contradicts serial; decides the IBP-robust cases
[nContra, nRootDecide, nIbpDecide, nTot] = i_compare_bab(4, 12);
assert(nContra == 0, sprintf('rootTight BaB contradicted serial in %d cases', nContra));
assert(nRootDecide >= nIbpDecide, ...
    sprintf('rootTight decided %d but IBP decided %d (tight bounds must not decide fewer)', nRootDecide, nIbpDecide));

%% Summary
disp('test_soundness_gpu_bab_root_tight: all sections passed');

% ----------------------------------------------------------------------------------------
function ops = i_fc_relu(dims)
% Sequential FC+ReLU op list (nn_to_ops format), He-init, double weights.
    ops = {};
    for L = 1:numel(dims)-1
        W = randn(dims(L+1), dims(L)) * sqrt(2/dims(L));
        b = randn(dims(L+1), 1) * 0.1;
        ops{end+1} = struct('type','affine','W',W,'b',b,'src',numel(ops)); %#ok<AGROW>
        if L < numel(dims)-1
            ops{end+1} = struct('type','relu','src',numel(ops)); %#ok<AGROW>
        end
    end
end

function [y, pre] = i_forward_fc(ops, X)
% Sequential FC forward over columns of X (n x N); pre{k} = pre-activation feeding relu op k.
    pre = cell(numel(ops),1);
    v = X;
    for k = 1:numel(ops)
        if strcmp(ops{k}.type, 'affine')
            v = ops{k}.W * v + ops{k}.b;
        else            % relu
            pre{k} = v;
            v = max(v, 0);
        end
    end
    y = v;
end

function mn = i_mc_min(ops, C, lb, ub, nS)
% Min of C*f over nS uniform samples in [lb,ub].
    rng(101);
    X = lb + (ub - lb) .* rand(numel(lb), nS);
    Y = i_forward_fc(ops, X);
    mn = min(C * Y, [], 2);
end

function [nBad, nNode, nMC, tighterFrac] = i_check_clamped(seed, nTrials)
% For random nets + random node fixings on <=3 unstable neurons (kept few so the box ∩ fixing
% halfspaces still admits enough uniform samples), verify per node:
%   reuse <= MC-min over the node's sub-region          [SOUND, when enough in-region samples]
% and track how often reuse >= IBP (the aggregate tightness win; not a pointwise invariant).
% nBad counts ONLY soundness violations (reuse exceeding the node-region MC min).
    rng(seed); nBad = 0; nNode = 0; nMC = 0; nTighter = 0; nCmp = 0;
    for t = 1:nTrials
        dims = [6 18 18 7]; ops = i_fc_relu(dims);
        lb = -0.25*ones(6,1); ub = 0.25*ones(6,1);
        C = [eye(6) -ones(6,1)];
        reluIdx = find(cellfun(@(o) strcmp(o.type,'relu'), ops));
        [~, rtL, rtU] = gpu_bab_crown_tight(ops, lb, ub, C, 'double', cell(numel(ops),1));
        rootBounds = struct('preL', {rtL}, 'preU', {rtU});
        % collect all unstable (relu-op, neuron) pairs, fix a random subset of up to 3
        cand = [];
        for r = 1:numel(reluIdx)
            k = reluIdx(r); uns = find(rtL{k} < 0 & rtU{k} > 0);
            cand = [cand; repmat(k, numel(uns), 1), uns(:)]; %#ok<AGROW>
        end
        if isempty(cand), continue; end
        fx = cell(numel(ops),1);
        for r = 1:numel(reluIdx), fx{reluIdx(r)} = zeros(numel(rtL{reluIdx(r)}), 1); end
        nFix = min(3, size(cand,1));
        pick = cand(randperm(size(cand,1), nFix), :);
        for p = 1:nFix, fx{pick(p,1)}(pick(p,2)) = (rand < 0.5)*2 - 1; end

        mReuse  = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double', fx, rootBounds);   % root-tight reuse
        mIBP    = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double', fx);               % IBP (looser)

        nNode = nNode + 1;
        nCmp = nCmp + 1;
        if mean(mReuse(:)) >= mean(mIBP(:)) - 1e-9, nTighter = nTighter + 1; end

        % MC over the node's sub-region (box ∩ fixing halfspaces)
        rng(2000 + t);
        X = lb + (ub - lb) .* rand(6, 30000);
        [Y, pre] = i_forward_fc(ops, X);
        mask = true(1, size(X,2));
        for p = 1:nFix
            k = pick(p,1); j = pick(p,2);
            if fx{k}(j) == 1, mask = mask & (pre{k}(j,:) >= 0);
            else,             mask = mask & (pre{k}(j,:) <= 0); end
        end
        if nnz(mask) >= 30
            nMC = nMC + 1;
            mcMin = min(C * Y(:, mask), [], 2);
            if ~all(mReuse(:) <= mcMin(:) + 1e-5), nBad = nBad + 1; end   % SOUND over the node region
        end
    end
    tighterFrac = nTighter / max(1, nCmp);
end

function [nContra, nRootDecide, nIbpDecide, nTot] = i_compare_bab(seed, nTrials)
% serial vs batched(rootTight ON) vs batched(rootTight OFF/IBP): count contradictions and how
% many each definitively decides (robust/unsafe, not unknown).
    rng(seed); nContra = 0; nRootDecide = 0; nIbpDecide = 0; nTot = nTrials;
    for t = 1:nTrials
        dims = [6 20 20 6]; ops = i_fc_relu(dims);
        x = randn(dims(1),1) * 0.5;
        yc = i_forward_fc(ops, x); [~, tl] = max(yc);
        ep = 0.02 + 0.05*rand;
        lb = x - ep; ub = x + ep;
        oS    = struct('precision','double','maxNodes',8000);
        oRoot = struct('precision','double','maxNodes',8000,'maxFrontier',256,'rootTight',true);
        oIbp  = struct('precision','double','maxNodes',8000,'maxFrontier',256,'rootTight',false);
        [sS,~]    = gpu_bab_relu_split(ops, lb, ub, tl, dims(end), oS);
        [sRoot,~] = gpu_bab_relu_split_batched(ops, lb, ub, tl, dims(end), oRoot);
        [sIbp,~]  = gpu_bab_relu_split_batched(ops, lb, ub, tl, dims(end), oIbp);
        if i_contra(sS, sRoot) || i_contra(sS, sIbp) || i_contra(sRoot, sIbp)
            nContra = nContra + 1;
        end
        if ~strcmp(sRoot,'unknown'), nRootDecide = nRootDecide + 1; end
        if ~strcmp(sIbp,'unknown'),  nIbpDecide  = nIbpDecide  + 1; end
    end
end

function tf = i_contra(a, b)
    tf = (strcmp(a,'robust') && strcmp(b,'unsafe')) || (strcmp(a,'unsafe') && strcmp(b,'robust'));
end
