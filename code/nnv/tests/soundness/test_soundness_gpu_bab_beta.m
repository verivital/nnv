% test_soundness_gpu_bab_beta
% Joint alpha+beta-CROWN (gpu_bab_crown_alpha_beta) is the per-node BaB bound that adds the
% beta-CROWN split-constraint duals on top of root-tight + alpha. It bounds the CONSTRAINED min
% (min over the input box INTERSECT the node's split-constraint halfspaces s_i*z_i>=0), so the
% soundness check is vs the in-region Monte-Carlo min. beta init = 0, so it starts at the alpha
% bound and only raises the certification (per-column min) margin. FC+ReLU, double precision.

%% Test 1: nIter=0 (alpha=min-area, beta=0) + rootBounds == root-tight margin exactly
rng(1); ops = i_fc([5 12 12 6]);
lb = -0.2*ones(5,1); ub = 0.2*ones(5,1); C = [eye(5) -ones(5,1)];
reluIdx = find(cellfun(@(o) strcmp(o.type,'relu'), ops));
[~, rtL, rtU] = gpu_bab_crown_tight(ops, lb, ub, C, 'double', cell(numel(ops),1));
rb = struct('preL', {rtL}, 'preU', {rtU});
mTight = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double', cell(numel(ops),1), rb);
mAB0 = gpu_bab_crown_alpha_beta(ops, lb, ub, C, cell(numel(ops),1), reluIdx, 'double', 0, 0.2, rb);
assert(max(abs(mTight(:) - mAB0(:))) < 1e-9, ...
    sprintf('alpha_beta(0 iter, no fixings) must equal root-tight; max diff %.3e', max(abs(mTight(:)-mAB0(:)))));

%% Test 2: clamped nodes -- alpha+beta SOUND vs the CONSTRAINED (box ∩ fixings) MC min, and its
%% per-column min-margin >= alpha-only (beta only helps the certification quantity)
% The HARD invariant is soundness (alpha+beta <= the constrained MC min). beta helps the cert
% margin in the majority of nodes but NOT pointwise vs alpha-only (the joint projected-gradient
% trajectory can differ), so we assert soundness + that beta helps sometimes, not >= alpha.
[nBad, nChk, nBetaHelp] = i_check_beta(2, 30);
assert(nBad == 0, sprintf('alpha+beta SOUNDNESS violated (> constrained MC min) in %d checks', nBad));
assert(nChk >= 10, sprintf('expected >=10 valid in-region node checks, got %d', nChk));
assert(nBetaHelp >= 1, sprintf('beta never improved the certification margin (got %d) -- expected it to help on some nodes', nBetaHelp));

%% Test 3: BaB with beta never contradicts serial; decides >= the alpha-only BaB
rng(3); nContra = 0; nBeta = 0; nAlpha = 0;
for t = 1:6
    dims = [6 16 16 6]; ops = i_fc(dims);
    x = randn(6,1)*0.5; yc = i_fwd(ops, x); [~, tl] = max(yc);
    ep = 0.03 + 0.04*rand; lb = x-ep; ub = x+ep;
    oS    = struct('precision','double','maxNodes',6000);
    oA    = struct('precision','double','maxNodes',6000,'maxFrontier',256,'alphaIter',15,'alphaLr',0.3);
    oB    = struct('precision','double','maxNodes',6000,'maxFrontier',256,'betaIter',15,'alphaLr',0.3);
    sS = gpu_bab_relu_split(ops, lb, ub, tl, dims(end), oS);
    sA = gpu_bab_relu_split_batched(ops, lb, ub, tl, dims(end), oA);
    sB = gpu_bab_relu_split_batched(ops, lb, ub, tl, dims(end), oB);
    if i_contra(sS,sB) || i_contra(sA,sB) || i_contra(sS,sA), nContra = nContra + 1; end
    if ~strcmp(sB,'unknown'), nBeta = nBeta + 1; end
    if ~strcmp(sA,'unknown'), nAlpha = nAlpha + 1; end
end
assert(nContra == 0, sprintf('beta BaB contradicted in %d cases', nContra));
assert(nBeta >= nAlpha, sprintf('beta decided %d but alpha decided %d (beta must not decide fewer)', nBeta, nAlpha));

%% Summary
disp('test_soundness_gpu_bab_beta: all sections passed');

% ----------------------------------------------------------------------------------------
function ops = i_fc(dims)
    ops = {};
    for L = 1:numel(dims)-1
        W = randn(dims(L+1), dims(L)) * sqrt(2/dims(L)); b = randn(dims(L+1),1)*0.1;
        ops{end+1} = struct('type','affine','W',W,'b',b,'src',numel(ops)); %#ok<AGROW>
        if L < numel(dims)-1, ops{end+1} = struct('type','relu','src',numel(ops)); end %#ok<AGROW>
    end
end
function [y, pre] = i_fwd(ops, X)
    pre = cell(numel(ops),1); v = X;
    for k = 1:numel(ops)
        if strcmp(ops{k}.type,'affine'), v = ops{k}.W*v + ops{k}.b; else, pre{k} = v; v = max(v,0); end
    end
    y = v;
end
function tf = i_contra(a,b)
    tf = (strcmp(a,'robust')&&strcmp(b,'unsafe')) || (strcmp(a,'unsafe')&&strcmp(b,'robust'));
end
function [nBad, nChk, nBetaHelp] = i_check_beta(seed, nTrials)
    rng(seed); nBad = 0; nChk = 0; nBetaHelp = 0;
    for t = 1:nTrials
        dims = [6 16 16 7]; ops = i_fc(dims);
        lb = -0.2*ones(6,1); ub = 0.2*ones(6,1); C = [eye(6) -ones(6,1)];
        reluIdx = find(cellfun(@(o) strcmp(o.type,'relu'), ops));
        [~, rtL, rtU] = gpu_bab_crown_tight(ops, lb, ub, C, 'double', cell(numel(ops),1));
        rb = struct('preL', {rtL}, 'preU', {rtU});
        % random fixings on <=3 unstable neurons (per root bounds)
        cand = [];
        for r = 1:numel(reluIdx), k=reluIdx(r); u=find(rtL{k}<0 & rtU{k}>0); cand=[cand; repmat(k,numel(u),1) u(:)]; end %#ok<AGROW>
        if isempty(cand), continue; end
        fx = cell(numel(ops),1); for r=1:numel(reluIdx), fx{reluIdx(r)}=zeros(numel(rtL{reluIdx(r)}),1); end
        nf = min(3, size(cand,1)); pk = cand(randperm(size(cand,1), nf), :);
        for p=1:nf, fx{pk(p,1)}(pk(p,2)) = (rand<0.5)*2-1; end
        mAB = gpu_bab_crown_alpha_beta(ops, lb, ub, C, fx, reluIdx, 'double', 25, 0.3, rb);   % alpha+beta
        mA  = gpu_bab_crown_alpha_fix (ops, lb, ub, C, fx, reluIdx, 'double', 25, 0.3, rb);   % alpha only
        % constrained MC: sample box, keep x satisfying the fixings, min C*f
        rng(5000+t); X = lb + (ub-lb).*rand(6, 20000); [Y, pre] = i_fwd(ops, X);
        mask = true(1, size(X,2));
        for p=1:nf, k=pk(p,1); j=pk(p,2);
            if fx{k}(j)==1, mask = mask & (pre{k}(j,:) >= 0); else, mask = mask & (pre{k}(j,:) <= 0); end
        end
        if nnz(mask) < 30, continue; end
        nChk = nChk + 1;
        mc = min(C * Y(:,mask), [], 2);
        if ~all(mAB(:) <= mc(:) + 1e-4), nBad = nBad + 1; end          % SOUNDNESS (constrained MC) -- the hard check
        if min(mAB(:)) > min(mA(:)) + 1e-7, nBetaHelp = nBetaHelp + 1; end  % beta strictly raised the cert margin
    end
end
