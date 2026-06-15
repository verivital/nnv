% test_soundness_gpu_bab_alpha
% alpha-CROWN node bounding (gpu_bab_crown_alpha_fix) with ROOT-TIGHT reuse must be SOUND and
% at least as tight as plain root-tight: it reuses the tight root intermediate bounds and
% optimizes the unstable lower slopes (alpha in [0,1]) to MAXIMIZE the spec lower bound, keeping
% the best iterate (init = the min-area slope = the root-tight bound). FC+ReLU, double precision.

%% Test 1: alpha with 0 iters + rootBounds == root-tight (gpu_bab_crown_spec_dag) margin exactly
rng(1); ops = i_fc([5 12 12 6]);
lb = -0.2*ones(5,1); ub = 0.2*ones(5,1); C = [eye(5) -ones(5,1)];
reluIdx = find(cellfun(@(o) strcmp(o.type,'relu'), ops));
[~, rtL, rtU] = gpu_bab_crown_tight(ops, lb, ub, C, 'double', cell(numel(ops),1));
rb = struct('preL', {rtL}, 'preU', {rtU});
mTight = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double', cell(numel(ops),1), rb);
mA0 = gpu_bab_crown_alpha_fix(ops, lb, ub, C, cell(numel(ops),1), reluIdx, 'double', 0, 0.2, rb);
assert(max(abs(mTight(:) - mA0(:))) < 1e-9, ...
    sprintf('alpha(0 iter)+root must equal root-tight margin; max diff %.3e', max(abs(mTight(:)-mA0(:)))));

%% Test 2: alpha (iters>0) is >= root-tight (optimization only raises the bound) AND sound
rng(2); ops = i_fc([6 16 16 7]);
lb = -0.15*ones(6,1); ub = 0.15*ones(6,1); C = [eye(6) -ones(6,1)];
reluIdx = find(cellfun(@(o) strcmp(o.type,'relu'), ops));
[~, rtL, rtU] = gpu_bab_crown_tight(ops, lb, ub, C, 'double', cell(numel(ops),1));
rb = struct('preL', {rtL}, 'preU', {rtU});
mTight = gpu_bab_crown_spec_dag(ops, lb, ub, C, 'double', cell(numel(ops),1), rb);
mAlpha = gpu_bab_crown_alpha_fix(ops, lb, ub, C, cell(numel(ops),1), reluIdx, 'double', 30, 0.3, rb);
% alpha maximizes the per-column MIN-margin (the BaB certification quantity, init = root-tight),
% so that min must not drop -- but individual spec margins MAY trade off (aggregate objective).
assert(min(mAlpha(:)) >= min(mTight(:)) - 1e-6, ...
    sprintf('alpha min-margin must be >= root-tight (%.4f vs %.4f)', min(mAlpha(:)), min(mTight(:))));
mc = i_mc_min(ops, C, lb, ub, 6000);
assert(all(mAlpha(:) <= mc(:) + 1e-5), 'alpha must be a SOUND lower bound (<= Monte-Carlo min)');

%% Test 3: BaB with alpha never contradicts serial, and decides >= the root-tight-only BaB
rng(3); nContra = 0; nAlpha = 0; nRoot = 0;
for t = 1:8
    dims = [6 18 18 6]; ops = i_fc(dims);
    x = randn(6,1)*0.5; yc = i_fwd(ops, x); [~, tl] = max(yc);
    ep = 0.03 + 0.04*rand; lb = x-ep; ub = x+ep;
    oS = struct('precision','double','maxNodes',8000);
    oRoot  = struct('precision','double','maxNodes',8000,'maxFrontier',256,'rootTight',true,'alphaIter',0);
    oAlpha = struct('precision','double','maxNodes',8000,'maxFrontier',256,'rootTight',true,'alphaIter',20,'alphaLr',0.3);
    sS  = gpu_bab_relu_split(ops, lb, ub, tl, dims(end), oS);
    sR  = gpu_bab_relu_split_batched(ops, lb, ub, tl, dims(end), oRoot);
    sA  = gpu_bab_relu_split_batched(ops, lb, ub, tl, dims(end), oAlpha);
    if i_contra(sS,sA) || i_contra(sR,sA) || i_contra(sS,sR), nContra = nContra + 1; end
    if ~strcmp(sA,'unknown'), nAlpha = nAlpha + 1; end
    if ~strcmp(sR,'unknown'), nRoot  = nRoot  + 1; end
end
assert(nContra == 0, sprintf('alpha BaB contradicted in %d cases', nContra));
assert(nAlpha >= nRoot, sprintf('alpha decided %d but root-tight decided %d (alpha must not decide fewer)', nAlpha, nRoot));

%% Summary
disp('test_soundness_gpu_bab_alpha: all sections passed');

% ----------------------------------------------------------------------------------------
function ops = i_fc(dims)
    ops = {};
    for L = 1:numel(dims)-1
        W = randn(dims(L+1), dims(L)) * sqrt(2/dims(L)); b = randn(dims(L+1),1)*0.1;
        ops{end+1} = struct('type','affine','W',W,'b',b,'src',numel(ops)); %#ok<AGROW>
        if L < numel(dims)-1, ops{end+1} = struct('type','relu','src',numel(ops)); end %#ok<AGROW>
    end
end
function y = i_fwd(ops, X)
    v = X;
    for k = 1:numel(ops)
        if strcmp(ops{k}.type,'affine'), v = ops{k}.W*v + ops{k}.b; else, v = max(v,0); end
    end
    y = v;
end
function mn = i_mc_min(ops, C, lb, ub, nS)
    rng(101); X = lb + (ub-lb).*rand(numel(lb), nS); mn = min(C * i_fwd(ops, X), [], 2);
end
function tf = i_contra(a,b)
    tf = (strcmp(a,'robust')&&strcmp(b,'unsafe')) || (strcmp(a,'unsafe')&&strcmp(b,'robust'));
end
