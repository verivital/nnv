function [status, info] = gpu_bab_relu_split_batched(ops, x_lb, x_ub, trueLabel, nClasses, opts)
% GPU_BAB_RELU_SPLIT_BATCHED  Batched-DFS ReLU-split branch-and-bound for FC+ReLU nets --
%   the GPU-saturating refinement of gpu_bab_relu_split. The serial version bounds ONE BaB
%   node per CROWN pass, so on the GPU the device is idle at batch=1. This keeps a STACK of
%   BaB nodes and, each round, pops a BATCH of up to maxFrontier of them and bounds them in
%   ONE gpu_bab_crown_spec call -- the B columns are nodes that share the input box but
%   partition the unstable neurons via per-relu fixing clamps. A LIFO stack keeps the search
%   depth-first (bounded live set + fast early termination on the convex barrier, like the
%   serial DFS) while the GPU still processes B nodes per kernel.
%
%   [status, info] = GPU_BAB_RELU_SPLIT_BATCHED(ops, x_lb, x_ub, trueLabel, nClasses, opts)
%     status : 'robust' | 'unsafe' (info.cex misclassifies) | 'unknown' (budget/barrier)
%     opts (all optional):
%       .precision   'single'(default)|'double'
%       .maxNodes    20000  total BaB nodes bounded before giving up (unknown)
%       .maxFrontier 512    nodes bounded per batched GPU call (the GPU batch width)
%       .maxStack    200000 live-stack cap (memory guard -> unknown if exceeded)
%       .margin      0      FP-soundness slack: require every spec margin > margin
%       .nSample     16     random samples per round for the concrete cex search
%
%   SOUNDNESS (sound-or-unknown; identical guarantees to gpu_bab_relu_split's IBP-clamped
%   path -- the bounding is bound-for-bound the same, just batched over node columns):
%     * each node is the input box plus a set of ReLU fixings; a split's active/inactive
%       children PARTITION that neuron (z>=0 or z<0), so 'robust' is returned only when the
%       stack empties with EVERY popped node certified (all spec margins > margin);
%     * clamping l up (active) / u down (inactive) is a sound, tighter bound on the child's
%       sub-domain, so each node's CROWN margin is a valid lower bound there;
%     * 'unsafe' only from a CONCRETE misclassifying input on the ORIGINAL box;
%     * 'unknown' on the node/stack budget, or when a popped node is undecided yet has no
%       unstable unfixed neuron left to split (the convex-relaxation barrier).
%   FC+ReLU only: conv/pool/normaffine need the tight intermediate bounds of
%   gpu_bab_relu_split (intermediate='tight'); the IBP-clamped bounding here would be
%   unsound for them, so we refuse loudly.

    if nargin < 6, opts = struct(); end
    precision   = i_get(opts, 'precision', 'single');
    maxNodes    = i_get(opts, 'maxNodes', 20000);
    maxFrontier = i_get(opts, 'maxFrontier', 512);
    maxStack    = i_get(opts, 'maxStack', 200000);
    margin      = cast(i_get(opts, 'margin', 0), precision);
    nSample     = i_get(opts, 'nSample', 16);

    % FC+ReLU guard: the batched bounding uses IBP intermediate bounds, sound ONLY for
    % affine+relu (a conv/pool/normaffine op would hit gpu_bab_crown_spec's error).
    if any(cellfun(@(o) ~any(strcmp(o.type, {'affine','relu'})), ops))
        error('gpu_bab_relu_split_batched:fcOnly', ...
            'affine+relu only; use gpu_bab_relu_split (intermediate=''tight'') for conv/pool nets.');
    end

    nOps    = numel(ops);
    reluIdx = find(cellfun(@(o) strcmp(o.type, 'relu'), ops));
    C = -eye(nClasses, precision);
    C(:, trueLabel) = C(:, trueLabel) + 1;
    C(trueLabel, :) = [];
    info = struct('nodes', 0, 'rounds', 0, 'maxStack', 1, 'cex', [], 'cexLabel', []);

    % concrete counterexample on the whole box first (cheap falsification)
    [cex, lab] = i_find_cex(ops, x_lb, x_ub, trueLabel, nSample, precision);
    if ~isempty(cex)
        status = 'unsafe'; info.cex = cex; info.cexLabel = lab; return;
    end

    % --- root: bound the un-split network once; certify, or seed the stack from its split ---
    [mRoot, preL, preU] = i_bound_batch(ops, reluIdx, x_lb, x_ub, C, precision, cell(nOps,1), 1);
    info.nodes = 1;
    if all(mRoot > margin, 1)
        status = 'robust'; return;
    end
    reluDim = zeros(1, nOps);
    for r = 1:numel(reluIdx), reluDim(reluIdx(r)) = size(preL{reluIdx(r)}, 1); end
    [sK, sJ, hasS] = i_pick_splits(reluIdx, preL, preU, cell(nOps,1), 1);
    if ~hasS
        status = 'unknown'; return;                 % root undecided, nothing to split
    end
    % stack as per-relu fixing matrices fixStack{k} = [dim_k x M]; seed with root's 2 children.
    fixStack = cell(nOps, 1);
    for r = 1:numel(reluIdx), k = reluIdx(r); fixStack{k} = zeros(reluDim(k), 2); end
    fixStack{sK}(sJ, 1) =  1;                        % active child
    fixStack{sK}(sJ, 2) = -1;                        % inactive child
    M = 2;

    % --- batched DFS over the node stack ---
    while M > 0
        info.rounds = info.rounds + 1;
        info.maxStack = max(info.maxStack, M);
        if M > maxStack
            status = 'unknown'; return;             % live set too large (memory guard)
        end
        B = min(maxFrontier, M);
        cols = (M - B + 1):M;                       % pop the top B nodes (LIFO -> DFS)
        fixc = cell(nOps, 1);
        for r = 1:numel(reluIdx), k = reluIdx(r); fixc{k} = fixStack{k}(:, cols); end
        for r = 1:numel(reluIdx), k = reluIdx(r); fixStack{k}(:, cols) = []; end
        M = M - B;
        info.nodes = info.nodes + B;
        if info.nodes > maxNodes
            status = 'unknown'; return;
        end

        % bound the popped batch
        [margins, preL, preU, infeas] = i_bound_batch(ops, reluIdx, x_lb, x_ub, C, precision, fixc, B);
        % An INFEASIBLE node (a fixing combination that made the clamped IBP box empty,
        % preL>preU at some neuron) represents an EMPTY sub-region: the property holds
        % vacuously, so it is certified. Without this, IBP+CROWN cannot bound such a node,
        % so the BaB splits it forever and degrades a robust query to a false 'unknown'.
        undec = find(~(all(margins > margin, 1) | infeas));
        if isempty(undec)
            continue;                               % whole batch certified; pop the next
        end

        % periodic concrete counterexample search (cheap; once per round on the original box)
        [cex, lab] = i_find_cex(ops, x_lb, x_ub, trueLabel, nSample, precision);
        if ~isempty(cex)
            status = 'unsafe'; info.cex = cex; info.cexLabel = lab; return;
        end

        % split each undecided node on its largest-gap unstable unfixed neuron
        [sK, sJ, hasS] = i_pick_splits(reluIdx, preL, preU, fixc, undec);
        if ~all(hasS)
            status = 'unknown'; return;             % an undecided node cannot be split
        end

        % push the two children of every undecided node onto the stack
        nu = numel(undec);
        for r = 1:numel(reluIdx)
            k = reluIdx(r);
            parent = fixc{k}(:, undec);             % dim_k x nu
            fixStack{k} = [fixStack{k}, parent, parent];   % [active block | inactive block]
        end
        for i = 1:nu
            k = sK(i); j = sJ(i);
            fixStack{k}(j, M + i)      =  1;        % active child  (column M+i)
            fixStack{k}(j, M + nu + i) = -1;        % inactive child (column M+nu+i)
        end
        M = M + 2 * nu;
    end
    status = 'robust';
end

% ------------------------------------------------------------------------------------
function [margins, preL, preU, infeasible] = i_bound_batch(ops, reluIdx, x_lb, x_ub, C, precision, fixc, B)
% Bound B frontier nodes (B <= maxFrontier) in one batched gpu_bab_crown_spec call. The
% backward CROWN (pagemtimes) runs on-device; margins and the per-relu pre-activation
% bounds are gathered to host (small) for the verdict + split-neuron selection.
% infeasible(j) is true when node j's fixings made the clamped IBP box empty (preL>preU at
% some neuron) -- an empty sub-region the caller certifies vacuously.
    LB = repmat(x_lb, 1, B); UB = repmat(x_ub, 1, B);
    [m, pL, pU] = gpu_bab_crown_spec(ops, LB, UB, C, precision, fixc);
    margins = gather(m);
    preL = cell(numel(ops), 1); preU = cell(numel(ops), 1);
    infeasible = false(1, B);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        preL{k} = gather(pL{k});
        preU{k} = gather(pU{k});
        infeasible = infeasible | any(preL{k} > preU{k} + 1e-6, 1);   % clamp made box empty
    end
end

% ------------------------------------------------------------------------------------
function [splitK, splitJ, hasSplit] = i_pick_splits(reluIdx, preL, preU, fixc, undec)
% For each undecided node, pick the unstable, unfixed neuron with the largest ReLU
% relaxation gap bu = u*(-l)/(u-l) (the upper-line slack a split removes) -- a BaBSR-lite
% heuristic that closes the bound far faster than first-unstable. Returns the chosen relu
% op index + neuron index per node, and a flag (false => no splittable neuron left).
    nu = numel(undec);
    bestGap = -inf(1, nu);
    splitK  = zeros(1, nu);
    splitJ  = zeros(1, nu);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        l = preL{k}(:, undec); u = preU{k}(:, undec);     % dim_k x nu
        uns = (l < 0) & (u > 0);
        if isempty(fixc{k}), free = true(size(l)); else, free = (fixc{k}(:, undec) == 0); end
        gap = u .* (-l) ./ (u - l);
        gap(~(uns & free)) = -inf;
        [g, j] = max(gap, [], 1);                          % 1 x nu
        better = g > bestGap;
        bestGap(better) = g(better);
        splitK(better)  = k;
        splitJ(better)  = j(better);
    end
    hasSplit = isfinite(bestGap);
end

% ------------------------------------------------------------------------------------
function [cex, lab] = i_find_cex(ops, LB, UB, trueLabel, nSample, precision)
% Exact forward eval at the box center + random in-box samples; return the first input that
% misclassifies (a genuine witness -> sound 'unsafe'). 'like' LB keeps the samples on the
% same device as the input (no host<->GPU transfer for the GPU path).
    cex = []; lab = [];
    X = (LB + UB) / 2;
    for s = 1:max(0, nSample)
        X = [X, LB + (UB - LB) .* rand(size(LB), 'like', LB)]; %#ok<AGROW>
    end
    Y = gpu_bab_ibp(ops, X, X, precision);
    [~, pred] = max(Y, [], 1);
    bad = find(pred ~= trueLabel, 1);
    if ~isempty(bad), cex = X(:, bad); lab = pred(bad); end
end

function v = i_get(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
