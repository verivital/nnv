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
%       .rootTight   true   compute the TIGHT intermediate bounds ONCE at the root
%                           (gpu_bab_crown_tight) and reuse them -- clamped per node -- for
%                           every batched bound, instead of the loose per-node IBP forward.
%                           Tight bounds at batched speed; set false for the old IBP path.
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
    rootTight   = i_get(opts, 'rootTight', true);   % root-tight intermediate-bound reuse (#1 tightness lever)
    alphaIter   = i_get(opts, 'alphaIter', 0);      % alpha-CROWN slope-opt iters per batched node bound (0 = off; FC only)
    alphaLr     = i_get(opts, 'alphaLr', 0.2);
    betaIter    = i_get(opts, 'betaIter', 0);       % beta-CROWN split-dual iters (>0 -> JOINT alpha+beta; FC only)

    % Supported ops: affine/relu (FC) + conv/normaffine/avgpool/add (SEQUENTIAL + residual-DAG
    % conv nets). The batched bounding (gpu_bab_crown_spec_dag) is sound for these -- conv/BN/
    % avgpool/add are LINEAR (exact interval forward + exact adjoint backward; 'add' routes its
    % backward coefficient unchanged to both summands, gpu_bab_crown_spec_dag handles the DAG via
    % op.src/op.inputs), only ReLU is relaxed. 'maxpool' needs a per-node window relaxation not
    % yet batched -> refuse (sound-by-refusal -> caller runs the serial tight path / Star).
    supported = {'affine','relu','conv','normaffine','avgpool','add'};
    badIdx = find(cellfun(@(o) ~any(strcmp(o.type, supported)), ops), 1);
    if ~isempty(badIdx)
        error('gpu_bab_relu_split_batched:unsupportedOp', ...
            'op "%s" unsupported (affine/relu/conv/normaffine/avgpool/add only; maxpool needs gpu_bab_relu_split tight).', ops{badIdx}.type);
    end

    nOps    = numel(ops);
    reluIdx = find(cellfun(@(o) strcmp(o.type, 'relu'), ops));
    % alpha/beta-CROWN node bounding (gpu_bab_crown_alpha_fix / _alpha_beta) is FC-only
    % (affine/relu); for conv/normaffine/avgpool nets fall back to the (root-tight) CROWN bound
    % so we never mis-bound.
    % alpha/beta route: FC (affine/relu) -> alpha_fix/alpha_beta; DAG (conv/normaffine/avgpool/add)
    % -> gpu_bab_crown_alpha_dag (alpha+beta over the full DAG). Only maxpool has no batched alpha
    % backward -> disable alpha/beta there (sound: falls back to fixed-slope spec_dag / tight path).
    if (alphaIter > 0 || betaIter > 0) && any(cellfun(@(o) strcmp(o.type, 'maxpool'), ops))
        alphaIter = 0; betaIter = 0;
    end
    % SPEC: argmax-robustness (default) builds C internally; a general-halfspace caller passes
    % opts.spec = struct('C',C,'g',g,'rowGroups',{rows-per-disjunct}). The node certification then
    % generalises "all margins > 0" to "every unsafe DISJUNCT has SOME row G_i*y-g_i > 0" (the
    % output avoids every unsafe polytope). argmax is the special case: each margin is its own
    % single-row disjunct, so the disjunctive test reduces EXACTLY to all-margins>0.
    spec = i_get(opts, 'spec', []);
    if isempty(spec)
        C = -eye(nClasses, precision);
        C(:, trueLabel) = C(:, trueLabel) + 1;
        C(trueLabel, :) = [];
        gOff = zeros(size(C,1), 1, precision);
        rowGroups = num2cell(1:size(C,1));
    else
        C = cast(spec.C, precision); gOff = cast(spec.g(:), precision); rowGroups = spec.rowGroups;
    end
    doCex = isempty(spec);             % i_find_cex is argmax (misclassification) only; skip for halfspace
    info = struct('nodes', 0, 'rounds', 0, 'maxStack', 1, 'cex', [], 'cexLabel', []);

    % concrete counterexample on the whole box first (cheap falsification; argmax only)
    if doCex
        [cex, lab] = i_find_cex(ops, x_lb, x_ub, trueLabel, nSample, precision);
        if ~isempty(cex)
            status = 'unsafe'; info.cex = cex; info.cexLabel = lab; return;
        end
    end

    % --- root: bound the un-split network once; certify, or seed the stack from its split.
    % With rootTight (default), compute the TIGHT intermediate bounds ONCE here via a backward
    % CROWN pass over the full box (gpu_bab_crown_tight) and reuse them -- clamped per node --
    % for every batched bound below. This gives tight bounds (vs the loose per-node IBP) at
    % batched speed (vs recomputing the tight pass per node). Falls back to the IBP forward
    % when rootTight is off. ---
    rootBounds = [];
    if rootTight
        [mRoot, rtL, rtU] = gpu_bab_crown_tight(ops, x_lb, x_ub, C, precision, cell(nOps,1));
        mRoot = gather(mRoot(:));
        rootBounds = struct('preL', {rtL}, 'preU', {rtU});
        preL = cell(nOps,1); preU = cell(nOps,1);
        for r = 1:numel(reluIdx), k = reluIdx(r); preL{k} = gather(rtL{k}); preU{k} = gather(rtU{k}); end
    else
        [mRoot, preL, preU] = i_bound_batch(ops, reluIdx, x_lb, x_ub, C, precision, cell(nOps,1), 1);
    end
    info.nodes = 1;
    if i_certify_dis(mRoot, gOff, rowGroups, margin)
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

        % bound the popped batch (reusing the tight root bounds, clamped per node, if rootTight)
        [margins, preL, preU, infeas] = i_bound_batch(ops, reluIdx, x_lb, x_ub, C, precision, fixc, B, rootBounds, alphaIter, alphaLr, betaIter);
        % An INFEASIBLE node (a fixing combination that made the clamped IBP box empty,
        % preL>preU at some neuron) represents an EMPTY sub-region: the property holds
        % vacuously, so it is certified. Without this, IBP+CROWN cannot bound such a node,
        % so the BaB splits it forever and degrades a robust query to a false 'unknown'.
        undec = find(~(i_certify_dis(margins, gOff, rowGroups, margin) | infeas));
        if isempty(undec)
            continue;                               % whole batch certified; pop the next
        end

        % periodic concrete counterexample search (cheap; argmax only -- halfspace sat is the dispatcher's job)
        if doCex
            [cex, lab] = i_find_cex(ops, x_lb, x_ub, trueLabel, nSample, precision);
            if ~isempty(cex)
                status = 'unsafe'; info.cex = cex; info.cexLabel = lab; return;
            end
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
function ok = i_certify_dis(margins, gOff, rowGroups, marginSlack)
% Disjunctive certification: a node is certified SAFE iff EVERY unsafe disjunct is avoided, and a
% disjunct {G*y<=g} is avoided iff SOME of its rows has a lower bound G_i*y-g_i > slack (the
% reachable set lies entirely outside that halfspace -> outside the polytope). margins: nRows x B
% (sound lower bound of G*y). Returns 1 x B. For argmax (each row its own disjunct, g=0) this is
% exactly all(margins>slack) -- the original argmax certification, bound-for-bound.
    av = (margins - gOff) > marginSlack;          % nRows x B : row i provably avoided
    ok = true(1, size(margins, 2));
    for d = 1:numel(rowGroups)
        ok = ok & any(av(rowGroups{d}, :), 1);    % disjunct d avoided iff some of its rows avoided
    end
end

% ------------------------------------------------------------------------------------
function [margins, preL, preU, infeasible] = i_bound_batch(ops, reluIdx, x_lb, x_ub, C, precision, fixc, B, rootBounds, alphaIter, alphaLr, betaIter)
% Bound B frontier nodes (B <= maxFrontier) in one batched CROWN call. The backward CROWN
% (pagemtimes) runs on-device; margins and the per-relu pre-activation bounds are gathered to
% host (small) for the verdict + split-neuron selection. infeasible(j) is true when node j's
% fixings made the clamped box empty (preL>preU at some neuron) -- an empty sub-region the
% caller certifies vacuously. rootBounds (optional) routes the tight root bounds in (root-tight
% reuse). Node-bound tightening (FC nets only), strongest first:
%   betaIter>0  -> gpu_bab_crown_alpha_beta (JOINT alpha-CROWN + beta-CROWN split duals)
%   alphaIter>0 -> gpu_bab_crown_alpha_fix  (alpha-CROWN slopes only)
%   else        -> gpu_bab_crown_spec_dag   (fixed-slope root-tight CROWN)
% All sound (alpha in [0,1], beta>=0 are valid relaxations/Lagrangian duals).
    if nargin < 9, rootBounds = []; end
    if nargin < 10 || isempty(alphaIter), alphaIter = 0; end
    if nargin < 11 || isempty(alphaLr), alphaLr = 0.2; end
    if nargin < 12 || isempty(betaIter), betaIter = 0; end
    LB = repmat(x_lb, 1, B); UB = repmat(x_ub, 1, B);
    % DAG nets (conv/normaffine/avgpool/add) route to gpu_bab_crown_alpha_dag (alpha+beta over the
    % full DAG); the FC-only alpha_fix/alpha_beta mis-bound a conv/add op (fail-open) so are used
    % ONLY for pure affine/relu. alphaIter==betaIter==0 -> fixed-slope spec_dag (the cheap bound).
    isDAG = any(cellfun(@(o) any(strcmp(o.type, {'conv','normaffine','avgpool','add'})), ops));
    if isDAG
        if alphaIter > 0 || betaIter > 0
            [m, ~, pL, pU] = gpu_bab_crown_alpha_dag(ops, LB, UB, C, fixc, reluIdx, precision, max(alphaIter,betaIter), alphaLr, rootBounds);
        else
            [m, pL, pU] = gpu_bab_crown_spec_dag(ops, LB, UB, C, precision, fixc, rootBounds);
        end
    elseif betaIter > 0
        [m, ~, pL, pU] = gpu_bab_crown_alpha_beta(ops, LB, UB, C, fixc, reluIdx, precision, max(alphaIter,betaIter), alphaLr, rootBounds);
    elseif alphaIter > 0
        [m, ~, pL, pU] = gpu_bab_crown_alpha_fix(ops, LB, UB, C, fixc, reluIdx, precision, alphaIter, alphaLr, rootBounds);
    else
        [m, pL, pU] = gpu_bab_crown_spec_dag(ops, LB, UB, C, precision, fixc, rootBounds);
    end
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
