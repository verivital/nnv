function [status, info] = gpu_bab_relu_split_batched(ops, x_lb, x_ub, trueLabel, nClasses, opts)
% GPU_BAB_RELU_SPLIT_BATCHED  Batched-frontier ReLU-split branch-and-bound for FC+ReLU
%   nets -- the GPU-saturating refinement of gpu_bab_relu_split. The serial version is a
%   DFS that bounds ONE BaB node per CROWN pass, so on the GPU the device sits idle at
%   batch=1. This version keeps the WHOLE frontier of BaB nodes and bounds them in ONE
%   batched gpu_bab_crown_spec call: the B columns are nodes that share the SAME input box
%   but partition the unstable neurons via per-relu fixing clamps. Throughput scales with
%   the frontier width (pagemtimes over the node dimension), not serial node count.
%
%   [status, info] = GPU_BAB_RELU_SPLIT_BATCHED(ops, x_lb, x_ub, trueLabel, nClasses, opts)
%     status : 'robust' | 'unsafe' (info.cex misclassifies) | 'unknown' (budget/barrier)
%     opts (all optional):
%       .precision   'single'(default)|'double'
%       .maxNodes    4096   total BaB nodes explored before giving up (unknown)
%       .maxFrontier 512    GPU frontier-width cap; wider frontiers are bounded in chunks
%       .margin      0      FP-soundness slack: require every margin > margin
%       .nSample     16     random samples per round for the concrete cex search
%
%   SOUNDNESS (sound-or-unknown; identical guarantees to gpu_bab_relu_split's IBP-clamped
%   path -- the bounding is bound-for-bound the same, just batched over node columns):
%     * each column is the input box plus a set of ReLU fixings; a split's active/inactive
%       children PARTITION that neuron (z>=0 or z<0), so 'robust' is returned only when
%       EVERY surviving leaf is certified (all spec margins > margin);
%     * clamping l up (active) / u down (inactive) is a sound, tighter bound on the child's
%       sub-domain, so each leaf's CROWN margin is a valid lower bound there;
%     * 'unsafe' only from a CONCRETE misclassifying input on the ORIGINAL box (a real
%       witness for the whole query), never inferred from a bound;
%     * 'unknown' on the node budget, or when an undecided node has no unstable unfixed
%       neuron left to split (the convex-relaxation barrier -- needs tighter bounds).
%   FC+ReLU only: a conv/pool/normaffine op needs the tight intermediate bounds of
%   gpu_bab_relu_split (intermediate='tight'); the IBP-clamped bounding here would be
%   unsound for them (it would silently treat them as a ReLU), so we refuse loudly.

    if nargin < 6, opts = struct(); end
    precision   = i_get(opts, 'precision', 'single');
    maxNodes    = i_get(opts, 'maxNodes', 4096);
    maxFrontier = i_get(opts, 'maxFrontier', 512);
    margin      = cast(i_get(opts, 'margin', 0), precision);
    nSample     = i_get(opts, 'nSample', 16);

    % FC+ReLU guard: the batched bounding uses IBP intermediate bounds, sound ONLY for
    % affine+relu. A conv/avgpool/maxpool/normaffine op would hit gpu_bab_crown_spec's
    % unsupported-op error -- refuse here with a clear message instead.
    if any(cellfun(@(o) ~any(strcmp(o.type, {'affine','relu'})), ops))
        error('gpu_bab_relu_split_batched:fcOnly', ...
            'affine+relu only; use gpu_bab_relu_split (intermediate=''tight'') for conv/pool nets.');
    end

    nOps    = numel(ops);
    reluIdx = find(cellfun(@(o) strcmp(o.type, 'relu'), ops));
    C = -eye(nClasses, precision);
    C(:, trueLabel) = C(:, trueLabel) + 1;
    C(trueLabel, :) = [];
    info = struct('nodes', 0, 'rounds', 0, 'maxFrontier', 1, 'cex', [], 'cexLabel', []);

    % concrete counterexample on the whole box first (cheap falsification)
    [cex, lab] = i_find_cex(ops, x_lb, x_ub, trueLabel, nSample, precision);
    if ~isempty(cex)
        status = 'unsafe'; info.cex = cex; info.cexLabel = lab; return;
    end

    % Frontier: per-relu fixing matrix fix{k} = [dim_k x Q] of -1/0/+1 (host: used only for
    % logical indexing + child construction). Start with one all-free node (Q=1).
    fix = cell(nOps, 1);
    Q = 1;

    while Q > 0
        info.rounds = info.rounds + 1;
        info.nodes  = info.nodes + Q;
        info.maxFrontier = max(info.maxFrontier, Q);
        if info.nodes > maxNodes
            status = 'unknown'; return;
        end

        % --- batched bounding of the whole frontier (chunked to maxFrontier columns) ---
        [margins, preL, preU] = i_bound_frontier(ops, reluIdx, x_lb, x_ub, C, precision, fix, Q, maxFrontier);
        certified = all(margins > margin, 1);              % 1 x Q (host)
        undec = find(~certified);
        if isempty(undec)
            status = 'robust'; return;
        end

        % --- periodic concrete counterexample search (once per round, original box) ---
        [cex, lab] = i_find_cex(ops, x_lb, x_ub, trueLabel, nSample, precision);
        if ~isempty(cex)
            status = 'unsafe'; info.cex = cex; info.cexLabel = lab; return;
        end

        % --- per-node split: largest ReLU relaxation gap among unstable, unfixed neurons ---
        [splitK, splitJ, hasSplit] = i_pick_splits(reluIdx, preL, preU, fix, undec);
        if ~all(hasSplit)
            status = 'unknown'; return;                    % an undecided node cannot be split
        end

        % --- new frontier: each undecided node -> active child (+1) and inactive child (-1) ---
        fix = i_make_children(reluIdx, fix, preL, undec, splitK, splitJ);
        Q = 2 * numel(undec);
    end
    status = 'robust';
end

% ------------------------------------------------------------------------------------
function [margins, preL, preU] = i_bound_frontier(ops, reluIdx, x_lb, x_ub, C, precision, fix, Q, maxFrontier)
% Bound all Q frontier nodes, chunking the columns to at most maxFrontier per batched call
% (bounds GPU memory). The backward CROWN (pagemtimes) runs on-device inside
% gpu_bab_crown_spec; margins and the per-relu pre-activation bounds are gathered to host
% (small) for the verdict decision and split-neuron selection.
    nOps  = numel(ops);
    nSpec = size(C, 1);
    margins = zeros(nSpec, Q, precision);
    preL = cell(nOps, 1); preU = cell(nOps, 1);
    for c0 = 1:maxFrontier:Q
        cols = c0:min(c0 + maxFrontier - 1, Q);
        nb = numel(cols);
        fixc = cell(nOps, 1);
        for r = 1:numel(reluIdx)
            k = reluIdx(r);
            if ~isempty(fix{k}), fixc{k} = fix{k}(:, cols); end
        end
        LB = repmat(x_lb, 1, nb); UB = repmat(x_ub, 1, nb);
        [m, pL, pU] = gpu_bab_crown_spec(ops, LB, UB, C, precision, fixc);
        margins(:, cols) = gather(m);
        for r = 1:numel(reluIdx)
            k = reluIdx(r);
            if c0 == 1
                preL{k} = zeros(size(pL{k}, 1), Q, precision);
                preU{k} = zeros(size(pU{k}, 1), Q, precision);
            end
            preL{k}(:, cols) = gather(pL{k});
            preU{k}(:, cols) = gather(pU{k});
        end
    end
end

% ------------------------------------------------------------------------------------
function [splitK, splitJ, hasSplit] = i_pick_splits(reluIdx, preL, preU, fix, undec)
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
        if isempty(fix{k}), free = true(size(l)); else, free = (fix{k}(:, undec) == 0); end
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
function fix = i_make_children(reluIdx, fix, preL, undec, splitK, splitJ)
% Build the next frontier: for each undecided node, two children (active +1 / inactive -1
% on its split neuron). Columns 1..nu are the active children, nu+1..2nu the inactive
% children, in the same node order as `undec`.
    nu   = numel(undec);
    nOps = numel(fix);
    newfix = cell(nOps, 1);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        dimk = size(preL{k}, 1);
        if isempty(fix{k})
            parent = zeros(dimk, nu);
        else
            parent = fix{k}(:, undec);
        end
        newfix{k} = [parent, parent];          % [active children | inactive children]
    end
    for i = 1:nu
        k = splitK(i); j = splitJ(i);
        newfix{k}(j, i)      = 1;               % active child  (column i)
        newfix{k}(j, nu + i) = -1;              % inactive child (column nu+i)
    end
    fix = newfix;
end

% ------------------------------------------------------------------------------------
function [cex, lab] = i_find_cex(ops, LB, UB, trueLabel, nSample, precision)
% Exact forward eval at the box center + random in-box samples; return the first input
% that misclassifies (a genuine witness -> sound 'unsafe'). 'like' LB keeps the samples on
% the same device as the input (no host<->GPU transfer for the GPU path).
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
