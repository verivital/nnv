function [status, info] = gpu_bab_relu_split(ops, x_lb, x_ub, trueLabel, nClasses, opts)
% GPU_BAB_RELU_SPLIT  Sound ReLU-split branch-and-bound robustness verifier
%   (beta-CROWN-style) -- branch on unstable NEURONS, not input dimensions, so it
%   scales to high-dim inputs (MNIST) where gpu_bab_bab's input-split cannot.
%
%   [status, info] = GPU_BAB_RELU_SPLIT(ops, x_lb, x_ub, trueLabel, nClasses, opts)
%     status : 'robust' | 'unsafe' (info.cex misclassifies) | 'unknown' (budget)
%     opts   : .precision 'single'(def)|'double'  .maxNodes 500  .margin 0  .nSample 16
%
%   Method: each BaB node is the original input box plus a set of neuron fixings.
%   A node is bounded by CROWN with the fixed neurons' pre-activation bounds CLAMPED
%   (active fix z>=0 -> l := max(l,0); inactive fix z<=0 -> u := min(u,0)). Clamping
%   both pins the neuron's relaxation to stable (exact on its sub-domain) AND propagates
%   the split constraint forward to tighten downstream bounds. An undecided node is split
%   on its first unfixed unstable neuron into an active child and an inactive child.
%
%   SOUNDNESS (sound-or-unknown; never a wrong verdict):
%     * the active/inactive children partition the node (every unstable neuron is z>=0
%       or z<0), so 'robust' requires EVERY leaf certified -- the property holds on the
%       box iff it holds on every leaf;
%     * clamping l up / u down is a valid (tighter) bound on the child's sub-domain, so
%       each leaf's CROWN margin is a sound lower bound there;
%     * 'unsafe' only from a concretely-evaluated misclassifying input on the ORIGINAL
%       box (a real witness for the whole query);
%     * 'unknown' on budget, or if a node has no unstable neurons left yet CROWN still
%       cannot certify it (the convex-relaxation barrier -- needs a tighter bound, not
%       more splitting).

    if nargin < 6, opts = struct(); end
    precision = i_get(opts, 'precision', 'single');
    maxNodes  = i_get(opts, 'maxNodes', 500);
    margin    = cast(i_get(opts, 'margin', 0), precision);
    nSample   = i_get(opts, 'nSample', 16);
    alphaIter = i_get(opts, 'alphaIter', 0);   % >0 -> alpha-CROWN node bounds (tighter -> fewer splits)
    alphaLr   = i_get(opts, 'alphaLr', 0.2);
    cexEvery  = i_get(opts, 'cexEvery', 25);   % per-node counterexample cadence (0 = off)
    intermediate = i_get(opts, 'intermediate', 'ibp');  % 'ibp' | 'tight' (CROWN intermediate bounds)

    C = -eye(nClasses, precision);
    C(:, trueLabel) = C(:, trueLabel) + 1;
    C(trueLabel, :) = [];

    reluIdx = find(cellfun(@(o) strcmp(o.type,'relu'), ops));
    info = struct('nodes', 0, 'maxDepth', 0, 'cex', [], 'cexLabel', []);

    % concrete counterexample on the whole box first (cheap falsification)
    [cex, lab] = i_find_cex(ops, x_lb, x_ub, trueLabel, nSample, precision);
    if ~isempty(cex)
        status = 'unsafe'; info.cex = cex; info.cexLabel = lab; return;
    end

    % DFS over fixings. A fixing is a cell (one per relu op index) of -1/0/+1 vectors.
    fix0 = cell(numel(ops), 1);
    stack = {struct('fix', {fix0}, 'depth', 0)};
    while ~isempty(stack)
        node = stack{end}; stack(end) = [];
        info.nodes = info.nodes + 1;
        info.maxDepth = max(info.maxDepth, node.depth);
        if info.nodes > maxNodes
            status = 'unknown'; return;
        end
        if strcmp(intermediate, 'tight')
            [margins, ~, ~, unstable] = gpu_bab_crown_tight(ops, x_lb, x_ub, C, precision, node.fix);
        elseif alphaIter > 0
            [margins, unstable] = gpu_bab_crown_alpha_fix(ops, x_lb, x_ub, C, node.fix, reluIdx, precision, alphaIter, alphaLr);
        else
            [margins, unstable] = i_crown_clamped(ops, x_lb, x_ub, C, precision, node.fix, reluIdx);
        end
        if all(margins > margin)
            continue;                                   % leaf certified
        end
        % periodic concrete counterexample search (falsify non-robust queries)
        if cexEvery > 0 && mod(info.nodes, cexEvery) == 0
            [cex, lab] = i_find_cex(ops, x_lb, x_ub, trueLabel, nSample, precision);
            if ~isempty(cex), status = 'unsafe'; info.cex = cex; info.cexLabel = lab; return; end
        end
        [kop, j] = i_pick_split(ops, unstable, node.fix, reluIdx);
        if isempty(kop)
            status = 'unknown'; return;                 % no unstable left, still not certified
        end
        fa = node.fix; if isempty(fa{kop}), fa{kop} = zeros(i_dim(unstable,kop),1); end; fa{kop}(j) = 1;
        fi = node.fix; if isempty(fi{kop}), fi{kop} = zeros(i_dim(unstable,kop),1); end; fi{kop}(j) = -1;
        stack{end+1} = struct('fix', {fa}, 'depth', node.depth+1); %#ok<AGROW>
        stack{end+1} = struct('fix', {fi}, 'depth', node.depth+1); %#ok<AGROW>
    end
    status = 'robust';
end

function [margins, unstable] = i_crown_clamped(ops, x_lb, x_ub, C, precision, fixings, reluIdx)
% Forward IBP with the fixed neurons' pre-activation bounds clamped, then backward CROWN.
    nOps = numel(ops);
    preL = cell(nOps,1); preU = cell(nOps,1); unstable = cell(nOps,1);
    lb = cast(x_lb, precision); ub = cast(x_ub, precision);
    for k = 1:nOps
        op = ops{k};
        if strcmp(op.type, 'affine')
            W = cast(op.W, precision); b = cast(op.b(:), precision);
            Wp = max(W,0); Wn = min(W,0);
            nlb = Wp*lb + Wn*ub + b; nub = Wp*ub + Wn*lb + b; lb = nlb; ub = nub;
        else
            r = find(reluIdx == k, 1);
            fx = fixings{k};
            if ~isempty(fx)
                lb(fx == 1)  = max(lb(fx == 1),  0);    % active: z >= 0
                ub(fx == -1) = min(ub(fx == -1), 0);    % inactive: z <= 0
            end
            preL{k} = lb; preU{k} = ub;
            unstable{k} = (lb < 0) & (ub > 0);
            lb = max(lb, 0); ub = max(ub, 0);
        end
    end
    % backward CROWN (single box, lower bound on C*f)
    nSpec = size(C,1);
    A = cast(C, precision); d = zeros(nSpec, 1, precision);
    for k = nOps:-1:1
        op = ops{k};
        if strcmp(op.type, 'affine')
            W = cast(op.W, precision); b = cast(op.b(:), precision);
            d = d + A*b; A = A*W;
        else
            l = preL{k}; u = preU{k}; m = numel(l);
            au = zeros(m,1,precision); bu = zeros(m,1,precision); al = zeros(m,1,precision);
            act = (l >= 0); au(act) = 1; al(act) = 1;
            uns = (l < 0) & (u > 0); dn = u(uns) - l(uns);
            au(uns) = u(uns)./dn; bu(uns) = -au(uns).*l(uns);
            al(uns) = cast(u(uns) >= -l(uns), precision);
            Apos = max(A,0); Aneg = min(A,0);
            d = d + Aneg*bu;
            A = Apos .* al.' + Aneg .* au.';
        end
    end
    Apos = max(A,0); Aneg = min(A,0);
    margins = Apos*cast(x_lb,precision) + Aneg*cast(x_ub,precision) + d;
end

function [kop, j] = i_pick_split(ops, unstable, fixings, reluIdx)
% First unfixed unstable neuron, earliest layer (simple, sound; smarter heuristics later).
    kop = []; j = [];
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        uns = unstable{k};
        fx = fixings{k};
        if isempty(fx), fx = zeros(numel(uns),1); end
        cand = find(uns(:) & (fx(:) == 0), 1);
        if ~isempty(cand)
            kop = k; j = cand; return;
        end
    end
end

function dimv = i_dim(unstable, kop)
    dimv = numel(unstable{kop});
end

function [cex, lab] = i_find_cex(ops, LB, UB, trueLabel, nSample, precision)
    cex = []; lab = [];
    X = (LB + UB) / 2;
    for s = 1:max(0, nSample)
        X = [X, LB + (UB - LB) .* rand(size(LB), precision)]; %#ok<AGROW>
    end
    Y = gpu_bab_ibp(ops, X, X, precision);
    [~, pred] = max(Y, [], 1);
    bad = find(pred ~= trueLabel, 1);
    if ~isempty(bad), cex = X(:, bad); lab = pred(bad); end
end

function v = i_get(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
