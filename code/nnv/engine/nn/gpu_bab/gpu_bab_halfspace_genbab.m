function [verdict, info] = gpu_bab_halfspace_genbab(ops, lb, ub, Gd, gd, opts)
% GPU_BAB_HALFSPACE_GENBAB  GenBaB-style branch-and-bound for BILINEAR general-halfspace specs
%   (lsnc_relu / linearizenn): targeted **product-input branching** instead of naive net-input
%   bisection. Each node carries the (fixed) net-input box plus `mulFix` = per-'product' value-range
%   overrides. Per node: gpu_bab_crown_tight(...,mulFix) gives the bound + each product input's
%   range; each unsafe disjunct is certified by single-row separation OR Clip-and-Verify. If
%   undecided, the WORST product element (largest McCormick gap (xu-xl)(yu-yl)) has its WIDER input's
%   value-range split at the midpoint into two children -> that element's McCormick gap HALVES, so
%   the bilinear looseness shrinks GEOMETRICALLY per split (vs the quadratic-in-dimension cost of
%   bisecting the 6-D net input). When all products are tight but a node is still undecided (residual
%   ReLU looseness), it falls back to net-input bisection.
%
%   [verdict, info] = ... ; verdict = 'robust' | 'unknown'.
%
%   SOUNDNESS (sound-or-unknown; never a wrong unsat): a product-input split partitions inputs by
%   v(x) <= mid / v(x) >= mid (every x falls in exactly one half), and crown_tight clamps that
%   product input's range to [lo,mid]/[mid,hi] per child -> the spec lower bound with the clamp is a
%   valid lower bound for the inputs whose v lies in that half (the McCormick over the narrower range
%   is a valid over-approx there). A net-input bisection split is the usual closed box split. 'robust'
%   is returned ONLY when every leaf node is certified (single-row OR clip), so the whole input box is
%   covered. FP64; never a false 'empty'/separation -> never a wrong unsat.
%
%   opts: .maxNodes (2e5) .timeCap (Inf) .tol (1e-6) .precision ('double') .minWidth (1e-9)
%         .gapTol (1e-9, switch to net-input bisection when the worst product gap is below this).

    if nargin < 6, opts = struct(); end
    maxNodes = i_g(opts, 'maxNodes', 200000);
    timeCap  = i_g(opts, 'timeCap',  Inf);
    tol      = i_g(opts, 'tol',      1e-6);
    prec     = i_g(opts, 'precision','double');
    minWidth = i_g(opts, 'minWidth', 1e-9);
    gapTol   = i_g(opts, 'gapTol',   1e-9);

    nOps = numel(ops);
    C    = cat(1, Gd{:}); gAll = cat(1, gd{:});
    rows = cell(1, numel(Gd)); r0 = 0;
    for d = 1:numel(Gd), nr = size(Gd{d}, 1); rows{d} = r0 + (1:nr); r0 = r0 + nr; end
    lb = double(lb(:)); ub = double(ub(:));
    prodIdx = find(cellfun(@(o) strcmp(o.type, 'product'), ops));   % bilinear ops to branch

    info = struct('nodes', 0, 'maxStack', 1, 'sep', 0, 'clip', 0, 'pruned', 0, 'prodSplit', 0, 'inpSplit', 0, 'reason', '');

    stack = {struct('lb', lb, 'ub', ub, 'mf', {cell(nOps, 1)})};   % root: full box, no product overrides
    t0 = tic;
    while ~isempty(stack)
        nd = stack{end}; stack(end) = [];
        info.nodes = info.nodes + 1;
        info.maxStack = max(info.maxStack, numel(stack) + 1);
        if info.nodes > maxNodes, verdict = 'unknown'; info.reason = 'maxNodes'; return; end
        if toc(t0) > timeCap,     verdict = 'unknown'; info.reason = 'timeCap'; return; end

        [m, preL, preU, ~, Ain, din] = gpu_bab_crown_tight(ops, nd.lb, nd.ub, C, prec, {}, nd.mf);
        m = double(m(:)); Ain = double(Ain); din = double(din(:));

        % Infeasible sub-problem: a stale mulFix override (set at an ancestor, before a later
        % net-input bisection) clamped a product value-range to a half its child box can't reach.
        % CROWN proves it empty (computed preL > clamped preU ==> the true value exceeds the clamp
        % everywhere in this box ==> no input satisfies the split) -> prune, vacuously safe. SOUND:
        % only fires when a real coordinate's lower bound exceeds its upper by a true margin, never
        % an eps-level FP wobble (tol 1e-9), so a non-empty node is never dropped.
        empty = false;
        for kk = 1:numel(prodIdx)
            if any(preL{prodIdx(kk)} > preU{prodIdx(kk)} + 1e-9), empty = true; break; end
        end
        if empty, info.pruned = info.pruned + 1; continue; end

        ok = true;
        for d = 1:numel(Gd)
            rr = rows{d};
            if any(m(rr) - gAll(rr) > tol), info.sep = info.sep + 1; continue; end
            [~, ~, isEmpty] = gpu_bab_clip(nd.lb, nd.ub, Ain(rr, :), din(rr) - gAll(rr));
            if isEmpty, info.clip = info.clip + 1; else, ok = false; break; end
        end
        if ok, continue; end                                       % node certified

        % --- pick the worst product element (largest McCormick gap), branch its wider input ---
        bestGap = -inf; bk = 0; bIdx = 0; vlo = 0; vhi = 0;
        for kk = 1:numel(prodIdx)
            k = prodIdx(kk); wa = ops{k}.sizes(1);
            la = preL{k}(1:wa); ua = preU{k}(1:wa); ly = preL{k}(wa+1:end); uy = preU{k}(wa+1:end);
            wx = ua - la; wy = uy - ly; gap = wx .* wy;
            [gmax, j] = max(gap);
            if gmax > bestGap
                bestGap = gmax; bk = k;
                if wx(j) >= wy(j), bIdx = j;      vlo = la(j); vhi = ua(j);     % branch input 1, elt j
                else,              bIdx = wa + j; vlo = ly(j); vhi = uy(j); end % branch input 2
            end
        end

        if bestGap > gapTol && (vhi - vlo) > minWidth
            % PRODUCT-INPUT split: halve the worst product input's value-range
            mid = (vlo + vhi) / 2;
            wstack = 2 * ops{bk}.sizes(1);
            base = nd.mf{bk};
            if isempty(base), base = struct('lo', -inf(wstack, 1), 'hi', inf(wstack, 1)); end
            cLo = nd; cLo.mf{bk} = base; cLo.mf{bk}.hi(bIdx) = min(base.hi(bIdx), mid);   % v in [.,mid]
            cHi = nd; cHi.mf{bk} = base; cHi.mf{bk}.lo(bIdx) = max(base.lo(bIdx), mid);   % v in [mid,.]
            stack{end+1} = cLo; stack{end+1} = cHi; %#ok<AGROW>
            info.prodSplit = info.prodSplit + 1;
        else
            % FALLBACK net-input bisection (products tight, residual ReLU looseness)
            w = nd.ub - nd.lb; [wmax, jdim] = max(w);
            if wmax < minWidth
                if ~isempty(getenv('NNV_GENBAB_DBG'))
                    bw = inf; for dd = 1:numel(Gd), bw = min(bw, max(m(rows{dd}) - gAll(rows{dd}))); end
                    fprintf('[genbab] BARRIER worstDisjunctBestRowMargin=%.6g boxWidth=%.2e prodGap=%.2e nodes=%d\n', bw, wmax, bestGap, info.nodes);
                end
                verdict = 'unknown'; info.reason = 'barrier (box+products tight, undecided)'; return;
            end
            md = (nd.lb(jdim) + nd.ub(jdim)) / 2;
            c1 = nd; c1.ub(jdim) = md; c2 = nd; c2.lb(jdim) = md;
            stack{end+1} = c1; stack{end+1} = c2; %#ok<AGROW>
            info.inpSplit = info.inpSplit + 1;
        end
    end
    verdict = 'robust'; info.reason = 'genbab covered (all leaves certified)';
end

function v = i_g(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
