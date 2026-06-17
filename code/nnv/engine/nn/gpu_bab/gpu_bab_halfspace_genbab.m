function [verdict, info] = gpu_bab_halfspace_genbab(ops, lb, ub, Gd, gd, opts)
% GPU_BAB_HALFSPACE_GENBAB  GenBaB branch-and-bound for BILINEAR general-halfspace specs
%   (lsnc_relu / linearizenn). Branches over a FIXED net-input box on two axes, the two sources of
%   relaxation looseness in a ReLU+bilinear net:
%     * PRODUCT-input value-ranges (`mulFix`): split the WORST product element (largest McCormick gap
%       (xu-xl)(yu-yl)) at its wider input's midpoint -> that element's gap HALVES per split
%       (geometric bilinear tightening), and the two children partition inputs by v(x)<=mid/>=mid.
%     * ReLU phases (`fixings`): split the unstable neuron with the largest CROWN relaxation gap
%       -pl*pu/(pu-pl) into active (+1) / inactive (-1) children (standard ReLU-split BaB).
%   Each node calls gpu_bab_crown_tight(ops, lb, ub, C, prec, fixings, mulFix) and certifies every
%   unsafe disjunct by single-row separation OR Clip-and-Verify. Products are tightened first
%   (until every gap < gapTol), then ReLUs. The net-input box NEVER changes, so a `mulFix` override
%   set at an ancestor stays valid at every descendant -> no stale-override churn (unlike net-input
%   bisection, which manufactured ~50% empty nodes).
%
%   [verdict, info] = ... ; verdict = 'robust' | 'unknown'.
%
%   SOUNDNESS (sound-or-unknown; never a wrong unsat). Two facts compose:
%     (1) REGION-VALIDITY. At a node, crown_tight's margin m and input-space lower plane (Ain,din)
%         under the node's overrides+fixings are a valid lower bound on C*y over the node's REGION
%         R = {x in box : v_k(x) in [lo_k,hi_k] for each product override k, and each fixed ReLU in
%         its phase} -- NOT over the whole box. A product override clamps the McCormick to [lo,hi],
%         a valid x*y envelope exactly where v(x) in [lo,hi]; a ReLU fixing clamps the relaxation to
%         its phase. So certify (single-row m(row)>g, OR clip-empty) certifies that R avoids the
%         unsafe disjunct -- a clip over the FULL box returning empty is a STRONGER claim that still
%         implies emptiness over R.
%     (2) COVERAGE. The leaves PARTITION the box: a product split gives v<=mid UNION v>=mid, a ReLU
%         split active UNION inactive, and the prune removes only empty regions.
%   Together: any input x* lies in exactly one leaf's region; if x* were a real violation it lands in
%   THAT leaf's region where the bound is valid, so margin <= C*y(x*) <= g and x* sits in the clip
%   set -> that leaf CANNOT certify (it splits or hits maxNodes -> 'unknown'). Hence the off-region
%   over-promise of a McCormick plane is harmless: a SIBLING leaf certifying its own region never
%   absolves the leaf that actually owns x*. 'robust' is returned ONLY when every leaf is certified
%   or pruned, so the whole box is covered by safe regions. (Empirically stress-tested by
%   test_gpu_bab_halfspace_genbab/test_fuzz_never_false_robust.)
%   The infeasibility prune is vacuously safe (CROWN proves the region empty) and fires only past a
%   PRECISION/MAGNITUDE-aware margin (unit roundoff * magnitude, floored at 1e-9), so a live node is
%   never dropped in single or double precision. FP64 by default; never a false 'empty'/separation.
%
%   opts: .maxNodes (3e5) .timeCap (Inf) .tol (1e-6) .precision ('double') .gapTol (1e-7)
%         .minWidth (1e-9, smallest product value-range still worth splitting).
%
%   STATUS / LIMITATION (empirical, lsnc_relu quadrotor2d, 2026-06-17): this is the SOUND GenBaB
%   FOUNDATION, not yet a converging lsnc verifier. With FIXED-slope crown_tight bounds the root
%   spec margin is ~ -41 (very loose), and product + ReLU-phase branching does NOT close it: full-node
%   certifications stay ~0 from 10k to 40k nodes while the tree grows (>1 child/node) -> diverges to
%   'unknown'. The looseness is the single-neuron CROWN relaxation (Salman/Mao barrier), the same wall
%   cersyve hits. Only INPUT-box bisection tightens it (box-shrinking tightens CROWN) but is
%   impractically slow (see gpu_bab_halfspace_input_bab). A converging lsnc verifier needs the full
%   abCROWN GenBaB: alpha/beta-CROWN optimized bounds (extend gpu_bab_crown_alpha_dag to 'product')
%   + alpha-optimized McCormick corners + this targeted product-input branching together. The
%   product-input override (crown_tight mulFix) and the infeasibility prune below are the validated,
%   sound building blocks for that; they are correct in isolation (test_gpu_bab_halfspace_genbab).

    if nargin < 6, opts = struct(); end
    maxNodes = i_g(opts, 'maxNodes', 300000);
    timeCap  = i_g(opts, 'timeCap',  Inf);
    tol      = i_g(opts, 'tol',      1e-6);
    prec     = i_g(opts, 'precision','double');
    gapTol   = i_g(opts, 'gapTol',   1e-7);
    minWidth = i_g(opts, 'minWidth', 1e-9);

    nOps = numel(ops);
    C    = cat(1, Gd{:}); gAll = cat(1, gd{:});
    rows = cell(1, numel(Gd)); r0 = 0;
    for d = 1:numel(Gd), nr = size(Gd{d}, 1); rows{d} = r0 + (1:nr); r0 = r0 + nr; end
    lb = double(lb(:)); ub = double(ub(:));
    prodIdx = find(cellfun(@(o) strcmp(o.type, 'product'), ops));   % bilinear ops
    reluIdx = find(cellfun(@(o) strcmp(o.type, 'relu'),    ops));   % ReLU ops

    info = struct('nodes', 0, 'maxStack', 1, 'sep', 0, 'clip', 0, 'pruned', 0, ...
                  'prodSplit', 0, 'reluSplit', 0, 'reason', '');

    stack = {struct('fixings', {cell(nOps, 1)}, 'mf', {cell(nOps, 1)})};   % root: free, no overrides
    t0 = tic;
    while ~isempty(stack)
        nd = stack{end}; stack(end) = [];
        info.nodes = info.nodes + 1;
        info.maxStack = max(info.maxStack, numel(stack) + 1);
        if info.nodes > maxNodes, verdict = 'unknown'; info.reason = 'maxNodes'; return; end
        if toc(t0) > timeCap,     verdict = 'unknown'; info.reason = 'timeCap'; return; end

        [m, preL, preU, unstable, Ain, din] = gpu_bab_crown_tight(ops, lb, ub, C, prec, nd.fixings, nd.mf);
        m = double(m(:)); Ain = double(Ain); din = double(din(:));

        % Infeasible node: a stale product override (after a later branch changed the bounds) OR a
        % ReLU fixing the box can't satisfy (e.g. forcing active when pu<0) makes CROWN give
        % preL{k}>preU{k}. Scan BOTH the product ops (override case) and the ReLU ops (fixing case);
        % either proves the node empty -> prune, vacuously safe. The empty-margin is PRECISION- and
        % MAGNITUDE-aware (unit roundoff u * 64 safety * bound magnitude, floored at 1e-9): the prune
        % fires only when the lower bound exceeds the upper by MORE than the worst-case accumulated
        % CROWN rounding error, so a genuinely non-empty node (where preL<=preU in exact arithmetic)
        % is NEVER dropped, in single OR double precision and at any magnitude. A fixed 1e-9 was
        % unsound for precision='single' / huge-magnitude nets.
        u = 1.1e-16; if strcmp(prec, 'single'), u = 6e-8; end
        empty = false;
        for k = [prodIdx(:); reluIdx(:)].'
            pl = preL{k}; pu = preU{k};
            if isempty(pl), continue; end
            emargin = max(1e-9, 64 * u * max([1; abs(pl(:)); abs(pu(:))]));
            if any(pl > pu + emargin), empty = true; break; end
        end
        if empty, info.pruned = info.pruned + 1; continue; end

        ok = true;
        for d = 1:numel(Gd)
            rr = rows{d};
            if any(m(rr) - gAll(rr) > tol), info.sep = info.sep + 1; continue; end
            [~, ~, isEmpty] = gpu_bab_clip(lb, ub, Ain(rr, :), din(rr) - gAll(rr));
            if isEmpty, info.clip = info.clip + 1; else, ok = false; break; end
        end
        if ok, continue; end                                       % node certified

        % --- branch: products first (until tight), then ReLU phases ---
        bestGap = -inf; bk = 0; bIdx = 0; vlo = 0; vhi = 0;
        for kk = 1:numel(prodIdx)
            k = prodIdx(kk); wa = ops{k}.sizes(1);
            la = preL{k}(1:wa); ua = preU{k}(1:wa); ly = preL{k}(wa+1:end); uy = preU{k}(wa+1:end);
            wx = ua - la; wy = uy - ly; gap = wx .* wy;
            [gmax, j] = max(gap);
            if gmax > bestGap
                bestGap = gmax; bk = k;
                if wx(j) >= wy(j), bIdx = j;      vlo = la(j); vhi = ua(j);
                else,              bIdx = wa + j; vlo = ly(j); vhi = uy(j); end
            end
        end

        if bestGap > gapTol && (vhi - vlo) > minWidth
            % PRODUCT-INPUT split: halve the worst product input's value-range
            mid = (vlo + vhi) / 2;
            wstack = sum(ops{bk}.sizes);            % stacked [in1;in2] width (inputs are equal-width)
            base = nd.mf{bk};
            if isempty(base), base = struct('lo', -inf(wstack, 1), 'hi', inf(wstack, 1)); end
            cLo = nd; cLo.mf{bk} = base; cLo.mf{bk}.hi(bIdx) = min(base.hi(bIdx), mid);
            cHi = nd; cHi.mf{bk} = base; cHi.mf{bk}.lo(bIdx) = max(base.lo(bIdx), mid);
            stack{end+1} = cLo; stack{end+1} = cHi; %#ok<AGROW>
            info.prodSplit = info.prodSplit + 1;
        else
            % RELU-PHASE split: the unstable neuron with the largest CROWN relaxation gap
            bestRG = -inf; rk = 0; rIdx = 0;
            for kk = 1:numel(reluIdx)
                k = reluIdx(kk);
                if isempty(unstable{k}) || ~any(unstable{k}), continue; end
                pl = preL{k}; pu = preU{k}; un = unstable{k};
                g = -inf(numel(pl), 1);
                g(un) = -pl(un) .* pu(un) ./ (pu(un) - pl(un));    % CROWN triangle max gap
                [gm, j] = max(g);
                if gm > bestRG, bestRG = gm; rk = k; rIdx = j; end
            end
            if rk == 0
                verdict = 'unknown'; info.reason = 'barrier (products+ReLUs tight, undecided)'; return;
            end
            base = nd.fixings{rk};
            if isempty(base), base = zeros(numel(preL{rk}), 1); end
            cA = nd; cA.fixings{rk} = base; cA.fixings{rk}(rIdx) =  1;   % active
            cI = nd; cI.fixings{rk} = base; cI.fixings{rk}(rIdx) = -1;   % inactive
            stack{end+1} = cA; stack{end+1} = cI; %#ok<AGROW>
            info.reluSplit = info.reluSplit + 1;
        end
    end
    verdict = 'robust'; info.reason = 'genbab covered (all leaves certified or pruned)';
end

function v = i_g(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
