function [verdict, info] = gpu_bab_halfspace_input_bab(ops, lb, ub, Gd, gd, opts)
% GPU_BAB_HALFSPACE_INPUT_BAB  Input-domain bisection branch-and-bound for general-halfspace
%   specs whose relaxation looseness shrinks with the INPUT box width -- in particular BILINEAR
%   / product nets (e.g. lsnc Lyapunov), where the McCormick gap is O((box width)^2) so halving
%   an input dimension quarters the gap. Per sub-box it computes a tight CROWN bound
%   (gpu_bab_crown_tight, which handles the 'product' McCormick op) and certifies each unsafe
%   disjunct by EITHER single-row separation (some row's sound lower bound G_i*y - g_i > tol) OR
%   the Clip-and-Verify joint emptiness test (gpu_bab_clip on the disjunct's conjunction). An
%   undecided sub-box is split on its widest input dimension.
%
%   [verdict, info] = ... ; verdict = 'robust' | 'unknown'.
%
%   SOUNDNESS (sound-or-unknown): the two children of a bisection PARTITION the parent box
%   (closed split [lb,mid] U [mid,ub]); 'robust' is returned ONLY when the stack empties with
%   EVERY popped sub-box certified avoided, so the whole input box is covered by certified-safe
%   sub-boxes. Both certification tests are sound (single-row = a sound CROWN lower bound;
%   gpu_bab_clip never reports a false "empty"). FP64. Never a wrong unsat.
%
%   WHEN TO USE: only when the gap shrinks with box width (bilinear/product). For pure-ReLU nets
%   (e.g. cersyve) the looseness lives in the ReLU relaxation, not the box, and input bisection
%   does NOT converge -- the caller should gate this on the presence of a 'product' op.
%
%   opts (all optional): .maxNodes (2e5) .timeCap (Inf, s) .tol (1e-6) .precision ('double')
%                        .minWidth (1e-9, barrier guard).

    if nargin < 6, opts = struct(); end
    maxNodes = i_g(opts, 'maxNodes', 200000);
    timeCap  = i_g(opts, 'timeCap',  Inf);
    tol      = i_g(opts, 'tol',      1e-6);
    prec     = i_g(opts, 'precision','double');
    minWidth = i_g(opts, 'minWidth', 1e-9);

    nOps = numel(ops);
    C    = cat(1, Gd{:});
    gAll = cat(1, gd{:});
    rows = cell(1, numel(Gd)); r0 = 0;
    for d = 1:numel(Gd), nr = size(Gd{d}, 1); rows{d} = r0 + (1:nr); r0 = r0 + nr; end
    lb = double(lb(:)); ub = double(ub(:)); n = numel(lb);

    info = struct('nodes', 0, 'maxStack', 1, 'sep', 0, 'clip', 0, 'reason', '');

    % LIFO stack of sub-boxes as columns of SL/SU (n x cap), amortized-doubling capacity.
    cap = 1024;
    SL = zeros(n, cap); SU = zeros(n, cap);
    SL(:, 1) = lb; SU(:, 1) = ub; top = 1;
    t0 = tic;

    while top > 0
        nl = SL(:, top); nu = SU(:, top); top = top - 1;
        info.nodes = info.nodes + 1;
        if info.nodes > maxNodes, verdict = 'unknown'; info.reason = 'maxNodes'; return; end
        if toc(t0) > timeCap,    verdict = 'unknown'; info.reason = 'timeCap'; return; end

        [m, ~, ~, ~, Ain, din] = gpu_bab_crown_tight(ops, nl, nu, C, prec, cell(nOps, 1));
        m = double(m(:)); Ain = double(Ain); din = double(din(:));

        ok = true;                                   % certify: every disjunct avoided over this sub-box
        for d = 1:numel(Gd)
            rr = rows{d};
            if any(m(rr) - gAll(rr) > tol), info.sep = info.sep + 1; continue; end   % single-row separation
            [~, ~, isEmpty] = gpu_bab_clip(nl, nu, Ain(rr, :), din(rr) - gAll(rr));   % joint clip emptiness
            if isEmpty, info.clip = info.clip + 1; else, ok = false; break; end
        end
        if ok, continue; end                         % sub-box safe -> discard

        w = nu - nl; [wmax, j] = max(w);
        if wmax < minWidth
            verdict = 'unknown'; info.reason = 'barrier (sub-box below minWidth still undecided)'; return;
        end
        mid = (nl(j) + nu(j)) / 2;
        if top + 2 > cap                             % grow (amortized doubling)
            cap = 2 * cap; SL(:, cap) = 0; SU(:, cap) = 0;
        end
        SL(:, top+1) = nl;       SU(:, top+1) = nu;       SU(j, top+1) = mid;   % child 1: [nl, mid]
        SL(:, top+2) = nl;       SU(:, top+2) = nu;       SL(j, top+2) = mid;   % child 2: [mid, nu]
        top = top + 2;
        info.maxStack = max(info.maxStack, top);
    end
    verdict = 'robust'; info.reason = 'input-bisection covered (all sub-boxes avoided)';
end

function v = i_g(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
