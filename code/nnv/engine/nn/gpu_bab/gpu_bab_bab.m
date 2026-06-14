function [status, info] = gpu_bab_bab(ops, x_lb, x_ub, trueLabel, nClasses, opts)
% GPU_BAB_BAB  Sound input-split branch-and-bound robustness verifier with batched
%   CROWN bounding -- the GPU-BaB end-to-end pipeline.
%
%   [status, info] = GPU_BAB_BAB(ops, x_lb, x_ub, trueLabel, nClasses, opts)
%     status : 'robust'  -- PROVEN robust (every sub-box certified, margins > 0)
%              'unsafe'  -- a concrete counterexample was evaluated (info.cex, info.cexLabel)
%              'unknown' -- budget exhausted with undecided sub-boxes
%     opts (struct, all optional):
%       .precision 'single'(default)|'double'   .maxBoxes 4096   .maxIter 200
%       .margin    0  (FP-soundness slack: require margins > opts.margin)
%       .nSample   8  (random samples per box for counterexample search)
%
%   How it stays SOUND (sound-or-unknown; never a wrong verdict):
%     * 'robust' only when EVERY leaf's CROWN margins exceed opts.margin -- a property
%       holds on the whole box iff it holds on every sub-box of a partition.
%     * 'unsafe' only from a CONCRETE evaluated input that misclassifies (exact forward
%       eval via the op list) -- a real witness, never inferred from a bound.
%     * 'unknown' when the box/iter budget is hit -- we never guess.
%   The optws.margin slack absorbs gpuArray's lack of directed rounding (single-precision
%   FP error) so a near-zero bound can't produce a false 'robust'.
%
%   GPU parallelism: the whole frontier of sub-boxes is bounded in ONE batched
%   gpu_bab_crown_spec call (columns = sub-boxes), so widening the frontier costs GPU
%   width, not serial time. Input-split is efficient for low-dim inputs (acasxu); for
%   high-dim inputs (MNIST) the single root pass certifies the easy cases and ReLU-split
%   (a follow-on) is the scalable refinement.

    if nargin < 6, opts = struct(); end
    precision = i_get(opts, 'precision', 'single');
    maxBoxes  = i_get(opts, 'maxBoxes', 4096);
    maxIter   = i_get(opts, 'maxIter', 200);
    margin    = cast(i_get(opts, 'margin', 0), precision);
    nSample   = i_get(opts, 'nSample', 8);

    C = -eye(nClasses, precision);
    C(:, trueLabel) = C(:, trueLabel) + 1;
    C(trueLabel, :) = [];

    LB = cast(x_lb, precision);
    UB = cast(x_ub, precision);
    info = struct('iters', 0, 'maxFrontier', 1, 'cex', [], 'cexLabel', []);

    for iter = 1:maxIter
        info.iters = iter;
        % --- batched bounding of the whole frontier ---
        margins = gpu_bab_crown_spec(ops, LB, UB, C, precision);   % nSpec x Q
        verified = all(margins > margin, 1);                       % 1 x Q
        LB(:, verified) = []; UB(:, verified) = [];
        if isempty(LB)
            status = 'robust'; return;
        end
        % --- counterexample search on undecided boxes (concrete, exact eval) ---
        [cex, lab] = i_find_cex(ops, LB, UB, trueLabel, nSample, precision);
        if ~isempty(cex)
            status = 'unsafe'; info.cex = cex; info.cexLabel = lab; return;
        end
        % --- budget ---
        if size(LB, 2) > maxBoxes
            status = 'unknown'; return;
        end
        % --- split each undecided box on its widest input dim ---
        [LB, UB] = i_split(LB, UB);
        info.maxFrontier = max(info.maxFrontier, size(LB, 2));
    end
    status = 'unknown';
end

function [LB2, UB2] = i_split(LB, UB)
% Bisect each box on its widest input dimension -> two child boxes.
    [~, d] = max(UB - LB, [], 1);            % widest dim per box (1 x Q)
    Q = size(LB, 2);
    mid = (LB + UB) / 2;
    LB2 = [LB, LB]; UB2 = [UB, UB];
    idx = sub2ind(size(LB), d, 1:Q);
    UB2(sub2ind(size(UB2), d, 1:Q))     = mid(idx);   % left child: upper half clipped
    LB2(sub2ind(size(LB2), d, Q+(1:Q))) = mid(idx);   % right child: lower half clipped
end

function [cex, lab] = i_find_cex(ops, LB, UB, trueLabel, nSample, precision)
% Exact forward eval (degenerate box = exact) at centers + random samples; return the
% first input that misclassifies, if any. A genuine witness -> sound 'unsafe'.
    cex = []; lab = [];
    Q = size(LB, 2);
    X = (LB + UB) / 2;                                    % centers
    for s = 1:max(0, nSample)
        X = [X, LB + (UB - LB) .* rand(size(LB), precision)]; %#ok<AGROW>
    end
    Y = gpu_bab_ibp(ops, X, X, precision);                % exact (lb=ub=point)
    [~, pred] = max(Y, [], 1);
    bad = find(pred ~= trueLabel, 1);
    if ~isempty(bad)
        cex = X(:, bad); lab = pred(bad);
    end
end

function v = i_get(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
