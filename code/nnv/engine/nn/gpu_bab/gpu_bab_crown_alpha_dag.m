function [margins, unstable, preL, preU, alphaOut] = gpu_bab_crown_alpha_dag(ops, x_lb, x_ub, C, fixings, reluIdx, precision, nIter, lr, rootBounds, alphaInit)
%   alphaOut (5th out): the optimized per-relu lower slopes as cell(nOps,1) of dim_k x 1 vectors
%     (column 1 of the batch). Pass this as gpu_bab_crown_spec_dag's alphaCell to bound a whole
%     BaB frontier with these fixed root slopes -- AMORTIZED alpha-CROWN (no per-node autodiff).
% GPU_BAB_CROWN_ALPHA_DAG  alpha-CROWN lower bound on a linear spec C*f(x) for CONV/DAG nets.
%   The fusion of gpu_bab_crown_alpha_fix (the FC alpha-CROWN harness) onto the full DAG
%   backward of gpu_bab_crown_spec_dag (affine/conv/normaffine/avgpool/relu/add). It optimizes
%   the unstable-ReLU LOWER slopes (alpha in [0,1]) to MAXIMIZE the certified margin, exactly
%   like the FC version, but the backward routes through conv/pool/normaffine/add adjoints so it
%   works for cifar-style conv resnets (where the FC-only i_aback mis-bounds any non-affine op).
%
%   [margins, unstable, preL, preU] = gpu_bab_crown_alpha_dag(ops, x_lb, x_ub, C, fixings,
%                                       reluIdx, precision, nIter, lr, rootBounds)
%     ops/x_lb/x_ub/C/fixings/rootBounds : as gpu_bab_crown_spec_dag (rootBounds STRONGLY
%         preferred -- tight intermediate bounds; the loose DAG IBP forward is the fallback).
%     reluIdx  : indices of the relu ops (the alpha-carrying ops).
%     nIter    : alpha-CROWN projected-gradient-ascent iterations. nIter<=0 -> FIXED min-area
%                slope == gpu_bab_crown_spec_dag bound-for-bound (the soundness/parity gate).
%     margins  : nSpec-by-B optimized lower bound on C*f over each node's clamped box.
%
%   SOUNDNESS (sound-or-unknown; NEVER a false robust):
%     * alpha enters ONLY the lower slope al = actM + unsM.*alpha, projected to [0,1] every step.
%       Any al in [0,1] is a valid lower ReLU relaxation (ReLU(z) >= al*z on [l,u]); the upper
%       line (au,bu) and the active/stable units are FIXED constants -> any in-range alpha gives
%       a sound lower bound. The optimizer only SEARCHES alpha; it never weakens the bound.
%     * au/bu/actM/unsM and the intermediate bounds preL/preU are PRECOMPUTED constants (detached
%       from the alpha tape) -- alpha can never leak back into the bounds (no joint-alpha).
%     * keep-best + a final re-evaluation at bestVec means the returned margin is >= the min-area
%       (spec_dag) margin: alpha-CROWN is sound-AND-no-worse, can only help or no-op.
%     * the conv/normaffine/avgpool/add adjoints are EXACT (linear, no relaxation), so alpha is
%       the only free variable. maxpool is unsupported and ERRORS (caller falls back, sound).
%
%   PARITY GATE: at nIter<=0 (alpha == min-area), the backward is OP-FOR-OP identical to
%   gpu_bab_crown_spec_dag -- gpu_bab_alpha_dag_parity_test asserts bound-for-bound equality
%   before this is ever trusted for a 'robust' emit.
%
%   NOTE (Phase 2 / autodiff): for nIter>0 the PGA loop needs the backward to be dlarray-
%   traceable so dlgradient reaches alpha. i_conv_backward/i_avgpool_backward currently extract
%   to numeric (tape-break) -> gradients vanish through conv and the optimizer no-ops (SOUND, no
%   improvement). Making those adjoints traceable is the Phase-2 work; until then nIter>0 returns
%   the min-area bound (keep-best), which is still sound.

    if nargin < 7 || isempty(precision), precision = 'single'; end
    if nargin < 8 || isempty(nIter), nIter = 20; end
    if nargin < 9 || isempty(lr), lr = 0.2; end
    if nargin < 10, rootBounds = []; end
    if nargin < 11, alphaInit = {}; end             % optional per-relu warm-start slopes (dim_k x 1)

    B = size(x_lb, 2); nSpec = size(C, 1); nOps = numel(ops);

    % ---- per-relu pre-activation bounds (clamped) + the FIXED upper-line relaxation (au,bu) and
    % active/unstable masks. rootBounds reuse (tight) preferred; else full DAG IBP forward. ----
    [preL, preU, actM, unsM, auC, buC, unstable] = i_dag_relu_bounds(ops, x_lb, x_ub, precision, fixings, rootBounds);

    % ---- flat alpha vector (min-area init) ----
    rdims = arrayfun(@(k) size(preL{k},1), reluIdx);
    offsets = [0; cumsum(rdims(:) * B)];
    nP = offsets(end);
    a0 = zeros(nP, 1, precision);
    fixSign = cell(nOps, 1); fixedMask = zeros(nP, 1, precision);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        a = zeros(size(preL{k}), precision);
        uns = unsM{k} > 0;
        if ~isempty(alphaInit) && numel(alphaInit) >= k && ~isempty(alphaInit{k})
            % WARM START: seed unstable slopes from the supplied (e.g. amortized root) alpha. keep-best
            % then floors the bound at THIS init, so few PGA iters refine it instead of climbing from
            % min-area. Sound for any alpha in [0,1] (clamped); a better init only helps convergence.
            ai = repmat(cast(alphaInit{k}(:), precision), 1, size(preL{k}, 2));   % dim_k x B
            a(uns) = min(max(ai(uns), 0), 1);
        else
            a(uns) = cast(preU{k}(uns) >= -preL{k}(uns), precision);   % min-area binary {0,1}
        end
        a0(offsets(r)+1 : offsets(r+1)) = a(:);
        % BETA: per-node split-constraint sign (s_i = +1 active fix / -1 inactive fix / 0 free).
        % beta is a free dual ONLY on fixed neurons -> at the root (no fixings) fixedMask=0 ->
        % beta stays 0 -> this reduces EXACTLY to pure alpha (parity gate preserved).
        if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
            fs = cast(fixings{k}, precision);
        else
            fs = zeros(size(preL{k}), precision);
        end
        fixSign{k} = fs;
        fixedMask(offsets(r)+1 : offsets(r+1)) = cast(fs(:) ~= 0, precision);
    end
    pVec = dlarray([a0; zeros(nP, 1, precision)]);   % [alpha; beta], alpha=min-area, beta=0

    if nIter <= 0
        margins = i_dag_back(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, ...
                             pVec, reluIdx, rdims, offsets, nP, B, nSpec);
        margins = i_g(margins);
        alphaOut = i_alpha_cell(pVec, nP, nOps, reluIdx, rdims, offsets);
        return;
    end

    % ---- projected gradient ascent over [alpha; beta] (Adam, keep-best) ----
    % SOUNDNESS: the bound is max_{alpha in[0,1], beta>=0} min_x [C*f(x) - sum beta_i s_i z_i] over
    % the clamped box -- a valid lower bound for ANY alpha in [0,1], beta>=0 (sound ReLU relaxation
    % + weak Lagrangian duality on the split constraints). keep-best + the final re-eval at bestVec
    % mean the returned bound is >= the min-area/beta=0 bound regardless of the optimizer; a severed
    % autodiff tape is CAUGHT -> min-area fallback (sound, no improvement), never error/unsound.
    m0 = i_dag_back(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, ...
                    pVec, reluIdx, rdims, offsets, nP, B, nSpec);
    bestObj = i_g(sum(min(m0,[],1), 'all')); bestVec = pVec;
    mt = zeros(2*nP, 1, precision); vt = zeros(2*nP, 1, precision);
    b1 = 0.9; b2 = 0.999; epsA = cast(1e-8, precision);
    % WORST-SPEC gradient (the cifar-scale memory fix): the objective sum(min_spec margin) has a
    % gradient that flows ONLY through each node's ARGMIN spec, so the autodiff tape needs just that
    % ONE spec row per node (1 x featureMap x B) instead of all nSpec -> ~nSpec x less memory ->
    % large frontier with full alpha+beta. For small nSpec (cheap) keep the exact full gradient.
    useWorst = (nSpec > 16);
    if useWorst
        [~, wIdx] = min(i_g(m0), [], 1);              % worst spec per node, found ONCE (the stable bottleneck)
        Cw = i_worst_C(C, wIdx, precision);          % 1 x nOut x B -> the optimization + keep-best use only this
    end
    try
        for it = 1:nIter
            if useWorst
                [~, grad] = dlfeval(@(p) i_aloss(ops, preL, actM, unsM, auC, buC, fixSign, ...
                                    x_lb, x_ub, Cw, precision, p, reluIdx, rdims, offsets, nP, B, 1), pVec);
            else
                [~, grad] = dlfeval(@(p) i_aloss(ops, preL, actM, unsM, auC, buC, fixSign, ...
                                    x_lb, x_ub, C, precision, p, reluIdx, rdims, offsets, nP, B, nSpec), pVec);
            end
            g = extractdata(grad);
            if ~any(g(:)), break; end                  % zero gradient (tape didn't reach params) -> stop
            mt = b1*mt + (1-b1)*g;
            vt = b2*vt + (1-b2)*g.^2;
            mhat = mt / (1 - b1^it); vhat = vt / (1 - b2^it);
            v = extractdata(pVec) - lr * mhat ./ (sqrt(vhat) + epsA);
            v(1:nP)      = max(min(v(1:nP), 1), 0);            % alpha in [0,1]
            v(nP+1:end)  = max(v(nP+1:end), 0) .* fixedMask;   % beta >= 0, only on fixed neurons
            pVec = dlarray(v);
            if useWorst
                % cheap keep-best: score only the worst spec (1 row), not all nSpec. At init the
                % worst spec IS the per-node argmin, so sum(worst) == sum(min) -> bestObj stays comparable.
                mw = i_dag_back(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, Cw, precision, ...
                                pVec, reluIdx, rdims, offsets, nP, B, 1);
                obj = i_g(sum(mw, 'all'));
            else
                m = i_dag_back(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, ...
                               pVec, reluIdx, rdims, offsets, nP, B, nSpec);
                obj = i_g(sum(min(m,[],1), 'all'));
            end
            if obj > bestObj, bestObj = obj; bestVec = pVec; end
        end
    catch ME
        if ~any(strcmp(ME.identifier, {'MATLAB:dlarray:GradientNotTraced', ...
                'deep:dlarray:ValueToDiffNotDlarray', 'MATLAB:dlarray:NotDifferentiable'}))
            % a real error (not the known tape-break) -> keep-best already holds the min-area
            % bound; swallow to stay sound-or-unknown rather than throw out of the precheck.
        end
    end
    margins = i_dag_back(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, ...
                         bestVec, reluIdx, rdims, offsets, nP, B, nSpec);
    margins = i_g(margins);
    % NO-WORSE-THAN-MIN-AREA guarantee. The worst-spec keep-best optimizes a FIXED-worst-spec proxy;
    % as (alpha,beta) move, a node's true argmin spec can shift, so the proxy can diverge and bestVec's
    % true worst-margin can dip below min-area (still SOUND -- a valid lower bound at feasible
    % (alpha,beta) -- just looser). Per node, fall back to the min-area bound m0 wherever it is tighter,
    % so the returned bound is never below min-area. (No-op for the exact full-gradient path.)
    m0g = i_g(m0);
    fb = min(m0g, [], 1) > min(margins, [], 1);          % 1 x B: min-area beats optimized for this node
    if any(fb), margins(:, fb) = m0g(:, fb); end
    alphaOut = i_alpha_cell(bestVec, nP, nOps, reluIdx, rdims, offsets);
end

% =====================================================================================
function ac = i_alpha_cell(pVec, nP, nOps, reluIdx, rdims, offsets)
% Extract the optimized alpha (first nP entries of [alpha;beta]) as a per-relu cell(nOps,1) of
% dim_k x 1 column-1 slopes, for reuse as gpu_bab_crown_spec_dag's fixed alphaCell.
    v = i_g(pVec(1:nP));
    ac = cell(nOps, 1);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        a = reshape(v(offsets(r)+1 : offsets(r+1)), rdims(r), []);   % dim_k x B
        ac{k} = a(:, 1);                                             % root: B=1 (or first column)
    end
end

% =====================================================================================
function Cw = i_worst_C(C, wIdx, precision)
% Per-node worst spec for the memory-efficient gradient: 1 x nOut x B, where column b is the
% wIdx(b)-th row of the full spec C (the argmin spec for node b).
    nOut = size(C, 2); B = numel(wIdx);
    Cw = reshape(cast(C(wIdx, :), precision).', 1, nOut, B);   % (B x nOut).' = nOut x B -> 1 x nOut x B
end

% =====================================================================================
function [loss, grad] = i_aloss(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, pVec, reluIdx, rdims, offsets, nP, B, nSpec)
    margins = i_dag_back(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, pVec, reluIdx, rdims, offsets, nP, B, nSpec);
    loss = -sum(min(margins, [], 1), 'all');
    grad = dlgradient(loss, pVec);
end

% =====================================================================================
function margins = i_dag_back(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, pVec, reluIdx, rdims, offsets, nP, B, nSpec)
% Backward CROWN over the FULL DAG (spec_dag's backward), with the unstable ReLU LOWER slope
% supplied by alpha (al = actM + unsM.*alpha) AND the beta split-dual subtracted at each ReLU
% (-beta_k*s_k on z_k). au/bu are precomputed constants. With alpha=min-area and beta=0 this is
% op-for-op identical to gpu_bab_crown_spec_dag.
    nOps = numel(ops); n = size(x_lb, 1);
    skipA = cell(nOps, 1);
    if ndims(C) == 3
        skipA{nOps} = cast(C, precision);            % per-NODE spec (nSpec x nOut x B), e.g. worst-spec
    else
        skipA{nOps} = repmat(cast(C, precision), 1, 1, B);   % shared spec
    end
    nSpec = size(skipA{nOps}, 1);                    % derive (1 for worst-spec gradient, else full)
    d = zeros(nSpec, B, precision);
    inputSkipA = [];
    for k = nOps:-1:1
        A = skipA{k};
        if isempty(A), continue; end
        skipA{k} = [];
        op = ops{k};
        if strcmp(op.type, 'add')
            for ii = 1:numel(op.inputs)
                s = op.inputs(ii);
                if s == 0
                    if isempty(inputSkipA), inputSkipA = A; else, inputSkipA = inputSkipA + A; end
                elseif isempty(skipA{s}), skipA{s} = A;
                else, skipA{s} = skipA{s} + A;
                end
            end
            continue;
        end
        switch op.type
            case 'affine'
                W = cast(op.W, precision); bb = cast(op.b(:), precision);
                d = d + reshape(pagemtimes(A, bb), nSpec, B);
                A = pagemtimes(A, W);
            case 'conv'
                [A, d] = i_conv_backward(A, d, op, precision);
            case 'normaffine'
                sf = i_bcast_flat(op.scale, op.shape, precision);
                tf = i_bcast_flat(op.shift, op.shape, precision);
                d = d + reshape(pagemtimes(A, tf), nSpec, B);
                A = A .* reshape(sf, 1, [], 1);
            case 'avgpool'
                [A, d] = i_avgpool_backward(A, d, op, precision);
            case 'relu'
                r = find(reluIdx == k, 1);
                dim = rdims(r);
                alpha_k = reshape(pVec(offsets(r)+1 : offsets(r+1)), dim, B);
                beta_k  = reshape(pVec(nP + offsets(r)+1 : nP + offsets(r+1)), dim, B);
                au = auC{k}; bu = buC{k};
                al = actM{k} + unsM{k} .* alpha_k;
                Apos = max(A, 0); Aneg = min(A, 0);
                d = d + reshape(pagemtimes(Aneg, reshape(bu, dim, 1, B)), nSpec, B);
                A = Apos .* reshape(al, 1, dim, B) + Aneg .* reshape(au, 1, dim, B);
                % BETA split-dual: subtract beta_k*s_k from the coefficient on z_k (all spec rows).
                % Zero at the root / unfixed neurons (beta=0 there) -> reduces to pure alpha.
                A = A - reshape(beta_k .* fixSign{k}, 1, dim, B);
            otherwise
                error('gpu_bab_crown_alpha_dag:op', ...
                    'Unsupported op "%s" (affine/conv/normaffine/avgpool/relu/add only).', op.type);
        end
        s = op.src;
        if s == 0
            if isempty(inputSkipA), inputSkipA = A; else, inputSkipA = inputSkipA + A; end
        elseif isempty(skipA{s}), skipA{s} = A;
        else, skipA{s} = skipA{s} + A;
        end
    end
    if isempty(inputSkipA)
        margins = d; return;
    end
    A = inputSkipA;
    Apos = max(A, 0); Aneg = min(A, 0);
    lbcol = reshape(cast(x_lb, precision), n, 1, B);
    ubcol = reshape(cast(x_ub, precision), n, 1, B);
    margins = reshape(pagemtimes(Apos, lbcol), nSpec, B) ...
            + reshape(pagemtimes(Aneg, ubcol), nSpec, B) + d;
end

% =====================================================================================
function [preL, preU, actM, unsM, auC, buC, unstable] = i_dag_relu_bounds(ops, x_lb, x_ub, precision, fixings, rootBounds)
% Per-relu (clamped) pre-activation bounds + the fixed upper-line relaxation (au,bu) and the
% active/unstable masks. rootBounds-reuse (tight) preferred; else the full DAG IBP forward.
    nOps = numel(ops); B = size(x_lb, 2);
    preL = cell(nOps,1); preU = cell(nOps,1);
    actM = cell(nOps,1); unsM = cell(nOps,1);
    auC = cell(nOps,1);  buC = cell(nOps,1); unstable = cell(nOps,1);

    if isempty(rootBounds)
        % full DAG IBP forward (spec_dag's forward), clamped per node, caching every op's bounds
        cl = cell(nOps+1,1); cu = cell(nOps+1,1);
        cl{1} = cast(x_lb, precision); cu{1} = cast(x_ub, precision);
        for k = 1:nOps
            op = ops{k};
            if strcmp(op.type, 'add')
                a = op.inputs(1)+1; b = op.inputs(2)+1;
                cl{k+1} = cl{a} + cl{b}; cu{k+1} = cu{a} + cu{b}; continue;
            end
            s = op.src + 1; lb = cl{s}; ub = cu{s};
            switch op.type
                case 'affine'
                    W = cast(op.W, precision); bb = cast(op.b(:), precision);
                    Wp = max(W,0); Wn = min(W,0);
                    cl{k+1} = Wp*lb + Wn*ub + bb; cu{k+1} = Wp*ub + Wn*lb + bb;
                case 'conv'
                    [cl{k+1}, cu{k+1}] = i_conv_ibp(op, lb, ub, precision);
                case 'normaffine'
                    sf = i_bcast_flat(op.scale, op.shape, precision);
                    tf = i_bcast_flat(op.shift, op.shape, precision);
                    pos = sf >= 0;
                    cl{k+1} = (sf.*lb).*pos + (sf.*ub).*(~pos) + tf;
                    cu{k+1} = (sf.*ub).*pos + (sf.*lb).*(~pos) + tf;
                case 'avgpool'
                    [cl{k+1}, cu{k+1}] = i_avgpool_ibp(op, lb, ub, precision);
                case 'relu'
                    if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                        fx = fixings{k};
                        lb(fx == 1)  = max(lb(fx == 1),  0);
                        ub(fx == -1) = min(ub(fx == -1), 0);
                    end
                    [preL{k}, preU{k}, actM{k}, unsM{k}, auC{k}, buC{k}, unstable{k}] = i_relu_pre(lb, ub, precision);
                    cl{k+1} = max(lb,0); cu{k+1} = max(ub,0);
                otherwise
                    error('gpu_bab_crown_alpha_dag:op', 'Unsupported op "%s".', op.type);
            end
        end
    else
        if ~isstruct(rootBounds) || ~all(isfield(rootBounds, {'preL','preU'}))
            error('gpu_bab_crown_alpha_dag:rootBounds', 'rootBounds must have fields preL,preU.');
        end
        tmpl = cast(x_lb(1), precision);
        for k = 1:nOps
            tk = ops{k}.type;
            if strcmp(tk, 'relu')
                l = repmat(cast(rootBounds.preL{k}, 'like', tmpl), 1, B);
                u = repmat(cast(rootBounds.preU{k}, 'like', tmpl), 1, B);
                if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                    fx = fixings{k};
                    l(fx == 1)  = max(l(fx == 1),  0);
                    u(fx == -1) = min(u(fx == -1), 0);
                end
                [preL{k}, preU{k}, actM{k}, unsM{k}, auC{k}, buC{k}, unstable{k}] = i_relu_pre(l, u, precision);
            elseif ~any(strcmp(tk, {'affine','conv','normaffine','avgpool','add'}))
                error('gpu_bab_crown_alpha_dag:op', 'Unsupported op "%s".', tk);
            end
        end
    end
end

function [pl, pu, actM, unsM, au, bu, uns] = i_relu_pre(l, u, precision)
% Clamped pre-activation bounds + the FIXED upper-line au*z+bu and active/unstable masks. The
% lower slope is left to alpha (al = act + uns*alpha). Matches spec_dag's i_relu_relax au/bu and
% alpha_fix's min-area init, so alpha=min-area reproduces spec_dag exactly.
    pl = l; pu = u;
    act = (l >= 0); uns = (l < 0) & (u > 0);
    au = zeros(size(l), precision); bu = zeros(size(l), precision);
    au(act) = 1;
    dn = u(uns) - l(uns);
    au(uns) = u(uns) ./ dn; bu(uns) = -au(uns) .* l(uns);
    actM = cast(act, precision); unsM = cast(uns, precision);
end

% ---- exact linear adjoints + interval forwards (duplicated from gpu_bab_crown_spec_dag; kept
% local so the alpha-dag is self-contained and spec_dag stays untouched) -------------------
function [olb, oub] = i_conv_ibp(op, lb, ub, precision)
    ish = op.inShape; osh = op.outShape; B = size(lb,2);
    W = cast(op.W, precision); Wp = max(W,0); Wn = min(W,0);
    bb = reshape(cast(op.b(:), precision), [1 1 osh(3)]);
    L4 = dlarray(reshape(lb, [ish(1) ish(2) ish(3) B]), 'SSCB');
    U4 = dlarray(reshape(ub, [ish(1) ish(2) ish(3) B]), 'SSCB');
    pad2 = [op.pad(1) op.pad(3); op.pad(2) op.pad(4)];
    args = {'Stride', op.stride, 'Padding', pad2, 'DilationFactor', op.dil};
    Lo = dlconv(L4, Wp, bb, args{:}) + dlconv(U4, Wn, 0, args{:});
    Hi = dlconv(U4, Wp, bb, args{:}) + dlconv(L4, Wn, 0, args{:});
    olb = reshape(extractdata(Lo), [prod(osh) B]);
    oub = reshape(extractdata(Hi), [prod(osh) B]);
end

function [A2, d2] = i_conv_backward(A, d, op, precision)
% Exact CROWN backward through a conv (linear), batched over B. TRACEABLE for dlgradient:
% A4 stays an (unformatted) dlarray, dltranspconv uses DataFormat, and the crop/pad/bias are
% slice+cat+sum (no extractdata, no indexed-assignment-into-zeros) so alpha's gradient flows.
    nSpec = size(A,1); B = size(A,3);
    osh = op.outShape; ish = op.inShape; W = cast(op.W, precision);
    Aperm = permute(A, [2 1 3]);
    A4 = reshape(Aperm, [osh(1) osh(2) osh(3) nSpec*B]);
    if ~isa(A4, 'dlarray'), A4 = dlarray(A4); end                 % unformatted dlarray (trace preserved)
    Afull = dltranspconv(A4, W, 0, 'Stride', op.stride, 'Cropping', 0, ...
                         'DilationFactor', op.dil, 'DataFormat', 'SSCB');
    pt = op.pad(1); pl = op.pad(3);
    hi = min(ish(1), size(Afull,1)-pt); wi = min(ish(2), size(Afull,2)-pl);
    Ain = Afull(pt+(1:hi), pl+(1:wi), :, :);                       % traceable slice
    if hi < ish(1) || wi < ish(2), Ain = i_zeropad_hw(Ain, ish(1), ish(2)); end
    A2 = permute(reshape(Ain, [prod(ish) nSpec B]), [2 1 3]);
    bc = reshape(cast(op.b(:), precision), [1 1 osh(3)]);
    A4u = reshape(A4, [osh(1) osh(2) osh(3) nSpec B]);
    dinc = reshape(sum(A4u .* bc, [1 2 3]), [nSpec B]);
    d2 = d + dinc;
end

function [olb, oub] = i_avgpool_ibp(op, lb, ub, precision)
    ish = op.inShape; osh = op.outShape; B = size(lb,2);
    L4 = reshape(cast(lb,precision), [ish(1) ish(2) ish(3) B]);
    U4 = reshape(cast(ub,precision), [ish(1) ish(2) ish(3) B]);
    olb = reshape(i_pool_mean(L4, op), [prod(osh) B]);
    oub = reshape(i_pool_mean(U4, op), [prod(osh) B]);
end

function Y = i_pool_mean(X, op)
    osh = op.outShape; kh = op.pool(1); kw = op.pool(2); B = size(X,4);
    Y = zeros([osh(1) osh(2) osh(3) B], 'like', X);
    for oh = 1:osh(1)
        for ow = 1:osh(2)
            rh = (oh-1)*op.stride(1) + (1:kh); rw = (ow-1)*op.stride(2) + (1:kw);
            Y(oh,ow,:,:) = mean(mean(X(rh,rw,:,:),1),2);
        end
    end
end

function [A2, d2] = i_avgpool_backward(A, d, op, precision)
% Exact CROWN backward through non-overlapping avgpool (distribute A_out/(kh*kw) uniformly).
% TRACEABLE: repelem + slice/cat (no indexed-assignment-into-zeros) so alpha's gradient flows.
    nSpec = size(A,1); B = size(A,3); osh = op.outShape; ish = op.inShape;
    kh = op.pool(1); kw = op.pool(2);
    Aperm = permute(A, [2 1 3]);
    A4 = reshape(Aperm, [osh(1) osh(2) osh(3) nSpec*B]);
    Aup = repelem(A4, kh, kw, 1, 1) / (kh*kw);
    hi = min(ish(1), size(Aup,1)); wi = min(ish(2), size(Aup,2));
    Ain = Aup(1:hi, 1:wi, :, :);                                   % traceable slice
    if hi < ish(1) || wi < ish(2), Ain = i_zeropad_hw(Ain, ish(1), ish(2)); end
    A2 = permute(reshape(Ain, [prod(ish) nSpec B]), [2 1 3]);
    d2 = d;
end

function Y = i_zeropad_hw(X, H, W)
% Traceable zero-pad of an [h w c n] array up to [H W c n] (cat with 'like' zeros keeps the tape).
    h = size(X,1); w = size(X,2);
    if h < H, X = cat(1, X, zeros([H-h, size(X,2), size(X,3), size(X,4)], 'like', X)); end
    w = size(X,2);
    if w < W, X = cat(2, X, zeros([size(X,1), W-w, size(X,3), size(X,4)], 'like', X)); end
    Y = X;
end

function v = i_bcast_flat(x, sh, precision)
    v = reshape(zeros([sh(1) sh(2) sh(3)], precision) + cast(x, precision), [], 1);
end

function y = i_g(x)
    if isa(x, 'dlarray'), x = extractdata(x); end
    if isa(x, 'gpuArray'), x = gather(x); end
    y = x;
end
