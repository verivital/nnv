function [margins, preL, preU, scoreCell] = gpu_bab_crown_spec_dag(ops, x_lb, x_ub, C, precision, fixings, rootBounds, alphaCell)
% AMORTIZED alpha-CROWN: optional alphaCell (cell(nOps,1) of per-relu lower-slope vectors,
% dim_k x 1, in [0,1]) lets the caller bound a whole BaB frontier with FIXED root-optimized
% slopes -- no autodiff, no per-node gradient tape -> large frontier. Unset -> min-area (default,
% bound-for-bound unchanged). Any alpha in [0,1] is a sound lower ReLU slope, so this is sound.
% GPU_BAB_CROWN_SPEC_DAG  Sound CROWN lower bound on a linear output spec C*f(x), batched
%   over B node columns, for feedforward conv nets INCLUDING residual DAGs
%   (affine/conv/normaffine/avgpool/relu/add). The generalisation of gpu_bab_crown_spec from
%   FC to conv: conv/BN/avgpool/add are LINEAR (exact interval forward + exact adjoint
%   backward), only ReLU is relaxed. The backward coefficient A is nSpec-by-dim-by-B (the SPEC
%   rows, nSpec ~ nClasses-1, NOT the layer width), so the memory is nSpec*HWC*B -- feasible
%   for conv, unlike batching the tight intermediate bounds' eye(nk) seed.
%
%   FULL DAG: the forward IBP caches each op's output bounds (so a residual 'add' fetches BOTH
%   its inputs and any op consuming a non-previous op bounds correctly) and the backward CROWN
%   routes each op's coefficient to op.src -- an 'add' routes it UNCHANGED to BOTH inputs. A
%   pure chain (every src==k-1, no 'add') reduces to the rolling forward + sequential backward
%   exactly (bound-for-bound unchanged). 'maxpool' is NOT handled here (per-node window
%   relaxation is not yet batched) and ERRORS -> sound-by-refusal (caller runs the tight path).
%
%   [margins, preL, preU] = GPU_BAB_CROWN_SPEC_DAG(ops, x_lb, x_ub, C, precision, fixings, rootBounds)
%     ops        : op list (nn_to_ops); affine/conv/normaffine/avgpool/relu/add (DAG via op.src,
%                  'add' via op.inputs=[a b]); NO 'maxpool'
%     x_lb,x_ub  : n-by-B input box columns (B = batch/node dim)
%     C          : nSpec-by-nOut spec
%     fixings    : optional cell(nOps,1) of per-relu dim_k-by-B node clamps (-1/0/+1)
%     rootBounds : optional struct with fields .preL,.preU (cell(nOps,1), dim_k-by-1 TIGHT
%                  pre-activation bounds computed ONCE at the root by gpu_bab_crown_tight).
%                  When given, the loose per-node IBP forward is SKIPPED and each node's
%                  pre-activation bounds are the (broadcast) root bounds clamped per node by
%                  the fixings -- tight bounds at batched speed (the #1 tightness lever).
%     margins    : nSpec-by-B lower bound on C*f over each node's clamped box
%     preL,preU  : per-relu (clamped) pre-activation bounds, dim_k-by-B
%
%   SOUNDNESS: interval-conv forward (Wp*lb+Wn*ub, the tightest linear interval) is sound; the
%   conv/avgpool/normaffine/add backward adjoints are EXACT (linear, no relaxation -- 'add'
%   routes the coefficient unchanged to both summands); the ReLU relaxation is the standard
%   sign-aware lower/upper line; the per-node clamps only tighten pre-activation bounds. With
%   rootBounds, the reused root bounds hold over the FULL input box -- a superset of every
%   node's sub-region -- and the per-neuron clamp (active: l>=0, inactive: u<=0) is the split's
%   domain restriction, so the bounds stay sound and are tighter than the per-node IBP. The DAG
%   backward is sound by induction: ops are topologically ordered, so every consumer of op k is
%   processed before op k, hence skipA{k} (the total coefficient on op k's output) is complete
%   when op k is bounded. For every x in node k's box, C*f(x) >= margins(:,k).

    if nargin < 5 || isempty(precision), precision = 'single'; end
    if nargin < 6, fixings = {}; end
    if nargin < 7, rootBounds = []; end
    if nargin < 8, alphaCell = {}; end
    B = size(x_lb, 2); nSpec = size(C, 1); nOps = numel(ops); n = size(x_lb, 1);
    preL = cell(nOps,1); preU = cell(nOps,1);

    if isempty(rootBounds)
        % ---- forward IBP (batched, FULL DAG), pre-activation bounds at each ReLU, with node
        % clamps. cl/cu cache each op's output bounds (index k+1 = op k; index 1 = op 0 = input)
        % so a residual 'add' fetches BOTH its inputs (out[a]+out[b], LINEAR -> exact) and any op
        % consuming a non-previous op (skip branch) bounds correctly. Sequential nets reduce to
        % the rolling bounds exactly. ----
        cl = cell(nOps+1,1); cu = cell(nOps+1,1);
        cl{1} = cast(x_lb, precision); cu{1} = cast(x_ub, precision);
        for k = 1:nOps
            op = ops{k};
            if strcmp(op.type, 'add')
                a = op.inputs(1)+1; b = op.inputs(2)+1;
                cl{k+1} = cl{a} + cl{b}; cu{k+1} = cu{a} + cu{b};
                continue;
            end
            if strcmp(op.type, 'concat')
                ic = op.inputs + 1;                       % stack input bounds (LINEAR -> exact)
                cl{k+1} = vertcat(cl{ic}); cu{k+1} = vertcat(cu{ic});
                continue;
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
                    preL{k} = lb; preU{k} = ub;
                    cl{k+1} = max(lb,0); cu{k+1} = max(ub,0);
                otherwise
                    error('gpu_bab_crown_spec_dag:op', ...
                        'Unsupported op "%s" (affine/conv/normaffine/avgpool/relu/add only).', op.type);
            end
        end
    else
        % ---- ROOT-TIGHT REUSE: skip the loose per-node IBP forward; build each ReLU's
        % pre-activation bounds from the TIGHT root bounds (broadcast to all B nodes), then
        % clamp per node by the fixings. The backward pass below only reads preL/preU at ReLU
        % ops, so no forward propagation is needed. SOUND (root bounds hold over the full box
        % >= each node's region; the own-neuron clamp is the split's domain restriction) and
        % much tighter than IBP. Non-relu ops (incl. 'add') carry no relu bounds; validate support.
        if ~isstruct(rootBounds) || ~all(isfield(rootBounds, {'preL','preU'}))
            error('gpu_bab_crown_spec_dag:rootBounds', ...
                'rootBounds must be a struct with fields preL and preU (per-op cell arrays).');
        end
        tmpl = cast(x_lb(1), precision); % scalar template for device+precision of the broadcasts (no n-by-B copy)
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
                preL{k} = l; preU{k} = u;
            elseif ~any(strcmp(tk, {'affine','conv','normaffine','avgpool','add','concat'}))
                error('gpu_bab_crown_spec_dag:op', ...
                    'Unsupported op "%s" (affine/conv/normaffine/avgpool/relu/add/concat only).', tk);
            end
        end
    end

    % ---- backward CROWN (batched over B, FULL DAG): lower bound on C*f ----
    % skipA{k} = accumulated nSpec-by-width_k-by-B coefficient on op k's OUTPUT. A single-input
    % op routes its coefficient (after its linear/relax backward) to op.src; an 'add' routes it
    % UNCHANGED to BOTH inputs (LINEAR -> exact). inputSkipA = total coefficient on op 0 (the
    % engine input). Topological order => every consumer of op k is visited before op k, so
    % skipA{k} is complete when op k is processed (sound by induction). For a pure chain this is
    % bound-for-bound the rolling pass (skipA{k-1} == the old running A).
    skipA = cell(nOps, 1);
    skipA{nOps} = repmat(cast(C, precision), 1, 1, B);   % spec sits on the network output (op nOps)
    d = zeros(nSpec, B, precision);
    inputSkipA = [];
    % BaBSR sensitivity score (computed only if the caller asks for the 4th output). At each ReLU,
    % neuron i's relaxation LOWERS the bound by exactly Aneg(s,i,b)*bu(i,b) for spec s; the total
    % slack a split there can recover is sum_s |Aneg(s,i,b)|*bu(i,b). Branching on the largest such
    % score (output-sensitivity x gap) closes the bound in far fewer nodes than largest-gap alone.
    wantScore = nargout >= 4;
    scoreCell = cell(nOps, 1);
    for k = nOps:-1:1
        A = skipA{k};
        if isempty(A), continue; end
        skipA{k} = [];                                   % free: keep only live branches resident
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
        if strcmp(op.type, 'concat')
            % out = [in_1; in_2; ...] -> SLICE the coefficient's dim axis (axis 2) back to each
            % input's block (transpose of stacking, LINEAR -> exact). A is nSpec x outDim x B.
            off = 0;
            for ii = 1:numel(op.inputs)
                sz = op.sizes(ii); s = op.inputs(ii);
                Ablk = A(:, off+(1:sz), :); off = off + sz;
                if s == 0
                    if isempty(inputSkipA), inputSkipA = Ablk; else, inputSkipA = inputSkipA + Ablk; end
                elseif isempty(skipA{s}), skipA{s} = Ablk;
                else, skipA{s} = skipA{s} + Ablk;
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
                l = preL{k}; u = preU{k};                  % dim x B
                ak = []; if ~isempty(alphaCell) && numel(alphaCell) >= k, ak = alphaCell{k}; end
                [au, bu, al] = i_relu_relax(l, u, precision, ak);
                dim = size(l, 1);
                Apos = max(A, 0); Aneg = min(A, 0);
                d = d + reshape(pagemtimes(Aneg, reshape(bu, dim, 1, B)), nSpec, B);
                if wantScore
                    % sum_s |Aneg(s,i,b)| * bu(i,b) -> dim x B (the BaBSR split score for this layer)
                    scoreCell{k} = reshape(sum(abs(Aneg), 1), dim, B) .* bu;
                end
                A = Apos .* reshape(al, 1, dim, B) + Aneg .* reshape(au, 1, dim, B);
            otherwise
                error('gpu_bab_crown_spec_dag:op', ...
                    'Unsupported op "%s" in backward (affine/conv/normaffine/avgpool/relu/add only).', op.type);
        end
        s = op.src;                                       % route coefficient to op.src
        if s == 0
            if isempty(inputSkipA), inputSkipA = A; else, inputSkipA = inputSkipA + A; end
        elseif isempty(skipA{s}), skipA{s} = A;
        else, skipA{s} = skipA{s} + A;
        end
    end

    if isempty(inputSkipA)
        margins = d;                                      % output independent of input (degenerate)
        return;
    end
    A = inputSkipA;
    Apos = max(A, 0); Aneg = min(A, 0);
    lbcol = reshape(cast(x_lb, precision), n, 1, B);
    ubcol = reshape(cast(x_ub, precision), n, 1, B);
    margins = reshape(pagemtimes(Apos, lbcol), nSpec, B) ...
            + reshape(pagemtimes(Aneg, ubcol), nSpec, B) + d;
end

% ---------------------------------------------------------------------------------------
function [olb, oub] = i_conv_ibp(op, lb, ub, precision)
% Interval conv (batched over B): the Wp/Wn split of the affine IBP with the matmul replaced
% by dlconv. Sound (tightest interval for the linear conv map). lb/ub are prod(inShape)-by-B.
    ish = op.inShape; osh = op.outShape; B = size(lb,2);
    W = cast(op.W, precision); Wp = max(W,0); Wn = min(W,0);
    bb = reshape(cast(op.b(:), precision), [1 1 osh(3)]);
    L4 = dlarray(reshape(lb, [ish(1) ish(2) ish(3) B]), 'SSCB');
    U4 = dlarray(reshape(ub, [ish(1) ish(2) ish(3) B]), 'SSCB');
    pad2 = [op.pad(1) op.pad(3); op.pad(2) op.pad(4)];      % [t l; b r]
    args = {'Stride', op.stride, 'Padding', pad2, 'DilationFactor', op.dil};
    Lo = dlconv(L4, Wp, bb, args{:}) + dlconv(U4, Wn, 0, args{:});
    Hi = dlconv(U4, Wp, bb, args{:}) + dlconv(L4, Wn, 0, args{:});
    olb = reshape(extractdata(Lo), [prod(osh) B]);
    oub = reshape(extractdata(Hi), [prod(osh) B]);
end

function [A2, d2] = i_conv_backward(A, d, op, precision)
% Exact CROWN backward through a conv (linear), batched over B: fold B into the
% dltranspconv batch dim (nSpec*B). A: nSpec x prod(outShape) x B.
    nSpec = size(A,1); B = size(A,3);
    osh = op.outShape; ish = op.inShape; W = cast(op.W, precision);
    Aperm = permute(A, [2 1 3]);                            % dim x nSpec x B
    A4 = dlarray(reshape(Aperm, [osh(1) osh(2) osh(3) nSpec*B]), 'SSCB');
    Afull = extractdata(dltranspconv(A4, W, 0, 'Stride', op.stride, 'Cropping', 0, 'DilationFactor', op.dil));
    pt = op.pad(1); pl = op.pad(3);
    Ain = zeros([ish(1) ish(2) ish(3) nSpec*B], 'like', Afull);
    hi = min(ish(1), size(Afull,1)-pt); wi = min(ish(2), size(Afull,2)-pl);
    if hi>0 && wi>0, Ain(1:hi,1:wi,:,:) = Afull(pt+(1:hi), pl+(1:wi), :, :); end
    A2 = reshape(Ain, [prod(ish) nSpec B]);                 % inDim x nSpec x B
    A2 = permute(A2, [2 1 3]);                              % nSpec x inDim x B
    bc = reshape(cast(op.b(:), precision), [1 1 osh(3)]);
    A4d = reshape(extractdata(A4), [osh(1) osh(2) osh(3) nSpec B]);
    dinc = reshape(sum(A4d .* bc, [1 2 3]), [nSpec B]);
    d2 = d + dinc;
end

function [olb, oub] = i_avgpool_ibp(op, lb, ub, precision)
% Interval avgpool (batched): mean over each non-overlapping window -> monotone, so the
% exact interval is the avgpool of the bounds. lb/ub prod(inShape)-by-B.
    ish = op.inShape; osh = op.outShape; B = size(lb,2);
    L4 = reshape(cast(lb,precision), [ish(1) ish(2) ish(3) B]);
    U4 = reshape(cast(ub,precision), [ish(1) ish(2) ish(3) B]);
    olb = reshape(i_pool_mean(L4, op), [prod(osh) B]);
    oub = reshape(i_pool_mean(U4, op), [prod(osh) B]);
end

function Y = i_pool_mean(X, op)
% Non-overlapping mean pool (stride==pool, unpadded) of X [H W C B].
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
% Exact CROWN backward through non-overlapping avgpool (batched over B): distribute
% A_out/(kh*kw) uniformly to each window's input cells. A: nSpec x prod(outShape) x B.
    nSpec = size(A,1); B = size(A,3); osh = op.outShape; ish = op.inShape;
    kh = op.pool(1); kw = op.pool(2);
    Aperm = permute(A, [2 1 3]);                            % dim x nSpec x B
    A4 = reshape(cast(Aperm,precision), [osh(1) osh(2) osh(3) nSpec*B]);
    Aup = repelem(A4, kh, kw, 1, 1) / (kh*kw);
    Ain = zeros([ish(1) ish(2) ish(3) nSpec*B], precision);
    hi = min(ish(1), size(Aup,1)); wi = min(ish(2), size(Aup,2));
    Ain(1:hi,1:wi,:,:) = Aup(1:hi,1:wi,:,:);
    A2 = permute(reshape(Ain, [prod(ish) nSpec B]), [2 1 3]);
    d2 = d;
end

function v = i_bcast_flat(x, sh, precision)
% Broadcast x ([1 1 C] / [H W C] / scalar) to a flat [prod(sh) x 1] column (column-major).
    v = reshape(zeros([sh(1) sh(2) sh(3)], precision) + cast(x, precision), [], 1);
end

function [au, bu, al] = i_relu_relax(l, u, precision, alpha)
% Per-neuron ReLU relaxation over [l,u] (elementwise, dim x B): stable-on (l>=0) identity,
% stable-off (u<=0) zero, unstable upper line au*z+bu / lower line al*z. al = min-area {0,1}
% by default; if alpha (dim x 1, the AMORTIZED root slopes) is given, al(unstable) = alpha
% (clamped [0,1] -> sound), broadcast over the B node columns.
    if nargin < 4, alpha = []; end
    au = zeros(size(l), precision); bu = zeros(size(l), precision); al = zeros(size(l), precision);
    act = (l >= 0); au(act) = 1; al(act) = 1;
    unst = (l < 0) & (u > 0); dn = u(unst) - l(unst);
    au(unst) = u(unst) ./ dn; bu(unst) = -au(unst) .* l(unst);
    if isempty(alpha)
        al(unst) = cast(u(unst) >= -l(unst), precision);     % min-area binary {0,1}
    else
        ab = repmat(cast(alpha(:), precision), 1, size(l,2));  % dim x 1 -> dim x B (root slopes)
        al(unst) = min(max(ab(unst), 0), 1);                   % clamp [0,1] (sound)
    end
end
