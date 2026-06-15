function [margins, preL, preU] = gpu_bab_crown_spec_dag(ops, x_lb, x_ub, C, precision, fixings, rootBounds)
% GPU_BAB_CROWN_SPEC_DAG  Sound CROWN lower bound on a linear output spec C*f(x), batched
%   over B node columns, for SEQUENTIAL conv nets (affine/conv/normaffine/avgpool/relu).
%   The generalisation of gpu_bab_crown_spec from FC to conv: conv/BN/avgpool are LINEAR
%   (exact interval forward + exact adjoint backward), only ReLU is relaxed. The backward
%   coefficient A is nSpec-by-dim-by-B (the SPEC rows, nSpec ~ nClasses-1, NOT the layer
%   width), so the memory is nSpec*HWC*B -- feasible for conv, unlike batching the tight
%   intermediate bounds' eye(nk) seed.
%
%   [margins, preL, preU] = GPU_BAB_CROWN_SPEC_DAG(ops, x_lb, x_ub, C, precision, fixings, rootBounds)
%     ops        : op list (nn_to_ops); SEQUENTIAL only (no 'add'/DAG, no 'maxpool')
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
%   SOUNDNESS: interval-conv forward (Wp*lb+Wn*ub, the tightest linear interval) is sound;
%   the conv/avgpool/normaffine backward adjoints are EXACT (linear, no relaxation); the
%   ReLU relaxation is the standard sign-aware lower/upper line; the per-node clamps only
%   tighten pre-activation bounds. With rootBounds, the reused root bounds hold over the FULL
%   input box -- a superset of every node's sub-region -- and the per-neuron clamp (active:
%   l>=0, inactive: u<=0) is the split's domain restriction, so the bounds stay sound and are
%   tighter than the per-node IBP. For every x in node k's box, C*f(x) >= margins(:,k).

    if nargin < 5 || isempty(precision), precision = 'single'; end
    if nargin < 6, fixings = {}; end
    if nargin < 7, rootBounds = []; end
    B = size(x_lb, 2); nSpec = size(C, 1); nOps = numel(ops);
    preL = cell(nOps,1); preU = cell(nOps,1);

    if isempty(rootBounds)
        % ---- forward IBP (batched), pre-activation bounds at each ReLU, with node clamps ----
        lb = cast(x_lb, precision); ub = cast(x_ub, precision);
        for k = 1:nOps
            op = ops{k};
            switch op.type
                case 'affine'
                    W = cast(op.W, precision); b = cast(op.b(:), precision);
                    Wp = max(W,0); Wn = min(W,0);
                    nlb = Wp*lb + Wn*ub + b; nub = Wp*ub + Wn*lb + b; lb = nlb; ub = nub;
                case 'conv'
                    [lb, ub] = i_conv_ibp(op, lb, ub, precision);
                case 'normaffine'
                    sf = i_bcast_flat(op.scale, op.shape, precision);
                    tf = i_bcast_flat(op.shift, op.shape, precision);
                    pos = sf >= 0;
                    nlb = (sf.*lb).*pos + (sf.*ub).*(~pos) + tf;
                    nub = (sf.*ub).*pos + (sf.*lb).*(~pos) + tf;
                    lb = nlb; ub = nub;
                case 'avgpool'
                    [lb, ub] = i_avgpool_ibp(op, lb, ub, precision);
                case 'relu'
                    if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                        fx = fixings{k};
                        lb(fx == 1)  = max(lb(fx == 1),  0);
                        ub(fx == -1) = min(ub(fx == -1), 0);
                    end
                    preL{k} = lb; preU{k} = ub;
                    lb = max(lb,0); ub = max(ub,0);
                otherwise
                    error('gpu_bab_crown_spec_dag:op', ...
                        'Unsupported op "%s" (sequential affine/conv/normaffine/avgpool/relu only).', op.type);
            end
        end
    else
        % ---- ROOT-TIGHT REUSE: skip the loose per-node IBP forward; build each ReLU's
        % pre-activation bounds from the TIGHT root bounds (broadcast to all B nodes), then
        % clamp per node by the fixings. The backward pass below only reads preL/preU at ReLU
        % ops, so no forward propagation is needed. SOUND (root bounds hold over the full box
        % >= each node's region; the own-neuron clamp is the split's domain restriction) and
        % much tighter than IBP. Non-relu ops are validated for support but carry no bounds.
        lb0 = cast(x_lb, precision);     % template fixing device + precision of the broadcasts
        for k = 1:nOps
            tk = ops{k}.type;
            if strcmp(tk, 'relu')
                l = repmat(cast(rootBounds.preL{k}, 'like', lb0), 1, B);
                u = repmat(cast(rootBounds.preU{k}, 'like', lb0), 1, B);
                if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                    fx = fixings{k};
                    l(fx == 1)  = max(l(fx == 1),  0);
                    u(fx == -1) = min(u(fx == -1), 0);
                end
                preL{k} = l; preU{k} = u;
            elseif ~any(strcmp(tk, {'affine','conv','normaffine','avgpool'}))
                error('gpu_bab_crown_spec_dag:op', ...
                    'Unsupported op "%s" (sequential affine/conv/normaffine/avgpool/relu only).', tk);
            end
        end
    end

    % ---- backward CROWN (batched over B): lower bound on C*f ----
    A = repmat(cast(C, precision), 1, 1, B);      % nSpec x nOut x B
    d = zeros(nSpec, B, precision);
    for k = nOps:-1:1
        op = ops{k};
        switch op.type
            case 'affine'
                W = cast(op.W, precision); b = cast(op.b(:), precision);
                d = d + reshape(pagemtimes(A, b), nSpec, B);
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
                [au, bu, al] = i_relu_relax(l, u, precision);
                dim = size(l, 1);
                Apos = max(A, 0); Aneg = min(A, 0);
                d = d + reshape(pagemtimes(Aneg, reshape(bu, dim, 1, B)), nSpec, B);
                A = Apos .* reshape(al, 1, dim, B) + Aneg .* reshape(au, 1, dim, B);
        end
    end

    n = size(x_lb, 1);
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

function [au, bu, al] = i_relu_relax(l, u, precision)
% Per-neuron ReLU relaxation over [l,u] (elementwise, dim x B): stable-on (l>=0) identity,
% stable-off (u<=0) zero, unstable upper line au*z+bu / lower line al*z (al in {0,1}).
    au = zeros(size(l), precision); bu = zeros(size(l), precision); al = zeros(size(l), precision);
    act = (l >= 0); au(act) = 1; al(act) = 1;
    unst = (l < 0) & (u > 0); dn = u(unst) - l(unst);
    au(unst) = u(unst) ./ dn; bu(unst) = -au(unst) .* l(unst);
    al(unst) = cast(u(unst) >= -l(unst), precision);
end
