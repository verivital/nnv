function [margins, preL, preU, unstable] = gpu_bab_crown_tight(ops, x_lb, x_ub, C, precision, fixings)
% GPU_BAB_CROWN_TIGHT  CROWN with TIGHT (backward) intermediate-layer bounds.
%
%   [margins, preL, preU] = GPU_BAB_CROWN_TIGHT(ops, x_lb, x_ub, C, precision)
%   computes a sound lower bound on C*f(x) over [x_lb,x_ub] using CROWN-tight
%   pre-activation bounds for EVERY ReLU layer -- each obtained by a backward CROWN
%   pass over the sub-network up to that layer, reusing the (already tight) bounds of
%   earlier layers. This is the O(L^2) core of real alpha-CROWN / auto_LiRPA and the fix
%   for the catastrophic looseness of IBP intermediate bounds on deep/wide nets
%   (measured: IBP gave a -40,000 margin where the true margin was +15).
%
%   Returns the final spec margins plus the tight preL/preU (so the ReLU-split BaB and
%   alpha-optimizer can reuse them instead of recomputing IBP bounds).
%
%   SOUNDNESS: every backward pass is a sound linear relaxation (sign-aware: a min uses
%   the lower line for positive coeffs / upper line for negative, a max the mirror), and
%   each layer's bound only uses SOUND bounds of strictly-earlier layers, so the
%   recursion is sound by induction. Single-box (the intermediate-bound recursion is
%   per input box); batching is a later optimization.

    if nargin < 5 || isempty(precision), precision = 'single'; end
    if nargin < 6, fixings = {}; end          % optional ReLU-split node fixings (-1/0/+1 per relu)
    nOps = numel(ops);
    preL = cell(nOps, 1);
    preU = cell(nOps, 1);
    unstable = cell(nOps, 1);

    % ---- tight intermediate bounds, layer by layer ----
    for k = 1:nOps
        if strcmp(ops{k}.type, 'relu')
            % z_k (the pre-activation feeding this ReLU) = output of ops[1..k-1].
            % Bound each of its neurons via a backward CROWN pass over ops[1..k-1],
            % which uses preL/preU of the strictly-earlier ReLUs (already computed).
            nk = i_layer_width(ops, k-1);
            Ck = eye(nk, precision);
            pu = i_backward(ops, k-1, Ck, x_lb, x_ub, preL, preU, precision, false);
            pl = i_backward(ops, k-1, Ck, x_lb, x_ub, preL, preU, precision, true);
            if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                fx = fixings{k};                    % BaB node: clamp fixed neurons + propagate
                pl(fx == 1)  = max(pl(fx == 1),  0);
                pu(fx == -1) = min(pu(fx == -1), 0);
            end
            preL{k} = pl; preU{k} = pu;
            unstable{k} = (pl < 0) & (pu > 0);
        end
    end

    % ---- final spec margin (lower bound on C*output) ----
    margins = i_backward(ops, nOps, cast(C, precision), x_lb, x_ub, preL, preU, precision, true);
end

function w = i_layer_width(ops, upto)
% # outputs of ops[1..upto] = rows of the last affine in that prefix.
    for k = upto:-1:1
        t = ops{k}.type;
        if strcmp(t, 'affine'),         w = size(ops{k}.W, 1);     return;
        elseif strcmp(t, 'conv'),       w = prod(ops{k}.outShape); return;
        elseif strcmp(t, 'normaffine'), w = prod(ops{k}.shape);    return;
        end
    end
    error('gpu_bab_crown_tight:nolinear', 'no affine/conv op before index %d', upto);
end

function bound = i_backward(ops, upto, A0, x_lb, x_ub, preL, preU, precision, lower)
% Backward CROWN over ops[1..upto] with initial coefficient A0, using preL/preU for the
% ReLU relaxations. lower=true -> lower bound (min); false -> upper bound (max).
    A = A0;
    nS = size(A0, 1);
    d = zeros(nS, 1, precision);
    for k = upto:-1:1
        op = ops{k};
        if strcmp(op.type, 'affine')
            W = cast(op.W, precision); b = cast(op.b(:), precision);
            d = d + A * b;
            A = A * W;
        elseif strcmp(op.type, 'conv')
            [A, d] = i_conv_backward(A, d, op, precision);
        elseif strcmp(op.type, 'normaffine')
            [A, d] = i_normaffine_backward(A, d, op, precision);
        else
            l = preL{k}; u = preU{k};
            [au, bu, al] = i_relax(l, u, precision);
            Apos = max(A, 0); Aneg = min(A, 0);
            if lower
                d = d + Aneg * bu;                 % min: +coeff->lower line(al), -coeff->upper line(au,bu)
                A = Apos .* al.' + Aneg .* au.';
            else
                d = d + Apos * bu;                 % max: +coeff->upper line(au,bu), -coeff->lower line(al)
                A = Apos .* au.' + Aneg .* al.';
            end
        end
    end
    Apos = max(A, 0); Aneg = min(A, 0);
    if lower
        bound = Apos * cast(x_lb, precision) + Aneg * cast(x_ub, precision) + d;
    else
        bound = Apos * cast(x_ub, precision) + Aneg * cast(x_lb, precision) + d;
    end
end

function [A2, d2] = i_conv_backward(A, d, op, precision)
% Exact CROWN backward through a conv (linear, NO relaxation): A_in = the conv ADJOINT
% (dltranspconv) of A_out, d += <A_out, b>. The adjoint identity <A_out,conv(W,x)> =
% <transpconv(A_out,W),x> holds EXACTLY -> sound (no over/under-approx). spec = batch dim.
    nS = size(A, 1);
    osh = op.outShape; ish = op.inShape;
    W = cast(op.W, precision);
    A4 = dlarray(reshape(A.', [osh(1) osh(2) osh(3) nS]), 'SSCB');
    % Exact conv adjoint: reconstruct the adjoint on the PADDED input (Cropping=0), then
    % crop to the original-input region (rows pad_t+1:pad_t+Hin, cols pad_l+1:pad_l+Win).
    % This is the exact transpose of "pad-then-conv" for any stride/pad/dilation (a
    % symmetric Cropping mis-sizes strided convs because the forward floor drops the
    % unused bottom pad). Tail rows the floor never reached carry zero coefficient -> 0.
    Afull = extractdata(dltranspconv(A4, W, 0, 'Stride', op.stride, 'Cropping', 0, 'DilationFactor', op.dil));
    pt = op.pad(1); pl = op.pad(3);
    Ain = zeros([ish(1) ish(2) ish(3) nS], 'like', Afull);
    hi = min(ish(1), size(Afull,1) - pt);
    wi = min(ish(2), size(Afull,2) - pl);
    if hi > 0 && wi > 0
        Ain(1:hi, 1:wi, :, :) = Afull(pt+(1:hi), pl+(1:wi), :, :);
    end
    A2 = reshape(Ain, [prod(ish) nS]).';                              % nS x prod(inShape)
    bc = reshape(cast(op.b(:), precision), [1 1 osh(3)]);
    dinc = squeeze(sum(extractdata(A4) .* bc, [1 2 3]));
    d2 = d + dinc(:);
end

function [A2, d2] = i_normaffine_backward(A, d, op, precision)
% Backward through y = s.*x + t (diagonal affine): A_in = A .* s_flat'; d += A * t_flat.
    sh = op.shape;
    sf = i_bcast_flat(op.scale, sh, precision);
    tf = i_bcast_flat(op.shift, sh, precision);
    A2 = A .* sf.';
    d2 = d + A * tf;
end

function v = i_bcast_flat(x, sh, precision)
% Broadcast x ([1 1 C] / [H W C] / scalar) to a flat [prod(sh) x 1] column (column-major).
    v = reshape(zeros([sh(1) sh(2) sh(3)], precision) + cast(x, precision), [], 1);
end

function [au, bu, al] = i_relax(l, u, precision)
% Per-neuron ReLU relaxation over [l,u]: stable-on (l>=0)->identity, stable-off (u<=0)->0,
% unstable-> upper line au*z+bu (au=u/(u-l), bu=-au*l>=0), lower line al*z (al min-area).
    m = numel(l);
    au = zeros(m, 1, precision); bu = zeros(m, 1, precision); al = zeros(m, 1, precision);
    act = (l >= 0); au(act) = 1; al(act) = 1;
    uns = (l < 0) & (u > 0);
    dn = u(uns) - l(uns);
    au(uns) = u(uns) ./ dn;
    bu(uns) = -au(uns) .* l(uns);
    al(uns) = cast(u(uns) >= -l(uns), precision);
end
