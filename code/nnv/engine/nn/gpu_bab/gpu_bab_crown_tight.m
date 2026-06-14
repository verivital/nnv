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
    w = [];
    for k = upto:-1:1
        if strcmp(ops{k}.type, 'affine'), w = size(ops{k}.W, 1); return; end
    end
    error('gpu_bab_crown_tight:noaffine', 'no affine op before index %d', upto);
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
