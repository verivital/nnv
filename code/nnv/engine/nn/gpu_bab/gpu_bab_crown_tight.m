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
    % Compute input bounds for every op whose backward relaxation needs them: ReLU (the
    % pre-activation) AND maxpool (the window inputs that decide the sound max relaxation).
    for k = 1:nOps
        tk = ops{k}.type;
        if strcmp(tk, 'relu') || strcmp(tk, 'maxpool')
            % z_k (the input feeding this op) = output of op `src` = ops{k}.src (the DAG input,
            % not assumed k-1). Bound each of its elements via a backward DAG-CROWN pass over
            % ops[1..src], reusing preL/preU of strictly-earlier ReLUs/maxpools (sound by induction).
            src = ops{k}.src;
            if src == 0, nk = numel(x_lb); else, nk = i_layer_width(ops, src); end
            Ck = eye(nk, precision);
            pu = i_backward(ops, src, Ck, x_lb, x_ub, preL, preU, precision, false);
            pl = i_backward(ops, src, Ck, x_lb, x_ub, preL, preU, precision, true);
            if strcmp(tk, 'relu') && ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                fx = fixings{k};                    % BaB node: clamp fixed neurons + propagate
                pl(fx == 1)  = max(pl(fx == 1),  0);
                pu(fx == -1) = min(pu(fx == -1), 0);
            end
            preL{k} = pl; preU{k} = pu;
            if strcmp(tk, 'relu')
                unstable{k} = (pl < 0) & (pu > 0);
            end
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
        elseif strcmp(t, 'avgpool'),    w = prod(ops{k}.outShape); return;
        elseif strcmp(t, 'maxpool'),    w = prod(ops{k}.outShape); return;
        elseif strcmp(t, 'add'),        w = prod(ops{k}.shape);    return;
        elseif strcmp(t, 'normaffine'), w = prod(ops{k}.shape);    return;
        end
    end
    error('gpu_bab_crown_tight:nolinear', 'no affine/conv op before index %d', upto);
end

function bound = i_backward(ops, upto, A0, x_lb, x_ub, preL, preU, precision, lower)
% Backward CROWN over the DAG ops[1..upto] with initial coefficient A0 (on op `upto`'s output).
% FULL DAG: each op routes its backward coefficient to op.src (its input op), accumulated in
% skipA{src}; an 'add' routes UNCHANGED to BOTH inputs (linear -> exact). skipA{k} = accumulated
% coefficient on op k's OUTPUT; inputSkipA = coefficient on op 0 (the engine input). Ops are
% topologically ordered, so k=upto:-1:1 visits every consumer of op k before op k itself, so
% skipA{k} is complete when op k is processed (sound by induction). lower=true -> lower bound.
    nS = size(A0, 1);
    d = zeros(nS, 1, precision);
    if upto == 0                              % A0 is already on the engine input (op 0)
        Apos = max(A0, 0); Aneg = min(A0, 0);
        if lower, bound = Apos*cast(x_lb,precision) + Aneg*cast(x_ub,precision);
        else,     bound = Apos*cast(x_ub,precision) + Aneg*cast(x_lb,precision); end
        return;
    end
    skipA = cell(upto, 1);
    skipA{upto} = A0;                         % seed: A0 is the coefficient on op `upto`'s output
    inputSkipA = 0;
    for k = upto:-1:1
        A = skipA{k};
        if isempty(A), continue; end          % nothing routed to op k (dead w.r.t. the output)
        op = ops{k};
        if strcmp(op.type, 'add')             % out = out[a]+out[b]: route UNCHANGED to both
            for ii = 1:numel(op.inputs)
                s = op.inputs(ii);
                if s == 0,                inputSkipA = inputSkipA + A;
                elseif isempty(skipA{s}), skipA{s} = A;
                else,                     skipA{s} = skipA{s} + A;
                end
            end
            continue;
        end
        % single-input op: A (on op's OUTPUT) -> A (on op's INPUT)
        if strcmp(op.type, 'affine')
            W = cast(op.W, precision); b = cast(op.b(:), precision);
            d = d + A * b;
            A = A * W;
        elseif strcmp(op.type, 'conv')
            [A, d] = i_conv_backward(A, d, op, precision);
        elseif strcmp(op.type, 'normaffine')
            [A, d] = i_normaffine_backward(A, d, op, precision);
        elseif strcmp(op.type, 'avgpool')
            [A, d] = i_avgpool_backward(A, d, op, precision);
        elseif strcmp(op.type, 'maxpool')
            [A, d] = i_maxpool_backward(A, d, op, preL{k}, preU{k}, precision, lower);
        else                                  % relu relaxation (sign-aware), preL/preU{k}
            l = preL{k}; u = preU{k};
            [au, bu, al] = i_relax(l, u, precision);
            Apos = max(A, 0); Aneg = min(A, 0);
            if lower
                d = d + Aneg * bu;
                A = Apos .* al.' + Aneg .* au.';
            else
                d = d + Apos * bu;
                A = Apos .* au.' + Aneg .* al.';
            end
        end
        s = op.src;                           % route the input coefficient to op.src
        if s == 0,                inputSkipA = inputSkipA + A;
        elseif isempty(skipA{s}), skipA{s} = A;
        else,                     skipA{s} = skipA{s} + A;
        end
    end
    A = inputSkipA;                           % total coefficient on the engine input
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

function [A2, d2] = i_avgpool_backward(A, d, op, precision)
% Exact CROWN backward through avgpool (LINEAR, no relaxation): each output cell is the
% mean of its kh x kw window, so the adjoint distributes A_out/(kh*kw) uniformly back to
% the window's input cells. Non-overlapping (stride==pool) -> every input cell lies in at
% most one window -> the adjoint is repelem(A_out, kh, kw)/(kh*kw), zero-padded to inShape
% (floor-dropped tail cells were in no window -> zero coefficient). No bias -> d unchanged.
    nS = size(A, 1);
    osh = op.outShape; ish = op.inShape;
    kh = op.pool(1); kw = op.pool(2);
    A4 = reshape(A.', [osh(1) osh(2) osh(3) nS]);           % [Hout Wout C nS]
    Aup = repelem(cast(A4, precision), kh, kw, 1, 1) / (kh*kw);   % distribute to window cells
    Ain = zeros([ish(1) ish(2) ish(3) nS], precision);
    hi = min(ish(1), size(Aup,1)); wi = min(ish(2), size(Aup,2));
    Ain(1:hi, 1:wi, :, :) = Aup(1:hi, 1:wi, :, :);
    A2 = reshape(Ain, [prod(ish) nS]).';                   % nS x prod(inShape)
    d2 = d;
end

function [A2, d2] = i_maxpool_backward(A, d, op, l, u, precision, lower)
% Sound CROWN backward through max pooling. For each output cell o = max over its window:
%   lower bound  y_o >= x_m        (m = argmax of the window's LOWER bounds; max >= any elem)
%   upper bound  y_o <= x_m        if m is DECIDED (l_m >= u_i for every other i) -- EXACT
%               y_o <= max_i u_i  otherwise (a constant; sound since max <= max of uppers)
% Lmat/Umat are n_out x n_in selection matrices (sparse); overlapping windows accumulate
% naturally in Apos*Lmat. Sign-aware: +coeff uses the relaxation that bounds the spec on the
% correct side. Undecided cells leave a constant gap (the BaB can split the feeding ReLUs).
    [mIdx, decided, umax] = i_maxpool_relax(op, l, u);
    n_in = prod(op.inShape); n_out = prod(op.outShape);
    Lmat = sparse(1:n_out, mIdx, 1, n_out, n_in);               % y_o >= x_{m_o}
    dec  = find(decided);
    Umat = sparse(dec, mIdx(dec), 1, n_out, n_in);              % decided: y_o <= x_{m_o}
    Ubias = zeros(n_out, 1); Ubias(~decided) = umax(~decided);  % undecided: y_o <= umax
    Ad = double(A); Apos = max(Ad, 0); Aneg = min(Ad, 0);
    if lower
        A2 = cast(Apos * Lmat + Aneg * Umat, precision);
        d2 = d + cast(Aneg * Ubias, precision);
    else
        A2 = cast(Apos * Umat + Aneg * Lmat, precision);
        d2 = d + cast(Apos * Ubias, precision);
    end
end

function [mIdx, decided, umax] = i_maxpool_relax(op, l, u)
% Per output cell: m = flat input index of the window's argmax-LOWER element; decided = that
% element dominates (its lower bound >= every other window upper => max is exactly it); umax =
% max window upper. Unpadded windows (enforced at extraction); overlap allowed.
    ish = op.inShape; osh = op.outShape;
    H = ish(1); W = ish(2); C = ish(3); Ho = osh(1); Wo = osh(2);
    kh = op.pool(1); kw = op.pool(2); sh = op.stride(1); sw = op.stride(2);
    L = reshape(double(l), [H W C]); U = reshape(double(u), [H W C]);
    n_out = prod(osh);
    mIdx = ones(n_out, 1); decided = false(n_out, 1); umax = zeros(n_out, 1);
    for c = 1:C
        for ow = 1:Wo
            rw = (ow-1)*sw + (1:kw);
            for oh = 1:Ho
                rh = (oh-1)*sh + (1:kh);
                lw = L(rh, rw, c); uw = U(rh, rw, c);
                [maxl, p] = max(lw(:));
                [pr, pc] = ind2sub([kh kw], p);
                o = sub2ind([Ho Wo C], oh, ow, c);
                mIdx(o) = sub2ind([H W C], rh(pr), rw(pc), c);
                uw2 = uw; uw2(p) = -inf;
                decided(o) = maxl >= max(uw2(:));               % argmax-lower dominates others
                umax(o) = max(uw(:));
            end
        end
    end
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
