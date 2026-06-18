function [margins, preL, preU, unstable, Ain, din, mulPlanes] = gpu_bab_crown_tight(ops, x_lb, x_ub, C, precision, fixings, mulFix)
% GPU_BAB_CROWN_TIGHT  CROWN with TIGHT (backward) intermediate-layer bounds.
%   Optional 5th/6th outputs Ain (nSpec x nIn), din (nSpec x 1) are the CROWN input-space LOWER
%   plane for the spec: C*f(x) >= Ain*x + din for every x in [x_lb,x_ub] (margins = min over the
%   box of that plane). Used by Clip-and-Verify domain clipping (gpu_bab_clip) for the joint
%   output-polytope avoidance that single-halfspace separation misses.
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
    if nargin < 7, mulFix = {}; end           % optional per-'product' input-range overrides (struct .lo/.hi, stacked [in1;in2]) for targeted product-input branching
    nOps = numel(ops);
    preL = cell(nOps, 1);
    preU = cell(nOps, 1);
    unstable = cell(nOps, 1);
    mulPlanes = cell(nOps, 1);   % per 'product' op: input value-planes (for targeted product-input branching)

    % SOUND-FP32 (M2): per-op output value-magnitude majorants, used by i_backward's derr to bound
    % the FP32 backward roundoff. Single precision only (FP64 rounding ~1e-16 is negligible -> derr=0,
    % the FP64 path stays the oracle). On any IBP failure -> vmag={} -> derr=0 (the single path falls
    % back to UNSOUND-screen behaviour, which is fine because it never emits a verdict; the sound EMIT
    % path requires vmag to succeed, gated in gpu_bab_try_verify). See FP32_SOUND_RECIPE.
    vmag = {};
    if strcmp(precision, 'single')
        try
            [~, ~, vmag] = gpu_bab_ibp(ops, x_lb, x_ub, precision);
        catch
            vmag = {};
        end
    end

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
            pu = i_backward(ops, src, Ck, x_lb, x_ub, preL, preU, precision, false, vmag);
            pl = i_backward(ops, src, Ck, x_lb, x_ub, preL, preU, precision, true, vmag);
            if strcmp(tk, 'relu') && ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                fx = fixings{k};                    % BaB node: clamp fixed neurons + propagate
                pl(fx == 1)  = max(pl(fx == 1),  0);
                pu(fx == -1) = min(pu(fx == -1), 0);
            end
            preL{k} = pl; preU{k} = pu;
            if strcmp(tk, 'relu')
                unstable{k} = (pl < 0) & (pu > 0);
            end
        elseif strcmp(tk, 'product')
            % BILINEAR product feeds the McCormick relaxation from the OUTPUT bounds of BOTH its
            % input ops. Bound each input op's output (backward CROWN over its prefix, sound by
            % induction) and store them STACKED [in1; in2] under preL{k}/preU{k} (numeric, so the
            % BaB infeasible/gather logic stays valid). Width per input = ops{k}.sizes(1)=sizes(2).
            ins = ops{k}.inputs;
            plS = []; puS = [];
            mp = struct('inputs', ins, 'sizes', ops{k}.sizes, ...
                        'AL', {cell(1,2)}, 'dL', {cell(1,2)}, 'AU', {cell(1,2)}, 'dU', {cell(1,2)});
            for ii = 1:2
                si = ins(ii);
                if si == 0, nki = numel(x_lb); else, nki = i_layer_width(ops, si); end
                Cki = eye(nki, precision);
                [pui, AinU, dinU] = i_backward(ops, si, Cki, x_lb, x_ub, preL, preU, precision, false, vmag);
                [pli, AinL, dinL] = i_backward(ops, si, Cki, x_lb, x_ub, preL, preU, precision, true, vmag);
                plS = [plS; pli]; puS = [puS; pui]; %#ok<AGROW>
                % input value-planes over the box: AL*x+dL <= value <= AU*x+dU (per input element).
                % Used by the BaB to branch a product input's range via gpu_bab_clip (constraint
                % value<=mid  =>  AL*x + (dL-mid) <= 0 ;  value>=mid  =>  -AU*x + (mid-dU) <= 0).
                mp.AL{ii} = AinL; mp.dL{ii} = dinL; mp.AU{ii} = AinU; mp.dU{ii} = dinU;
            end
            preL{k} = plS; preU{k} = puS;
            mulPlanes{k} = mp;
            % TARGETED BRANCHING override: clamp a product input's value-range to a branched
            % sub-interval. SOUND: the two children partition by v(x)<=mid / v(x)>=mid (every input
            % falls in one half), so bounding the spec with v clamped to [lo,hi] is a valid lower
            % bound for the inputs whose v lies there; the McCormick over the narrower range is
            % tighter, so the gap shrinks geometrically per split.
            if ~isempty(mulFix) && numel(mulFix) >= k && ~isempty(mulFix{k})
                preL{k} = max(preL{k}, cast(mulFix{k}.lo(:), precision));
                preU{k} = min(preU{k}, cast(mulFix{k}.hi(:), precision));
            end
        end
    end

    % ---- final spec margin (lower bound on C*output) + the input-space lower plane ----
    [margins, Ain, din] = i_backward(ops, nOps, cast(C, precision), x_lb, x_ub, preL, preU, precision, true, vmag);
end

function w = i_layer_width(ops, upto)
% Flat output width of op `upto`. Prefer the recorded nOut (exact for ALL ops incl. flat
% add/normaffine, where prod(shape)=prod([])=1 was a latent width=1 bug); fall back to the
% structural walk for ops emitted without it.
    if isfield(ops{upto}, 'nOut') && ~isempty(ops{upto}.nOut)
        w = ops{upto}.nOut; return;
    end
    for k = upto:-1:1
        t = ops{k}.type;
        if strcmp(t, 'affine'),         w = size(ops{k}.W, 1);     return;
        elseif strcmp(t, 'conv'),       w = prod(ops{k}.outShape); return;
        elseif strcmp(t, 'avgpool'),    w = prod(ops{k}.outShape); return;
        elseif strcmp(t, 'maxpool'),    w = prod(ops{k}.outShape); return;
        elseif strcmp(t, 'add'),        w = prod(ops{k}.shape);    return;
        elseif strcmp(t, 'normaffine'), w = prod(ops{k}.shape);    return;
        elseif strcmp(t, 'concat'),     w = sum(ops{k}.sizes);     return;
        elseif strcmp(t, 'product'),    w = ops{k}.sizes(1);       return;  % elementwise: out width = each input width
        end
    end
    error('gpu_bab_crown_tight:nolinear', 'no affine/conv op before index %d', upto);
end

function rad = i_outward_rad(A, x_lb, x_ub, precision)
% SOUND-FP32 outward error radius (Higham running-error). At a CROWN concretization
% bound = Apos*x_lb + Aneg*x_ub + d (an nS-vector, computed in `precision` over the box), the
% FP rounding error of the length-n input contraction is bounded by gamma_n*(|A|*|x_mag|) with
% gamma_n = n*u/(1-n*u) and u the unit roundoff. Returning this nS x 1 radius (a TRANSIENT vector,
% never an A-shaped tensor -> memory-flat, fits 11 GB) lets the caller widen OUTWARD: lower bounds
% -= rad, upper += rad. SOUND: deterministic worst-case gamma_n (never the probabilistic sqrt(n)u),
% pre-inflated by k (>=2) so the radius GEMM's own rounding can't make rad an under-estimate, plus
% n*realmin to absorb subnormal flush. Only fires for 'single' (FP64 rounding ~1e-16 is negligible
% and the FP64 path stays the oracle). See research/FP32_SOUND_RECIPE_2026-06-17.md.
%   NOTE (M1): k is deliberately huge (10) and only the final contraction is widened -- this proves
%   the outward plumbing + the speed/memory envelope; full soundness (per-op backward error, Edits
%   B-D) follows in M2. NOT a certified EMIT path until M2 + the G1 parity gate pass.
    if ~strcmp(precision, 'single'), rad = zeros(size(A,1), 1, 'like', A); return; end
    n = numel(x_lb);
    u = single(eps('single') / 2);
    k = single(2);                                    % M2: 2x pre-inflation (the per-op derr now covers the backward error)
    gbar = k * (n * u) / (1 - n * u);
    xmag = max(abs(single(x_lb(:))), abs(single(x_ub(:))));
    rad = gbar * (abs(A) * xmag) + single(n) * realmin('single');
end

function g = i_gamma(m, precision)
% Pre-inflated (2x) deterministic worst-case Higham factor for a length-m FP contraction (the 2x
% absorbs the derr computation's own rounding). gamma_m = m*u/(1-m*u); m over-estimated where cheap.
% SOUND only with the DETERMINISTIC m (never the probabilistic sqrt(m) -- a single tail = a -150).
    u = single(eps('single') / 2);
    m = single(m);
    g = cast(2, precision) * (m * u) / (1 - m * u);
end

function [bound, Ain, din] = i_backward(ops, upto, A0, x_lb, x_ub, preL, preU, precision, lower, vmag)
% Backward CROWN over the DAG ops[1..upto] with initial coefficient A0 (on op `upto`'s output).
% FULL DAG: each op routes its backward coefficient to op.src (its input op), accumulated in
% skipA{src}; an 'add' routes UNCHANGED to BOTH inputs (linear -> exact). skipA{k} = accumulated
% coefficient on op k's OUTPUT; inputSkipA = coefficient on op 0 (the engine input). Ops are
% topologically ordered, so k=upto:-1:1 visits every consumer of op k before op k itself, so
% skipA{k} is complete when op k is processed (sound by induction). lower=true -> lower bound.
% Optional Ain (nS x nIn), din (nS x 1): the input-space affine form, so for lower=true the
% bounded quantity >= Ain*x + din for all x in the box (bound = min over the box of that plane).
%
% SOUND-FP32 (M2): when vmag is supplied (single precision), accumulate derr (nS x 1) = a sound
% bound on the FP32 ROUNDING error each op's backward arithmetic injects into the final bound. Each
% op's matmul/relaxation roundoff <= gamma_m * |coeff| * |operand|; contracted RIGHT HERE with the
% op input's value-magnitude majorant vmag{src+1}, it collapses to nS x 1 (memory-flat). The matmul
% ops (affine/conv) use the |W|-amplified magnitude (|A_in| understates under cancellation); the
% elementwise/monotone ops (normaffine/avgpool/relu/maxpool) use |A_in| (cancellation-free), plus
% the relaxation INTERCEPT (relu bu / maxpool umax) error. The final bound widens OUTWARD by
% (rad + derr). See research/FP32_SOUND_RECIPE_2026-06-17.md.
    if nargin < 10, vmag = {}; end
    nS = size(A0, 1);
    d = zeros(nS, 1, precision);
    doErr = strcmp(precision, 'single') && ~isempty(vmag);
    derr = zeros(nS, 1, precision);
    if upto == 0                              % A0 is already on the engine input (op 0)
        Apos = max(A0, 0); Aneg = min(A0, 0);
        if lower, bound = Apos*cast(x_lb,precision) + Aneg*cast(x_ub,precision);
        else,     bound = Apos*cast(x_ub,precision) + Aneg*cast(x_lb,precision); end
        rad = i_outward_rad(A0, x_lb, x_ub, precision);      % sound-FP32: widen OUTWARD
        if lower, bound = bound - rad; else, bound = bound + rad; end
        Ain = A0; din = zeros(nS, 1, precision);
        return;
    end
    skipA = cell(upto, 1);
    skipA{upto} = A0;                         % seed: A0 is the coefficient on op `upto`'s output
    inputSkipA = zeros(nS, numel(x_lb), precision);   % coefficient on the engine input (op 0); always nS x nIn
    for k = upto:-1:1
        A = skipA{k};
        if isempty(A), continue; end          % nothing routed to op k (dead w.r.t. the output)
        op = ops{k};
        if strcmp(op.type, 'add')             % out = out[a]+out[b]: route UNCHANGED to both
            for ii = 1:numel(op.inputs)
                s = op.inputs(ii);
                if doErr, derr = derr + i_gamma(1, precision) * (abs(A) * vmag{s + 1}); end   % accumulation roundoff (conservative)
                if s == 0,                inputSkipA = inputSkipA + A;
                elseif isempty(skipA{s}), skipA{s} = A;
                else,                     skipA{s} = skipA{s} + A;
                end
            end
            continue;
        end
        if strcmp(op.type, 'concat')          % out = [in_1; in_2; ...]: SLICE columns back to each input
            off = 0;
            for ii = 1:numel(op.inputs)
                sz = op.sizes(ii); s = op.inputs(ii);
                Ablk = A(:, off+(1:sz)); off = off + sz;
                if s == 0,                inputSkipA = inputSkipA + Ablk;
                elseif isempty(skipA{s}), skipA{s} = Ablk;
                else,                     skipA{s} = skipA{s} + Ablk;
                end
            end
            continue;
        end
        if strcmp(op.type, 'product')         % out = in[a].*in[b], bilinear: McCormick, route to BOTH inputs
            wa = op.sizes(1); lz = preL{k}; uz = preU{k};
            la = lz(1:wa);       ua = uz(1:wa);          % input a (x) output bounds
            lyy = lz(wa+1:end);  uyy = uz(wa+1:end);     % input b (y) output bounds
            [aL,bL,cL,aU,bU,cU] = gpu_bab_mul_relax(la, ua, lyy, uyy, [], [], precision);
            Apos = max(A, 0); Aneg = min(A, 0);          % sign-aware: +coeff -> LOWER plane (lower bnd)
            if doErr   % McCormick plane mults to both inputs + the cL/cU intercept (conservative)
                va = vmag{op.inputs(1) + 1}; vb = vmag{op.inputs(2) + 1};
                derr = derr + i_gamma(2, precision) * (abs(A) * (va + vb)) ...
                            + i_gamma(wa, precision) * (abs(A) * (abs(cL) + abs(cU)));
            end
            if lower
                Ax = Apos .* aL.' + Aneg .* aU.';
                Ay = Apos .* bL.' + Aneg .* bU.';
                d  = d + Apos * cL + Aneg * cU;
            else
                Ax = Apos .* aU.' + Aneg .* aL.';
                Ay = Apos .* bU.' + Aneg .* bL.';
                d  = d + Apos * cU + Aneg * cL;
            end
            sa = op.inputs(1);
            if sa == 0,                inputSkipA = inputSkipA + Ax;
            elseif isempty(skipA{sa}), skipA{sa} = Ax;
            else,                      skipA{sa} = skipA{sa} + Ax;
            end
            sb = op.inputs(2);
            if sb == 0,                inputSkipA = inputSkipA + Ay;
            elseif isempty(skipA{sb}), skipA{sb} = Ay;
            else,                      skipA{sb} = skipA{sb} + Ay;
            end
            continue;
        end
        % single-input op: A (on op's OUTPUT) -> A (on op's INPUT)
        if doErr, vin = vmag{op.src + 1}; end           % input value-magnitude majorant (op 0 = input box)
        if strcmp(op.type, 'affine')
            W = cast(op.W, precision); b = cast(op.b(:), precision);
            if doErr   % A*W contracts over out-width; output value-mag majorant (no cancel) = |W|*vin+|b|
                derr = derr + i_gamma(size(A,2), precision) * (abs(A) * (abs(W) * vin + abs(b)));
            end
            d = d + A * b;
            A = A * W;
        elseif strcmp(op.type, 'conv')
            if doErr   % magnitude adjoint transpconv(|A|,|W|) + |b| pass; m = kh*kw*out-ch
                opm = op; opm.W = abs(op.W); opm.b = abs(op.b);
                [Amag, dmag] = i_conv_backward(abs(A), zeros(nS,1,precision), opm, precision);
                mC = size(op.W,1) * size(op.W,2) * size(op.W,4);
                derr = derr + i_gamma(mC, precision) * (Amag * vin + dmag);
            end
            [A, d] = i_conv_backward(A, d, op, precision);
        elseif strcmp(op.type, 'normaffine')
            [A, d] = i_normaffine_backward(A, d, op, precision);
            if doErr, derr = derr + i_gamma(2, precision) * (abs(A) * vin); end   % A now = A_in (elementwise, no cancel)
        elseif strcmp(op.type, 'avgpool')
            [A, d] = i_avgpool_backward(A, d, op, precision);
            if doErr, derr = derr + i_gamma(prod(op.pool), precision) * (abs(A) * vin); end
        elseif strcmp(op.type, 'maxpool')
            [A, d] = i_maxpool_backward(A, d, op, preL{k}, preU{k}, precision, lower);
            if doErr   % selection exact (0/1); conservative term for the umax relaxation intercept
                derr = derr + i_gamma(2*prod(op.pool), precision) * (abs(A) * (vin + vmag{k + 1}));
            end
        else                                  % relu relaxation (sign-aware), preL/preU{k}
            l = preL{k}; u = preU{k};
            [au, bu, al] = i_relax(l, u, precision);
            Apos = max(A, 0); Aneg = min(A, 0);
            if doErr   % slope mults (|A_in|<=|A|, no cancel) + the bu intercept (d += Aneg*bu)
                derr = derr + i_gamma(2, precision) * (abs(A) * vin) ...
                            + i_gamma(numel(bu), precision) * (abs(Aneg) * abs(bu));
            end
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
    A = inputSkipA;                           % total coefficient on the engine input (nS x nIn; zero if no input dependence)
    Ain = A; din = d;                         % input-space affine form: bounded >= Ain*x + din (lower)
    Apos = max(A, 0); Aneg = min(A, 0);
    if lower
        bound = Apos * cast(x_lb, precision) + Aneg * cast(x_ub, precision) + d;
    else
        bound = Apos * cast(x_ub, precision) + Aneg * cast(x_lb, precision) + d;
    end
    rad = i_outward_rad(A, x_lb, x_ub, precision);           % final-contraction roundoff (M1)
    if lower, bound = bound - rad - derr; else, bound = bound + rad + derr; end   % widen OUTWARD by rad + the per-op backward error (M2)
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
