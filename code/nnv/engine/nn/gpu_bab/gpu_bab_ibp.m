function [out_lb, out_ub] = gpu_bab_ibp(ops, in_lb, in_ub, precision)
% GPU_BAB_IBP  Sound interval bound propagation for a feedforward ReLU net.
%
%   [out_lb, out_ub] = GPU_BAB_IBP(ops, in_lb, in_ub, precision) forward-
%   propagates the input box [in_lb, in_ub] through the op list to SOUND output
%   bounds. This is the LP-free, GPU-native foundation of NNV's GPU-BaB engine
%   (Phase 1): pure linear algebra (no linprog), so the *identical* code runs on
%   host arrays or gpuArray, and the bound columns can carry a BATCH dimension
%   (the GPU-BaB parallel dimension over sub-domains). It is intentionally the
%   simplest sound bound (interval arithmetic); the tighter backward-CROWN pass
%   layers on top of this same op list / boundary.
%
%   Inputs:
%     ops       - cell array of op structs (see nn_to_ops):
%                   struct('type','affine','W',W,'b',b)   % y = W*x + b
%                   struct('type','relu')                 % y = max(0,x)
%     in_lb,in_ub - input box, n-by-1 (or n-by-BATCH for batched sub-domains)
%     precision - 'single' (default; fast, T4/A10G-native) | 'double'
%                 (tighter FP, slower; weak on T4/A10G, full-rate on V100).
%                 R1 knob: one cast at the boundary; everything downstream
%                 (the matmuls + the max) inherits the class, gpuArray included.
%
%   Output: [out_lb, out_ub] -- for EVERY x in [in_lb,in_ub], net(x) in [out_lb,out_ub].
%
%   SOUNDNESS: interval arithmetic is a guaranteed over-approximation. Affine:
%   split W into W+ = max(W,0), W- = min(W,0); the extremal outputs are
%   W+*lb + W-*ub + b (lower) and W+*ub + W-*lb + b (upper). ReLU is monotone, so
%   bounds map elementwise through max(0,.). No input is dropped.
%
%   NOTE: caller controls device. To run on GPU, pass in_lb/in_ub as gpuArray and
%   ensure ops' weights are gpuArray (or this function will gather to host via the
%   first matmul). For sound CERTIFIED runs, pair 'double' with outward rounding
%   (round lb down / ub up per op) -- a follow-on; this version is dev-precision.

    if nargin < 4 || isempty(precision)
        precision = 'single';
    end
    if ~ismember(precision, {'single', 'double'})
        error('gpu_bab_ibp:precision', "precision must be 'single' or 'double'");
    end

    lb = cast(in_lb, precision);
    ub = cast(in_ub, precision);
    nOps = numel(ops);
    hasAdd = any(cellfun(@(o) strcmp(o.type, 'add'), ops));

    if ~hasAdd
        % sequential (chain) net -- rolling bounds, no per-op cache.
        for k = 1:nOps
            [lb, ub] = i_apply_ibp(ops{k}, lb, ub, precision);
        end
        out_lb = lb; out_ub = ub;
        return;
    end

    % DAG net (residual): cache each op's output bounds (index k+1 = op k; index 1 = op 0 =
    % input), so an 'add' can sum two non-adjacent ops. Interval add is EXACT (sound + tight).
    cl = cell(nOps + 1, 1); cu = cell(nOps + 1, 1);
    cl{1} = lb; cu{1} = ub;
    for k = 1:nOps
        op = ops{k};
        if strcmp(op.type, 'add')
            a = op.inputs(1) + 1; b = op.inputs(2) + 1;
            cl{k+1} = cl{a} + cl{b};
            cu{k+1} = cu{a} + cu{b};
        else
            [cl{k+1}, cu{k+1}] = i_apply_ibp(op, cl{k}, cu{k}, precision);
        end
    end
    out_lb = cl{nOps + 1};
    out_ub = cu{nOps + 1};
end

function [nlb, nub] = i_apply_ibp(op, lb, ub, precision)
% Sound IBP for a single (single-input) op. 'add' is multi-input and handled by the caller.
    switch op.type
        case 'affine'
            W = cast(op.W, precision); b = cast(op.b(:), precision);
            Wp = max(W, 0); Wn = min(W, 0);
            nlb = Wp * lb + Wn * ub + b;
            nub = Wp * ub + Wn * lb + b;
        case 'relu'
            nlb = max(lb, 0); nub = max(ub, 0);
        case 'conv'
            [nlb, nub] = i_conv_ibp(op, lb, ub, precision);
        case 'normaffine'
            [nlb, nub] = i_normaffine_ibp(op, lb, ub, precision);
        case 'avgpool'
            [nlb, nub] = i_avgpool_ibp(op, lb, ub, precision);
        case 'maxpool'
            [nlb, nub] = i_maxpool_ibp(op, lb, ub, precision);
        otherwise
            error('gpu_bab_ibp:unsupportedOp', ...
                'Unsupported op type "%s" (affine/relu/conv/normaffine/avgpool/maxpool/add).', op.type);
    end
end

function [olb, oub] = i_conv_ibp(op, lb, ub, precision)
% Conv IBP = interval arithmetic over the LINEAR conv map: the SAME W+/W- split as the
% affine op, with the matmul replaced by dlconv. Sound (tightest interval for a linear
% map). Flat columns <-> [H W C B] (column-major) at the boundary; spec/batch = dim 4.
    B = size(lb, 2);
    ish = op.inShape; osh = op.outShape;
    W = cast(op.W, precision);
    bb = cast(op.b(:), precision);
    Wp = max(W, 0); Wn = min(W, 0);
    L4 = dlarray(reshape(lb, [ish(1) ish(2) ish(3) B]), 'SSCB');
    U4 = dlarray(reshape(ub, [ish(1) ish(2) ish(3) B]), 'SSCB');
    pad2 = [op.pad(1) op.pad(3); op.pad(2) op.pad(4)];   % [t l; b r] for dlconv
    args = {'Stride', op.stride, 'Padding', pad2, 'DilationFactor', op.dil};
    Lo = dlconv(L4, Wp, bb, args{:}) + dlconv(U4, Wn, 0, args{:});   % lower
    Hi = dlconv(U4, Wp, bb, args{:}) + dlconv(L4, Wn, 0, args{:});   % upper
    olb = reshape(extractdata(Lo), [prod(osh) B]);
    oub = reshape(extractdata(Hi), [prod(osh) B]);
end

function [olb, oub] = i_normaffine_ibp(op, lb, ub, precision)
% Per-element affine y = s.*x + t (s,t broadcast over [H W C]); sound interval map.
    s = cast(op.scale, precision); t = cast(op.shift, precision);
    sh = op.shape; B = size(lb, 2);
    L4 = reshape(lb, [sh(1) sh(2) sh(3) B]);
    U4 = reshape(ub, [sh(1) sh(2) sh(3) B]);
    a = s .* L4 + t; c = s .* U4 + t;
    olb = reshape(min(a, c), [prod(sh) B]);
    oub = reshape(max(a, c), [prod(sh) B]);
end

function [olb, oub] = i_avgpool_ibp(op, lb, ub, precision)
% Average pool = per-channel mean over each window: all +weights -> MONOTONE, so the
% bounds map directly (avgpool(lb), avgpool(ub)) with no W+/W- split. Non-overlapping,
% unpadded (enforced at extraction). dlarray avgpool with Stride=pool.
    B = size(lb, 2); ish = op.inShape; osh = op.outShape;
    L4 = dlarray(reshape(lb, [ish(1) ish(2) ish(3) B]), 'SSCB');
    U4 = dlarray(reshape(ub, [ish(1) ish(2) ish(3) B]), 'SSCB');
    Lo = avgpool(L4, op.pool, 'Stride', op.stride);
    Hi = avgpool(U4, op.pool, 'Stride', op.stride);
    olb = reshape(extractdata(Lo), [prod(osh) B]);
    oub = reshape(extractdata(Hi), [prod(osh) B]);
end

function [olb, oub] = i_maxpool_ibp(op, lb, ub, precision)
% Max pool = per-channel max over each window: MONOTONE, so bounds map directly
% (maxpool(lb), maxpool(ub)) -- exact interval, no relaxation. dlarray maxpool.
    B = size(lb, 2); ish = op.inShape; osh = op.outShape;
    L4 = dlarray(reshape(cast(lb, precision), [ish(1) ish(2) ish(3) B]), 'SSCB');
    U4 = dlarray(reshape(cast(ub, precision), [ish(1) ish(2) ish(3) B]), 'SSCB');
    Lo = maxpool(L4, op.pool, 'Stride', op.stride);
    Hi = maxpool(U4, op.pool, 'Stride', op.stride);
    olb = reshape(extractdata(Lo), [prod(osh) B]);
    oub = reshape(extractdata(Hi), [prod(osh) B]);
end
