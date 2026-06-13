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

    for k = 1:numel(ops)
        op = ops{k};
        switch op.type
            case 'affine'
                W = cast(op.W, precision);
                b = cast(op.b(:), precision);
                Wp = max(W, 0);
                Wn = min(W, 0);
                new_lb = Wp * lb + Wn * ub + b;
                new_ub = Wp * ub + Wn * lb + b;
                lb = new_lb;
                ub = new_ub;
            case 'relu'
                lb = max(lb, 0);
                ub = max(ub, 0);
            otherwise
                error('gpu_bab_ibp:unsupportedOp', ...
                    'Unsupported op type "%s" (Phase 1 supports affine + relu).', op.type);
        end
    end

    out_lb = lb;
    out_ub = ub;
end
