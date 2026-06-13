function [out_lb, out_ub] = gpu_bab_crown(ops, x_lb, x_ub, precision)
% GPU_BAB_CROWN  Sound CROWN-IBP output bounds for a feedforward ReLU net.
%
%   [out_lb, out_ub] = GPU_BAB_CROWN(ops, x_lb, x_ub, precision) computes SOUND
%   output bounds over the input box [x_lb, x_ub] using backward linear
%   relaxation (CROWN) with IBP-derived pre-activation bounds (CROWN-IBP). This
%   is the tight, LP-free core of NNV's GPU-BaB engine: a forward IBP pass for
%   the per-neuron ReLU relaxation, then a single backward pass of linear
%   coefficients to the input, concretised by the box's dual norm. Tighter than
%   gpu_bab_ibp by design, and still pure linear algebra -> gpuArray-ready.
%   Configurable precision (R1): 'single' (default) | 'double'.
%
%   This version is SINGLE-BOX (x_lb,x_ub are n-by-1) to validate the relaxation
%   math; the batched form (coefficients n_out-by-n-by-BATCH via pagemtimes, the
%   GPU-BaB sub-domain dimension) layers on top with identical relaxations.
%
%   SOUNDNESS: for every x in [x_lb,x_ub], net(x) in [out_lb,out_ub].
%     - affine z = W h + b: linear, propagates exactly (A := A*W, d += A*b).
%     - ReLU h = max(0,z) over pre-activation [l,u]:
%         u <= 0 (stable off): h = 0.
%         l >= 0 (stable on) : h = z.
%         l < 0 < u (unstable): upper line h <= au*z + bu, au = u/(u-l),
%           bu = -au*l >= 0; lower line h >= al*z, al in {0,1} adaptive
%           (al = 1 iff u >= -l, the min-area CROWN slope).
%       Sign-aware selection: for the UPPER output bound a positive coefficient
%       takes the upper line, a negative coefficient the lower line (and the
%       mirror for the LOWER output bound) -- this is what keeps it sound.
%
%   See gpu_bab_ibp (the looser foundation) and nn_to_ops (the op list).

    if nargin < 4 || isempty(precision)
        precision = 'single';
    end
    if ~ismember(precision, {'single', 'double'})
        error('gpu_bab_crown:precision', "precision must be 'single' or 'double'");
    end

    nOps = numel(ops);

    % ---- Forward IBP pass: capture pre-activation bounds at each ReLU ----
    preL = cell(nOps, 1);
    preU = cell(nOps, 1);
    lb = cast(x_lb, precision);
    ub = cast(x_ub, precision);
    for k = 1:nOps
        op = ops{k};
        switch op.type
            case 'affine'
                W = cast(op.W, precision); b = cast(op.b(:), precision);
                Wp = max(W, 0); Wn = min(W, 0);
                nlb = Wp*lb + Wn*ub + b;
                nub = Wp*ub + Wn*lb + b;
                lb = nlb; ub = nub;
            case 'relu'
                preL{k} = lb;   % pre-activation lower bound feeding this ReLU
                preU{k} = ub;   % pre-activation upper bound
                lb = max(lb, 0); ub = max(ub, 0);
            otherwise
                error('gpu_bab_crown:unsupportedOp', ...
                    'Unsupported op type "%s".', op.type);
        end
    end

    nOut = numel(lb);
    % ---- Backward CROWN passes: upper and lower output bounds ----
    out_ub = crown_backward(ops, preL, preU, cast(x_lb, precision), ...
                            cast(x_ub, precision), nOut, precision, true);
    out_lb = crown_backward(ops, preL, preU, cast(x_lb, precision), ...
                            cast(x_ub, precision), nOut, precision, false);
end

function out = crown_backward(ops, preL, preU, x_lb, x_ub, nOut, precision, upper)
% Backward pass of linear coefficients A (nOut-by-current_dim) + bias d to the
% input, then concretise. upper=true -> output upper bound; false -> lower bound.
    A = eye(nOut, precision);
    d = zeros(nOut, 1, precision);
    for k = numel(ops):-1:1
        op = ops{k};
        switch op.type
            case 'affine'
                W = cast(op.W, precision); b = cast(op.b(:), precision);
                d = d + A * b;
                A = A * W;
            case 'relu'
                l = preL{k}; u = preU{k};
                m = numel(l);
                au = zeros(m, 1, precision);   % upper-line slope
                bu = zeros(m, 1, precision);   % upper-line intercept (>= 0)
                al = zeros(m, 1, precision);   % lower-line slope (no intercept)
                act = (l >= 0);                % stable-on -> identity
                au(act) = 1; al(act) = 1;
                unst = (l < 0) & (u > 0);      % unstable -> relaxed
                denom = u(unst) - l(unst);
                au(unst) = u(unst) ./ denom;
                bu(unst) = -au(unst) .* l(unst);            % = -u*l/(u-l) >= 0
                al(unst) = cast(u(unst) >= -l(unst), precision); % min-area slope
                % stable-off (u<=0) keeps au=al=bu=0
                Apos = max(A, 0); Aneg = min(A, 0);
                if upper
                    % positive coeff -> upper line; negative coeff -> lower line
                    d = d + Apos * bu;
                    A = Apos .* au.' + Aneg .* al.';
                else
                    % mirror: positive coeff -> lower line; negative -> upper line
                    d = d + Aneg * bu;
                    A = Apos .* al.' + Aneg .* au.';
                end
            otherwise
                error('gpu_bab_crown:unsupportedOp', ...
                    'Unsupported op type "%s".', op.type);
        end
    end
    Apos = max(A, 0); Aneg = min(A, 0);
    if upper
        out = Apos * x_ub + Aneg * x_lb + d;   % max over the box
    else
        out = Apos * x_lb + Aneg * x_ub + d;   % min over the box
    end
end
