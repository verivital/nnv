function [margins, info] = gpu_bab_crown_alpha_joint(ops, x_lb, x_ub, C, precision, nIter, lr, fixings)
% GPU_BAB_CROWN_ALPHA_JOINT  JOINT alpha-CROWN: optimize the unstable-ReLU lower
%   slopes (the "alpha" parameters) INSIDE EVERY intermediate backward pass, not
%   just the final spec. This is the piece gpu_bab_crown_alpha was missing.
%
%   [margins, info] = GPU_BAB_CROWN_ALPHA_JOINT(ops, x_lb, x_ub, C, precision, nIter, lr)
%   maximizes a sound lower bound on C*f(x) over [x_lb,x_ub] by gradient ascent over
%   ALL free lower slopes alpha in [0,1] -- one per unstable neuron per ReLU layer.
%   The SAME alpha of an early layer feeds both (a) that layer's relaxation when the
%   pre-activation bounds preL/preU of LATER layers are computed (the O(L^2) backward
%   recursion of gpu_bab_crown_tight), and (b) the final-spec backward pass. So
%   tightening an early alpha tightens every downstream intermediate bound AND the
%   spec at once -- the joint optimization that makes CROWN competitive (auto_LiRPA
%   alpha-CROWN). gpu_bab_crown_alpha only optimized (b).
%
%   Difference from gpu_bab_crown_alpha:
%     - crown_alpha: intermediate preL/preU are computed ONCE (IBP or fixed-slope
%       crown_tight) and FROZEN; only the final spec's alpha is optimized.
%     - crown_alpha_joint: preL/preU are RECOMPUTED from the current alpha on every
%       forward, so the relaxation pieces (au,bu) themselves move with alpha.
%
%   Optional `fixings` (6th arg, cell per op of -1/0/+1 ReLU-split node fixings) clamps
%   fixed neurons (active l>=0 / inactive u<=0) inside the intermediate bounds, so the
%   ReLU-split BaB can call this per sub-domain (matches gpu_bab_crown_tight's signature).
%
%   SOUNDNESS (-150 rule): for ANY alpha in [0,1] every relaxation line is a valid
%   sound ReLU over-approximation, so every intermediate bound and the final margin
%   are sound lower bounds at every iterate. The reported margin is RE-EVALUATED with
%   the best alpha (i_joint_margin, untraced) and kept-best, so it is a genuine sound
%   lower bound regardless of gradient quality. Masks (which neurons are stable vs
%   unstable) are taken from the DETACHED bounds (constant gates, like ReLU's own
%   subgradient); the au/bu/al VALUES are traced through alpha. Single-box (B==1):
%   the intermediate-bound recursion is per input box (batching is a later step).

    if nargin < 5 || isempty(precision), precision = 'single'; end
    if nargin < 6 || isempty(nIter), nIter = 20; end
    if nargin < 7 || isempty(lr), lr = 0.5; end
    if nargin < 8, fixings = {}; end

    nOps = numel(ops);
    reluIdx = find(cellfun(@(o) strcmp(o.type,'relu'), ops));

    % ---- per-relu width + flat alpha layout (one dlarray for dlgradient) ----
    rdims = zeros(numel(reluIdx),1);
    for r = 1:numel(reluIdx)
        rdims(r) = i_layer_width(ops, reluIdx(r)-1);
    end
    offsets = [0; cumsum(rdims)];

    % ---- min-area init (reproduces gpu_bab_crown_tight exactly at iter 0) ----
    % al_init(unstable) = 1 if u >= -l else 0; computed from the FIXED-slope tight
    % bounds so iter-0 margin == crown_tight margin, then ascent only tightens.
    [~, preL0, preU0] = gpu_bab_crown_tight(ops, x_lb, x_ub, C, precision, fixings);
    a0 = zeros(offsets(end), 1, precision);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        l = preL0{k}; u = preU0{k};
        uns = (l < 0) & (u > 0);
        a = zeros(rdims(r), 1, precision);
        a(uns) = cast(u(uns) >= -l(uns), precision);
        a0(offsets(r)+1 : offsets(r+1)) = a;
    end
    alphaVec = dlarray(a0);

    m0 = i_joint_margin(ops, x_lb, x_ub, C, alphaVec, reluIdx, offsets, rdims, precision, fixings);
    bestObj = i_gather(sum(min(m0,[],1), 'all'));
    bestVec = alphaVec;
    info = struct('base_minmargin', bestObj, 'iters', nIter);

    % ---- projected gradient ascent (normalized step + keep-best -> monotone) ----
    for it = 1:nIter
        [~, grad] = dlfeval(@(av) i_joint_loss(ops, x_lb, x_ub, C, av, ...
                            reluIdx, offsets, rdims, precision, fixings), alphaVec);
        g = extractdata(grad);
        g = g / (max(abs(g(:))) + eps(precision));     % normalize -> step ~ lr per coord
        v = extractdata(alphaVec) - lr * g;            % descend loss = ascend margin
        v = max(min(v, 1), 0);                         % project to [0,1]
        alphaVec = dlarray(v);
        m = i_joint_margin(ops, x_lb, x_ub, C, alphaVec, reluIdx, offsets, rdims, precision, fixings);
        obj = i_gather(sum(min(m,[],1), 'all'));
        if obj > bestObj, bestObj = obj; bestVec = alphaVec; end   % never regress
    end

    margins = i_joint_margin(ops, x_lb, x_ub, C, bestVec, reluIdx, offsets, rdims, precision, fixings);
    info.alpha_minmargin = i_gather(sum(min(margins,[],1), 'all'));
end

% =========================================================================

function [loss, grad] = i_joint_loss(ops, x_lb, x_ub, C, alphaVec, reluIdx, offsets, rdims, precision, fixings)
    margins = i_joint_margin(ops, x_lb, x_ub, C, alphaVec, reluIdx, offsets, rdims, precision, fixings);
    loss = -sum(min(margins, [], 1), 'all');
    grad = dlgradient(loss, alphaVec);
end

function margins = i_joint_margin(ops, x_lb, x_ub, C, alphaVec, reluIdx, offsets, rdims, precision, fixings)
% One forward: (1) tight intermediate bounds preL/preU from the CURRENT alpha, then
% (2) the final spec margin -- both via the same alpha-parameterised backward pass.
    nOps = numel(ops);
    preL = cell(nOps,1); preU = cell(nOps,1);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        nk = rdims(r);
        Ck = eye(nk, precision);
        pu = i_alpha_backward(ops, k-1, Ck, x_lb, x_ub, preL, preU, ...
                              alphaVec, reluIdx, offsets, rdims, precision, false);
        pl = i_alpha_backward(ops, k-1, Ck, x_lb, x_ub, preL, preU, ...
                              alphaVec, reluIdx, offsets, rdims, precision, true);
        if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
            fx = fixings{k};
            pl = pl + max(0 - pl, 0) .* cast(fx(:) == 1, precision);   % active -> l >= 0
            pu = pu - max(pu - 0, 0) .* cast(fx(:) == -1, precision);  % inactive -> u <= 0
        end
        preL{k} = pl; preU{k} = pu;
    end
    margins = i_alpha_backward(ops, nOps, cast(C, precision), x_lb, x_ub, preL, preU, ...
                              alphaVec, reluIdx, offsets, rdims, precision, true);
end

function bound = i_alpha_backward(ops, upto, A0, x_lb, x_ub, preL, preU, alphaVec, reluIdx, offsets, rdims, precision, lower)
% Backward CROWN over ops[1..upto] with initial coeff A0. The unstable lower-line
% slope is alpha (traced from alphaVec); the upper line (au,bu) is derived from the
% (alpha-dependent, traced) preL/preU. Masks come from the DETACHED bounds.
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
            r = find(reluIdx == k, 1);
            alpha_k = alphaVec(offsets(r)+1 : offsets(r+1));     % rdims(r) x 1 (traced)
            l = preL{k}; u = preU{k};                            % traced dlarrays
            ld = i_data(l); ud = i_data(u);                      % detached for masks
            act = cast(ld >= 0, precision);
            uns = cast((ld < 0) & (ud > 0), precision);
            dn = u - l;
            safe = dn + (1 - uns);                               % avoid 0/0 off the unstable set
            auU = u ./ safe;                                     % traced upper slope on unstable
            au = act + uns .* auU;                               % active->1, else->0 (overwritten on uns)
            bu = uns .* (-(auU) .* l);                           % -au*l on unstable, 0 elsewhere
            al = act + uns .* alpha_k(:);                        % active->1, inactive->0, unstable->alpha
            Apos = max(A, 0); Aneg = min(A, 0);
            if lower
                d = d + Aneg * bu;
                A = Apos .* al.' + Aneg .* au.';
            else
                d = d + Apos * bu;
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

function w = i_layer_width(ops, upto)
    w = [];
    for k = upto:-1:1
        if strcmp(ops{k}.type, 'affine'), w = size(ops{k}.W, 1); return; end
    end
    error('gpu_bab_crown_alpha_joint:noaffine', 'no affine op before index %d', upto);
end

function y = i_data(x)
    if isa(x, 'dlarray'), y = extractdata(x); else, y = x; end
end

function y = i_gather(x)
    if isa(x, 'dlarray'), x = extractdata(x); end
    if isa(x, 'gpuArray'), x = gather(x); end
    y = x;
end
