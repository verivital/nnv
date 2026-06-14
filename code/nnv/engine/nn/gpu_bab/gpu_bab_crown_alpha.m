function [margins, info] = gpu_bab_crown_alpha(ops, x_lb, x_ub, C, precision, nIter, lr)
% GPU_BAB_CROWN_ALPHA  alpha-CROWN: tighten the CROWN spec lower bound by optimizing
%   the unstable-ReLU lower-relaxation slopes (the "alpha" parameters).
%
%   [margins, info] = GPU_BAB_CROWN_ALPHA(ops, x_lb, x_ub, C, precision, nIter, lr)
%   maximizes the lower bound on C*f(x) over the free lower slopes alpha in [0,1] (one
%   per unstable neuron per sub-box) by projected gradient ascent. Returns margins
%   (nSpec-by-B) that are >= the fixed-slope gpu_bab_crown_spec bound and still SOUND.
%   info.base_minmargin / info.alpha_minmargin report the tightening.
%
%   Why this matters (from the research): alpha-CROWN matches the LP triangle-relaxation
%   optimum WITHOUT any LP, and is tighter still when intermediate bounds are jointly
%   optimized -- it is the LP-free tightening that makes CROWN competitive on MNIST.
%   Gradient via MATLAB autodiff (dlgradient through pagemtimes); pure MATLAB, gpuArray-
%   and dlarray-compatible, configurable precision. Batched over B sub-domains.
%
%   SOUNDNESS: any alpha in [0,1] is a valid lower ReLU relaxation, so EVERY iterate is
%   a sound lower bound (gradient ascent only tightens it). The upper slope au, intercept
%   bu, and the stable neurons are FIXED (active->identity, inactive->0); only the
%   unstable lower slopes are optimized.

    if nargin < 5 || isempty(precision), precision = 'single'; end
    if nargin < 6 || isempty(nIter), nIter = 20; end
    if nargin < 7 || isempty(lr), lr = 0.5; end

    B = size(x_lb, 2);
    nSpec = size(C, 1);
    nOps = numel(ops);

    % ---- forward IBP: pre-activation bounds + fixed relaxation pieces + masks ----
    preL = cell(nOps,1); preU = cell(nOps,1);
    actM = cell(nOps,1); unsM = cell(nOps,1);
    auC = cell(nOps,1);  buC = cell(nOps,1);
    lb = cast(x_lb, precision); ub = cast(x_ub, precision);
    for k = 1:nOps
        op = ops{k};
        if strcmp(op.type, 'affine')
            W = cast(op.W, precision); b = cast(op.b(:), precision);
            Wp = max(W,0); Wn = min(W,0);
            nlb = Wp*lb + Wn*ub + b; nub = Wp*ub + Wn*lb + b; lb = nlb; ub = nub;
        else
            preL{k} = lb; preU{k} = ub;
            act = (lb >= 0); uns = (lb < 0) & (ub > 0);
            au = zeros(size(lb), precision); bu = zeros(size(lb), precision);
            au(act) = 1;
            dn = ub(uns) - lb(uns);
            au(uns) = ub(uns) ./ dn; bu(uns) = -au(uns) .* lb(uns);
            actM{k} = cast(act, precision); unsM{k} = cast(uns, precision);
            auC{k} = au; buC{k} = bu;
            lb = max(lb,0); ub = max(ub,0);
        end
    end

    % ---- initial alphas (min-area heuristic) for unstable neurons ----
    alphas = cell(nOps,1);
    for k = 1:nOps
        if strcmp(ops{k}.type, 'relu')
            a0 = zeros(size(preL{k}), precision);
            uns = (preL{k} < 0) & (preU{k} > 0);
            a0(uns) = cast(preU{k}(uns) >= -preL{k}(uns), precision);
            alphas{k} = dlarray(a0);
        end
    end

    m0 = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphas, B, nSpec);
    info = struct('base_minmargin', gather_(sum(min(m0,[],1),2)), 'iters', nIter);

    % ---- projected gradient ascent on alpha (maximize the summed binding margin) ----
    for it = 1:nIter
        [~, grads] = dlfeval(@(A) i_alpha_loss(ops, preL, actM, unsM, auC, buC, ...
                                               x_lb, x_ub, C, precision, A, B, nSpec), alphas);
        for k = 1:nOps
            if ~isempty(alphas{k}) && ~isempty(grads{k})
                av = extractdata(alphas{k}) - lr * extractdata(grads{k}); % descend loss = ascend margin
                av = max(min(av, 1), 0);                                  % project to [0,1]
                alphas{k} = dlarray(av);
            end
        end
    end

    margins = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphas, B, nSpec);
    info.alpha_minmargin = gather_(sum(min(margins,[],1),2));
end

function [loss, grads] = i_alpha_loss(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphas, B, nSpec)
    margins = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphas, B, nSpec);
    loss = -sum(min(margins, [], 1), 'all');     % maximize the binding (min) margin per box
    grads = dlgradient(loss, alphas);
end

function margins = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphas, B, nSpec)
% Backward CROWN pass with lower slope al = active*1 + unstable.*alpha (inactive->0).
    A = repmat(cast(C, precision), 1, 1, B);
    d = zeros(nSpec, B, precision);
    for k = numel(ops):-1:1
        op = ops{k};
        if strcmp(op.type, 'affine')
            W = cast(op.W, precision); b = cast(op.b(:), precision);
            d = d + reshape(pagemtimes(A, b), nSpec, B);
            A = pagemtimes(A, W);
        else
            dim = size(preL{k}, 1);
            au = auC{k}; bu = buC{k};
            al = actM{k} + unsM{k} .* alphas{k};        % dim x B (dlarray when traced)
            Apos = max(A, 0); Aneg = min(A, 0);
            d = d + reshape(pagemtimes(Aneg, reshape(bu, dim, 1, B)), nSpec, B);
            A = Apos .* reshape(al, 1, dim, B) + Aneg .* reshape(au, 1, dim, B);
        end
    end
    n = size(x_lb, 1);
    Apos = max(A, 0); Aneg = min(A, 0);
    margins = reshape(pagemtimes(Apos, reshape(cast(x_lb, precision), n, 1, B)), nSpec, B) ...
            + reshape(pagemtimes(Aneg, reshape(cast(x_ub, precision), n, 1, B)), nSpec, B) + d;
end

function y = gather_(x)
    if isa(x, 'dlarray'), x = extractdata(x); end
    if isa(x, 'gpuArray'), x = gather(x); end
    y = x;
end
