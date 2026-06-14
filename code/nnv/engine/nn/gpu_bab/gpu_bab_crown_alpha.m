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
%   optimized -- the LP-free tightening that makes CROWN competitive on MNIST. Gradient
%   via MATLAB autodiff (dlgradient through pagemtimes); pure MATLAB, gpuArray/dlarray-
%   compatible, configurable precision, batched over B sub-domains.
%
%   SOUNDNESS: any alpha in [0,1] is a valid lower ReLU relaxation, so EVERY iterate is a
%   sound lower bound (ascent only tightens it). The upper slope au, intercept bu, and the
%   stable neurons are FIXED; only the unstable lower slopes are optimized. All free alphas
%   are concatenated into ONE dlarray vector (dlgradient needs a single variable, not a
%   cell) and sliced per layer inside the traced backward pass.

    if nargin < 5 || isempty(precision), precision = 'single'; end
    if nargin < 6 || isempty(nIter), nIter = 20; end
    if nargin < 7 || isempty(lr), lr = 0.5; end

    B = size(x_lb, 2);
    nSpec = size(C, 1);
    nOps = numel(ops);

    % ---- forward IBP: pre-activation bounds + fixed relaxation pieces + masks ----
    preL = cell(nOps,1); preU = cell(nOps,1); actM = cell(nOps,1); unsM = cell(nOps,1);
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

    % ---- flat alpha vector layout (one variable for dlgradient) ----
    reluIdx = find(cellfun(@(o) strcmp(o.type,'relu'), ops));
    rdims = arrayfun(@(k) size(preL{k},1), reluIdx);
    offsets = [0; cumsum(rdims(:) * B)];
    a0 = zeros(offsets(end), 1, precision);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        a = zeros(size(preL{k}), precision);
        uns = unsM{k} > 0;
        % min-area init: al = 1 if u >= -l, else 0 (only unstable matter)
        a(uns) = cast(preU{k}(uns) >= -preL{k}(uns), precision);
        a0(offsets(r)+1 : offsets(r+1)) = a(:);
    end
    alphaVec = dlarray(a0);

    m0 = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, ...
                      alphaVec, reluIdx, rdims, offsets, B, nSpec);
    info = struct('base_minmargin', i_gather(sum(min(m0,[],1),2)), 'iters', nIter);

    % ---- projected gradient ascent ----
    for it = 1:nIter
        [~, grad] = dlfeval(@(av) i_alpha_loss(ops, preL, actM, unsM, auC, buC, ...
                            x_lb, x_ub, C, precision, av, reluIdx, rdims, offsets, B, nSpec), alphaVec);
        v = extractdata(alphaVec) - lr * extractdata(grad);   % descend loss = ascend margin
        v = max(min(v, 1), 0);
        alphaVec = dlarray(v);
    end

    margins = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, ...
                           alphaVec, reluIdx, rdims, offsets, B, nSpec);
    info.alpha_minmargin = i_gather(sum(min(margins,[],1),2));
end

function [loss, grad] = i_alpha_loss(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphaVec, reluIdx, rdims, offsets, B, nSpec)
    margins = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, ...
                           alphaVec, reluIdx, rdims, offsets, B, nSpec);
    loss = -sum(min(margins, [], 1), 'all');
    grad = dlgradient(loss, alphaVec);
end

function margins = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphaVec, reluIdx, rdims, offsets, B, nSpec)
    A = repmat(cast(C, precision), 1, 1, B);
    d = zeros(nSpec, B, precision);
    for k = numel(ops):-1:1
        op = ops{k};
        if strcmp(op.type, 'affine')
            W = cast(op.W, precision); b = cast(op.b(:), precision);
            d = d + reshape(pagemtimes(A, b), nSpec, B);
            A = pagemtimes(A, W);
        else
            r = find(reluIdx == k, 1);
            dim = rdims(r);
            alpha_k = reshape(alphaVec(offsets(r)+1 : offsets(r+1)), dim, B);  % dim x B (traced)
            au = auC{k}; bu = buC{k};
            al = actM{k} + unsM{k} .* alpha_k;          % active->1, inactive->0, unstable->alpha
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

function y = i_gather(x)
    if isa(x, 'dlarray'), x = extractdata(x); end
    if isa(x, 'gpuArray'), x = gather(x); end
    y = x;
end
