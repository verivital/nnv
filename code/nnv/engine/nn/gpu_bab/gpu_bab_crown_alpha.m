function [margins, info] = gpu_bab_crown_alpha(ops, x_lb, x_ub, C, precision, nIter, lr, useTight)
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
    if nargin < 8 || isempty(useTight), useTight = false; end  % true -> CROWN-tight intermediate bounds

    B = size(x_lb, 2);
    nSpec = size(C, 1);
    nOps = numel(ops);

    % ---- intermediate bounds: IBP (default) or CROWN-tight (useTight; single-box) ----
    preL = cell(nOps,1); preU = cell(nOps,1);
    if useTight
        [~, preL, preU] = gpu_bab_crown_tight(ops, x_lb, x_ub, C, precision);  % assumes B==1
    else
        lb = cast(x_lb, precision); ub = cast(x_ub, precision);
        for k = 1:nOps
            op = ops{k};
            if strcmp(op.type, 'affine')
                W = cast(op.W, precision); b = cast(op.b(:), precision);
                Wp = max(W,0); Wn = min(W,0);
                nlb = Wp*lb + Wn*ub + b; nub = Wp*ub + Wn*lb + b; lb = nlb; ub = nub;
            else
                preL{k} = lb; preU{k} = ub;
                lb = max(lb,0); ub = max(ub,0);
            end
        end
    end
    % ---- derive fixed relaxation pieces + masks from preL/preU (tight or IBP) ----
    actM = cell(nOps,1); unsM = cell(nOps,1); auC = cell(nOps,1); buC = cell(nOps,1);
    for k = 1:nOps
        if strcmp(ops{k}.type, 'relu')
            l = preL{k}; u = preU{k};
            act = (l >= 0); uns = (l < 0) & (u > 0);
            au = zeros(size(l), precision); bu = zeros(size(l), precision);
            au(act) = 1;
            dn = u(uns) - l(uns);
            au(uns) = u(uns) ./ dn; bu(uns) = -au(uns) .* l(uns);
            actM{k} = cast(act, precision); unsM{k} = cast(uns, precision);
            auC{k} = au; buC{k} = bu;
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
    bestObj = i_gather(sum(min(m0,[],1), 'all'));
    bestVec = alphaVec;
    info = struct('base_minmargin', bestObj, 'iters', nIter);

    % ---- projected gradient ascent (normalized step + keep-best -> monotone) ----
    for it = 1:nIter
        [~, grad] = dlfeval(@(av) i_alpha_loss(ops, preL, actM, unsM, auC, buC, ...
                            x_lb, x_ub, C, precision, av, reluIdx, rdims, offsets, B, nSpec), alphaVec);
        g = extractdata(grad);
        g = g / (max(abs(g(:))) + eps(precision));        % normalize -> step ~ lr per coord
        v = extractdata(alphaVec) - lr * g;               % descend loss = ascend margin
        v = max(min(v, 1), 0);
        alphaVec = dlarray(v);
        m = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, ...
                         alphaVec, reluIdx, rdims, offsets, B, nSpec);
        obj = i_gather(sum(min(m,[],1), 'all'));
        if obj > bestObj, bestObj = obj; bestVec = alphaVec; end   % keep best (never regress)
    end

    margins = i_alpha_back(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, ...
                           bestVec, reluIdx, rdims, offsets, B, nSpec);
    info.alpha_minmargin = i_gather(sum(min(margins,[],1), 'all'));
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
