function [margins, unstable, preL, preU] = gpu_bab_crown_alpha_beta(ops, x_lb, x_ub, C, fixings, reluIdx, precision, nIter, lr, rootBounds)
% GPU_BAB_CROWN_ALPHA_BETA  Joint alpha-CROWN + beta-CROWN bound for a ReLU-split BaB node.
%
%   The per-node bound the batched ReLU-split BaB uses, combining THREE LP-free tightenings:
%     1. ROOT-TIGHT intermediate bounds (rootBounds, reused + clamped per node) for the upper
%        ReLU line (au,bu) -- tight pre-activation bounds at batched speed;
%     2. ALPHA: optimize the unstable-ReLU lower slopes (alpha in [0,1]) to maximize the bound;
%     3. BETA: Lagrangian duals (beta >= 0) for the node's SPLIT CONSTRAINTS s_i*z_i >= 0
%        (s_i=+1 active fix / -1 inactive fix). This is the per-node constraint PROPAGATION the
%        plain clamp lacks: it couples the split neuron's pre-activation (hence earlier layers /
%        the input) into the bound, the piece that closes the gap on hard BaB nodes.
%   alpha + beta are optimized JOINTLY by projected gradient ascent (keep-best). Batched over the
%   B node columns. FC+ReLU (affine/relu) only.
%
%   SOUNDNESS: the bound is max_{alpha in [0,1], beta>=0} min_x [C*f(x) - sum_i beta_i s_i z_i(x)]
%   over the (clamped) input box. For ANY alpha in [0,1] and beta>=0 this is a valid lower bound
%   on min_x C*f(x) s.t. the split constraints (weak Lagrangian duality + a sound ReLU
%   relaxation), so EVERY iterate -- and the returned best -- is sound. beta init = 0 (so the
%   start == the alpha/root-tight bound; optimization only raises it). nIter<=0 -> fixed slope,
%   beta=0 (== plain clamped/root-tight CROWN).

    if nargin < 7 || isempty(precision), precision = 'single'; end
    if nargin < 8 || isempty(nIter), nIter = 20; end
    if nargin < 9 || isempty(lr), lr = 0.2; end
    if nargin < 10, rootBounds = []; end

    B = size(x_lb, 2); nSpec = size(C, 1); nOps = numel(ops);

    % ---- per-relu (clamped) pre-activation bounds + fixed upper line (au,bu) + masks. Same as
    %      gpu_bab_crown_alpha_fix: root-tight reuse when rootBounds given, else IBP forward. ----
    preL = cell(nOps,1); preU = cell(nOps,1);
    actM = cell(nOps,1); unsM = cell(nOps,1);
    auC = cell(nOps,1);  buC = cell(nOps,1);
    unstable = cell(nOps,1); fixSign = cell(nOps,1);
    if isempty(rootBounds)
        lb = cast(x_lb, precision); ub = cast(x_ub, precision);
        for k = 1:nOps
            op = ops{k};
            if strcmp(op.type, 'affine')
                W = cast(op.W, precision); b = cast(op.b(:), precision);
                Wp = max(W,0); Wn = min(W,0);
                nlb = Wp*lb + Wn*ub + b; nub = Wp*ub + Wn*lb + b; lb = nlb; ub = nub;
            else
                fx = i_fx(fixings, k, size(lb));
                lb(fx == 1)  = max(lb(fx == 1),  0);
                ub(fx == -1) = min(ub(fx == -1), 0);
                [preL{k}, preU{k}, actM{k}, unsM{k}, auC{k}, buC{k}, unstable{k}] = i_relu_pre(lb, ub, precision);
                fixSign{k} = cast(fx, precision);
                lb = max(lb,0); ub = max(ub,0);
            end
        end
    else
        tmpl = cast(x_lb(1), precision);
        for k = 1:nOps
            if strcmp(ops{k}.type, 'relu')
                l = repmat(cast(rootBounds.preL{k}, 'like', tmpl), 1, B);
                u = repmat(cast(rootBounds.preU{k}, 'like', tmpl), 1, B);
                fx = i_fx(fixings, k, size(l));
                l(fx == 1)  = max(l(fx == 1),  0);
                u(fx == -1) = min(u(fx == -1), 0);
                [preL{k}, preU{k}, actM{k}, unsM{k}, auC{k}, buC{k}, unstable{k}] = i_relu_pre(l, u, precision);
                fixSign{k} = cast(fx, precision);
            end
        end
    end

    % ---- flat parameter vector p = [alpha; beta]; alpha init = min-area, beta init = 0 ----
    rdims = arrayfun(@(k) size(preL{k},1), reluIdx);
    offsets = [0; cumsum(rdims(:) * B)];
    nP = offsets(end);
    a0 = zeros(nP, 1, precision);
    fixedMask = zeros(nP, 1, precision);          % beta is free only on FIXED neurons
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        a = zeros(size(preL{k}), precision);
        uns = unsM{k} > 0;
        a(uns) = cast(preU{k}(uns) >= -preL{k}(uns), precision);
        a0(offsets(r)+1 : offsets(r+1)) = a(:);
        fm = (fixSign{k} ~= 0);
        fixedMask(offsets(r)+1 : offsets(r+1)) = cast(fm(:), precision);
    end
    pVec = dlarray([a0; zeros(nP,1,precision)]);  % [alpha; beta], beta=0

    if nIter <= 0
        margins = i_g(i_back_ab(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, pVec, reluIdx, rdims, offsets, nP, B, nSpec));
        return;
    end

    % ---- projected gradient ascent over [alpha; beta] (normalized step, keep-best) ----
    m0 = i_back_ab(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, pVec, reluIdx, rdims, offsets, nP, B, nSpec);
    bestObj = i_g(sum(min(m0,[],1), 'all')); bestVec = pVec;
    for it = 1:nIter
        [~, grad] = dlfeval(@(p) i_loss_ab(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, p, reluIdx, rdims, offsets, nP, B, nSpec), pVec);
        g = extractdata(grad); g = g / (max(abs(g(:))) + eps(precision));
        v = extractdata(pVec) - lr * g;
        v(1:nP)      = max(min(v(1:nP), 1), 0);   % alpha in [0,1]
        v(nP+1:end)  = max(v(nP+1:end), 0) .* fixedMask;   % beta >= 0, only on fixed neurons
        pVec = dlarray(v);
        m = i_back_ab(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, pVec, reluIdx, rdims, offsets, nP, B, nSpec);
        obj = i_g(sum(min(m,[],1), 'all'));
        if obj > bestObj, bestObj = obj; bestVec = pVec; end
    end
    margins = i_g(i_back_ab(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, bestVec, reluIdx, rdims, offsets, nP, B, nSpec));
end

function [loss, grad] = i_loss_ab(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, pVec, reluIdx, rdims, offsets, nP, B, nSpec)
    margins = i_back_ab(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, pVec, reluIdx, rdims, offsets, nP, B, nSpec);
    loss = -sum(min(margins, [], 1), 'all');
    grad = dlgradient(loss, pVec);
end

function margins = i_back_ab(ops, preL, actM, unsM, auC, buC, fixSign, x_lb, x_ub, C, precision, pVec, reluIdx, rdims, offsets, nP, B, nSpec)
% Backward CROWN of [C*f(x) - sum_i beta_i s_i z_i(x)]: at each ReLU the standard sign-aware
% relaxation (lower slope al = act + uns*alpha), THEN the beta dual adds -(beta_k.*s_k) to the
% coefficient on that ReLU's pre-activation z_k (broadcast over the nSpec spec rows), which then
% propagates back through the preceding affine -- the split-constraint term.
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
            alpha_k = reshape(pVec(offsets(r)+1 : offsets(r+1)), dim, B);
            beta_k  = reshape(pVec(nP + offsets(r)+1 : nP + offsets(r+1)), dim, B);
            au = auC{k}; bu = buC{k};
            al = actM{k} + unsM{k} .* alpha_k;
            Apos = max(A, 0); Aneg = min(A, 0);
            d = d + reshape(pagemtimes(Aneg, reshape(bu, dim, 1, B)), nSpec, B);
            A = Apos .* reshape(al, 1, dim, B) + Aneg .* reshape(au, 1, dim, B);
            % beta split-dual: subtract beta_k*s_k from the coefficient on z_k (all spec rows)
            A = A - reshape(beta_k .* fixSign{k}, 1, dim, B);
        end
    end
    n = size(x_lb, 1);
    Apos = max(A, 0); Aneg = min(A, 0);
    margins = reshape(pagemtimes(Apos, reshape(cast(x_lb, precision), n, 1, B)), nSpec, B) ...
            + reshape(pagemtimes(Aneg, reshape(cast(x_ub, precision), n, 1, B)), nSpec, B) + d;
end

function fx = i_fx(fixings, k, sz)
% Per-relu fixing vector (dim x B), 0 if none provided.
    if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k}), fx = fixings{k};
    else, fx = zeros(sz); end
end

function y = i_g(x)
    if isa(x, 'dlarray'), x = extractdata(x); end
    if isa(x, 'gpuArray'), x = gather(x); end
    y = x;
end

function [pl, pu, actM, unsM, au, bu, uns] = i_relu_pre(l, u, precision)
    pl = l; pu = u;
    act = (l >= 0); uns = (l < 0) & (u > 0);
    au = zeros(size(l), precision); bu = zeros(size(l), precision);
    au(act) = 1;
    dn = u(uns) - l(uns);
    au(uns) = u(uns) ./ dn; bu(uns) = -au(uns) .* l(uns);
    actM = cast(act, precision); unsM = cast(uns, precision);
end
