function [margins, unstable, preL, preU] = gpu_bab_crown_alpha_fix(ops, x_lb, x_ub, C, fixings, reluIdx, precision, nIter, lr, rootBounds)
% GPU_BAB_CROWN_ALPHA_FIX  alpha-CROWN bound for a ReLU-split BaB node.
%
%   [margins, unstable] = GPU_BAB_CROWN_ALPHA_FIX(ops, x_lb, x_ub, C, fixings, reluIdx,
%                                                 precision, nIter, lr)
%   Forward IBP with the node's neuron fixings CLAMPED into the pre-activation bounds
%   (active fix z>=0 -> l:=max(l,0); inactive z<=0 -> u:=min(u,0)), then optimize the
%   unstable lower slopes (alpha in [0,1]) to MAXIMIZE the spec lower bound. Returns the
%   tightened margins (nSpec-by-1) and the per-relu unstable mask (for split selection).
%   nIter<=0 -> fixed min-area slope (no optimization), i.e. plain clamped CROWN.
%
%   This is the bound the ReLU-split BaB uses per node: tighter node bounds -> far fewer
%   splits needed to certify. SOUND: any alpha in [0,1] is a valid lower ReLU relaxation,
%   and clamping l up / u down is a valid (tighter) bound on the node's sub-domain, so
%   every iterate is a sound lower bound. Single-box (B assumed 1 for a node).

    if nargin < 7 || isempty(precision), precision = 'single'; end
    if nargin < 8 || isempty(nIter), nIter = 20; end
    if nargin < 9 || isempty(lr), lr = 0.2; end
    if nargin < 10, rootBounds = []; end

    B = size(x_lb, 2);
    nSpec = size(C, 1);
    nOps = numel(ops);

    % ---- per-relu pre-activation bounds (clamped by the node fixings) + the fixed upper-line
    % relaxation (au,bu). With rootBounds, REUSE the tight root intermediate bounds (clamped per
    % node) instead of the loose IBP forward -- tight bounds AND optimized lower slopes (the
    % strongest LP-free node bound). SOUND either way: clamping l up / u down is a valid tighter
    % bound on the node sub-domain, and the reused root bounds hold over the full box >= it. ----
    preL = cell(nOps,1); preU = cell(nOps,1);
    actM = cell(nOps,1); unsM = cell(nOps,1);
    auC = cell(nOps,1);  buC = cell(nOps,1);
    unstable = cell(nOps,1);
    if isempty(rootBounds)
        lb = cast(x_lb, precision); ub = cast(x_ub, precision);
        for k = 1:nOps
            op = ops{k};
            if strcmp(op.type, 'affine')
                W = cast(op.W, precision); b = cast(op.b(:), precision);
                Wp = max(W,0); Wn = min(W,0);
                nlb = Wp*lb + Wn*ub + b; nub = Wp*ub + Wn*lb + b; lb = nlb; ub = nub;
            else
                fx = fixings{k};
                if ~isempty(fx)
                    lb(fx == 1)  = max(lb(fx == 1),  0);
                    ub(fx == -1) = min(ub(fx == -1), 0);
                end
                [preL{k}, preU{k}, actM{k}, unsM{k}, auC{k}, buC{k}, unstable{k}] = i_relu_pre(lb, ub, precision);
                lb = max(lb,0); ub = max(ub,0);
            end
        end
    else
        tmpl = cast(x_lb(1), precision);   % scalar device/precision template (no n-by-B copy)
        for k = 1:nOps
            if strcmp(ops{k}.type, 'relu')
                l = repmat(cast(rootBounds.preL{k}, 'like', tmpl), 1, B);
                u = repmat(cast(rootBounds.preU{k}, 'like', tmpl), 1, B);
                fx = fixings{k};
                if ~isempty(fx)
                    l(fx == 1)  = max(l(fx == 1),  0);
                    u(fx == -1) = min(u(fx == -1), 0);
                end
                [preL{k}, preU{k}, actM{k}, unsM{k}, auC{k}, buC{k}, unstable{k}] = i_relu_pre(l, u, precision);
            end
        end
    end

    % ---- flat alpha vector (min-area init) ----
    rdims = arrayfun(@(k) size(preL{k},1), reluIdx);
    offsets = [0; cumsum(rdims(:) * B)];
    a0 = zeros(offsets(end), 1, precision);
    for r = 1:numel(reluIdx)
        k = reluIdx(r);
        a = zeros(size(preL{k}), precision);
        uns = unsM{k} > 0;
        a(uns) = cast(preU{k}(uns) >= -preL{k}(uns), precision);
        a0(offsets(r)+1 : offsets(r+1)) = a(:);
    end
    alphaVec = dlarray(a0);

    if nIter <= 0
        margins = i_aback(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, ...
                          alphaVec, reluIdx, rdims, offsets, B, nSpec);
        margins = i_g(margins);    % extract from dlarray -> numeric (callers index/find on it)
        return;
    end

    % ---- projected gradient ascent (normalized step + keep-best) ----
    m0 = i_aback(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, ...
                 alphaVec, reluIdx, rdims, offsets, B, nSpec);
    bestObj = i_g(sum(min(m0,[],1), 'all')); bestVec = alphaVec;
    for it = 1:nIter
        [~, grad] = dlfeval(@(av) i_aloss(ops, preL, actM, unsM, auC, buC, ...
                            x_lb, x_ub, C, precision, av, reluIdx, rdims, offsets, B, nSpec), alphaVec);
        g = extractdata(grad); g = g / (max(abs(g(:))) + eps(precision));
        v = extractdata(alphaVec) - lr * g; v = max(min(v,1),0);
        alphaVec = dlarray(v);
        m = i_aback(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, ...
                    alphaVec, reluIdx, rdims, offsets, B, nSpec);
        obj = i_g(sum(min(m,[],1), 'all'));
        if obj > bestObj, bestObj = obj; bestVec = alphaVec; end
    end
    margins = i_aback(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, ...
                      bestVec, reluIdx, rdims, offsets, B, nSpec);
    margins = i_g(margins);        % extract from dlarray -> numeric (callers index/find on it)
end

function [loss, grad] = i_aloss(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphaVec, reluIdx, rdims, offsets, B, nSpec)
    margins = i_aback(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphaVec, reluIdx, rdims, offsets, B, nSpec);
    loss = -sum(min(margins, [], 1), 'all');
    grad = dlgradient(loss, alphaVec);
end

function margins = i_aback(ops, preL, actM, unsM, auC, buC, x_lb, x_ub, C, precision, alphaVec, reluIdx, rdims, offsets, B, nSpec)
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
            alpha_k = reshape(alphaVec(offsets(r)+1 : offsets(r+1)), dim, B);
            au = auC{k}; bu = buC{k};
            al = actM{k} + unsM{k} .* alpha_k;
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

function y = i_g(x)
    if isa(x, 'dlarray'), x = extractdata(x); end
    if isa(x, 'gpuArray'), x = gather(x); end
    y = x;
end

function [pl, pu, actM, unsM, au, bu, uns] = i_relu_pre(l, u, precision)
% Per-relu (clamped) pre-activation bounds + the FIXED upper-line relaxation au*z+bu and the
% active/unstable masks. The lower slope is left to the alpha optimizer (al = act + uns*alpha).
    pl = l; pu = u;
    act = (l >= 0); uns = (l < 0) & (u > 0);
    au = zeros(size(l), precision); bu = zeros(size(l), precision);
    au(act) = 1;
    dn = u(uns) - l(uns);
    au(uns) = u(uns) ./ dn; bu(uns) = -au(uns) .* l(uns);
    actM = cast(act, precision); unsM = cast(uns, precision);
end
