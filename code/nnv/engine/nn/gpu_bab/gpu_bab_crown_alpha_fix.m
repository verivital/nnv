function [margins, unstable] = gpu_bab_crown_alpha_fix(ops, x_lb, x_ub, C, fixings, reluIdx, precision, nIter, lr)
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

    B = size(x_lb, 2);
    nSpec = size(C, 1);
    nOps = numel(ops);

    % ---- forward IBP with the node's fixings clamped in ----
    preL = cell(nOps,1); preU = cell(nOps,1);
    actM = cell(nOps,1); unsM = cell(nOps,1);
    auC = cell(nOps,1);  buC = cell(nOps,1);
    unstable = cell(nOps,1);
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
            preL{k} = lb; preU{k} = ub;
            act = (lb >= 0); uns = (lb < 0) & (ub > 0);
            unstable{k} = uns;
            au = zeros(size(lb), precision); bu = zeros(size(lb), precision);
            au(act) = 1;
            dn = ub(uns) - lb(uns);
            au(uns) = ub(uns) ./ dn; bu(uns) = -au(uns) .* lb(uns);
            actM{k} = cast(act, precision); unsM{k} = cast(uns, precision);
            auC{k} = au; buC{k} = bu;
            lb = max(lb,0); ub = max(ub,0);
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
