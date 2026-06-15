function [margins, preL, preU] = gpu_bab_crown_spec(ops, x_lb, x_ub, C, precision, fixings)
% GPU_BAB_CROWN_SPEC  Sound CROWN lower bound on a linear output spec C*f(x),
%   batched over input sub-domains -- the GPU-BaB bounding primitive.
%
%   [margins, preL, preU] = GPU_BAB_CROWN_SPEC(ops, x_lb, x_ub, C, precision, fixings)
%     ops        : affine/relu op list (see nn_to_ops)
%     x_lb, x_ub : input box(es), n-by-B  (B sub-domains = the GPU BATCH dimension)
%     C          : spec matrix, nSpec-by-nOut (same spec applied to every sub-box)
%     precision  : 'single' (default) | 'double'
%     fixings    : OPTIONAL cell(nOps,1); for a relu op k, fixings{k} is a dim_k-by-B
%                  matrix of per-node ReLU-split clamps (-1 inactive / 0 free / +1 active).
%                  This is what lets the B columns be ReLU-split BaB NODES that share the
%                  SAME input box but partition the unstable neurons -- the clamp pins a
%                  fixed neuron's pre-activation (active: l:=max(l,0); inactive: u:=min(u,0))
%                  exactly as the serial gpu_bab_relu_split does, and propagates it forward.
%                  Omit (or pass {}) for plain input-split bounding (gpu_bab_bab).
%     margins    : nSpec-by-B;  margins(i,k) <= min_{x in box k} C(i,:)*f(x)
%     preL, preU : OPTIONAL cell(nOps,1); the (clamped) pre-activation bounds at each relu
%                  (dim_k-by-B), for the caller's split-neuron selection / bound reuse.
%
%   For robustness: C rows = e_true - e_j (j ~= true); all margins(:,k) > 0 proves
%   sub-box k robust. Bounding a SPEC (not per-output ranges) is what lets a single
%   backward pass certify the property directly -- the standard CROWN verification
%   form, and the shape beta-CROWN batches over branch-and-bound sub-domains.
%
%   GPU-BaB core: the backward coefficient tensor A is nSpec-by-dim-by-B (one page
%   per sub-box, because each box's ReLU relaxation differs), propagated with
%   PAGEMTIMES -- one batched kernel over all sub-boxes, no host loop. Pure MATLAB;
%   pass x_lb/x_ub and the ops' weights as gpuArray to run on device (pagemtimes is
%   gpuArray-overloaded). The bound columns are the GPU batch dimension.
%
%   SOUNDNESS: identical sign-aware ReLU relaxation as gpu_bab_crown (forward IBP for
%   per-neuron pre-activation bounds; backward linear bounds; for a LOWER bound a
%   positive coeff takes the lower line, a negative coeff the upper line), concretised
%   by the box dual norm. For every x in box k, C*f(x) >= margins(:,k).

    if nargin < 5 || isempty(precision)
        precision = 'single';
    end
    if ~ismember(precision, {'single','double'})
        error('gpu_bab_crown_spec:precision', "precision must be 'single' or 'double'");
    end
    if nargin < 6, fixings = {}; end   % optional per-relu BaB node clamps (-1/0/+1), dim_k-by-B

    B = size(x_lb, 2);
    nSpec = size(C, 1);
    nOps = numel(ops);

    % ---- forward IBP (batched): pre-activation bounds at each ReLU ----
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
                % Per-node ReLU-split clamp (BaB): a fixed-active neuron has z>=0 so its
                % lower pre-activation bound rises to 0 (the neuron is now stable-on); a
                % fixed-inactive neuron has z<=0 so its upper bound drops to 0 (stable-off).
                % Pinning both propagates the split forward exactly as gpu_bab_relu_split's
                % i_crown_clamped does, just batched over the B node columns.
                if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                    fx = fixings{k};                         % dim_k x B
                    lb(fx == 1)  = max(lb(fx == 1),  0);     % active fix:   z >= 0
                    ub(fx == -1) = min(ub(fx == -1), 0);     % inactive fix: z <= 0
                end
                preL{k} = lb; preU{k} = ub;
                lb = max(lb, 0); ub = max(ub, 0);
            otherwise
                error('gpu_bab_crown_spec:op', 'Unsupported op "%s".', op.type);
        end
    end

    % ---- backward (batched): lower bound on C*f ----
    A = repmat(cast(C, precision), 1, 1, B);     % nSpec x nOut x B
    d = zeros(nSpec, B, precision);
    for k = nOps:-1:1
        op = ops{k};
        switch op.type
            case 'affine'
                W = cast(op.W, precision); b = cast(op.b(:), precision);
                d = d + reshape(pagemtimes(A, b), nSpec, B);   % A*b per page
                A = pagemtimes(A, W);                          % nSpec x prevdim x B
            case 'relu'
                l = preL{k}; u = preU{k};                      % dim x B
                [au, bu, al] = i_relu_relax(l, u, precision);  % each dim x B
                dim = size(l, 1);
                Apos = max(A, 0); Aneg = min(A, 0);
                % lower bound: positive coeff -> lower line (al, no intercept);
                %              negative coeff -> upper line (au, bu)
                d = d + reshape(pagemtimes(Aneg, reshape(bu, dim, 1, B)), nSpec, B);
                A = Apos .* reshape(al, 1, dim, B) + Aneg .* reshape(au, 1, dim, B);
            otherwise
                error('gpu_bab_crown_spec:op', 'Unsupported op "%s".', op.type);
        end
    end

    n = size(x_lb, 1);
    Apos = max(A, 0); Aneg = min(A, 0);
    lbcol = reshape(cast(x_lb, precision), n, 1, B);
    ubcol = reshape(cast(x_ub, precision), n, 1, B);
    margins = reshape(pagemtimes(Apos, lbcol), nSpec, B) ...
            + reshape(pagemtimes(Aneg, ubcol), nSpec, B) + d;   % min over each box
end

function [au, bu, al] = i_relu_relax(l, u, precision)
% Per-neuron ReLU relaxation over pre-activation [l,u] (elementwise, dim x B):
%   stable-on (l>=0): identity. stable-off (u<=0): zero. unstable: upper line
%   au*z+bu (au=u/(u-l), bu=-au*l>=0); lower line al*z (al in {0,1}, min-area).
    au = zeros(size(l), precision);
    bu = zeros(size(l), precision);
    al = zeros(size(l), precision);
    act = (l >= 0);
    au(act) = 1; al(act) = 1;
    unst = (l < 0) & (u > 0);
    dn = u(unst) - l(unst);
    au(unst) = u(unst) ./ dn;
    bu(unst) = -au(unst) .* l(unst);
    al(unst) = cast(u(unst) >= -l(unst), precision);
end
