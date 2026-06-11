function [cex, found] = pgd_falsify(net, lb, ub, Hs, inputSize, inputFormat, opts)
%PGD_FALSIFY  Gradient-based (FGSM warm-start + PGD) adversarial falsification.
%   Searches for an input x in the box [lb,ub] whose network output lands in one
%   of the unsafe HalfSpaces in Hs -- i.e. a concrete COUNTEREXAMPLE / SAT witness
%   -- using projected gradient descent on a violation-margin loss, with the
%   gradient obtained by dlnetwork autodiff. This is the gradient-directed
%   complement to NNV's existing random sampling (which found 3x fewer SAT than
%   the field in VNN-COMP 2025); it is meant to be run FIRST, then fall back to
%   random sampling if it fails. See VNNCOMP2026_STRATEGY.md (Pillar 1).
%
%   SOUNDNESS: the returned point is a concrete evaluation, but the caller MUST
%   still validate it with validate_witness() before emitting `sat` (the -150
%   penalty makes a near-miss catastrophic).
%
%   Inputs:
%     net         dlnetwork (required for autodiff). If not a dlnetwork, returns
%                 found=false (the caller should fall back to random sampling).
%     lb, ub      input lower/upper bounds, flattened column vectors.
%     Hs          array of HalfSpace objects (the unsafe output region; a
%                 counterexample y satisfies Hs(h).contains(y) for some h).
%     inputSize   network input size used to reshape x (e.g. [1 784] or [32 32 3]).
%     inputFormat 'default' or a dlarray dims label ("CB","SSCB","BC",...).
%     opts        struct (all optional): n_restarts(20), n_steps(40), lr(0.1),
%                 fgsm(true), seed([] = leave RNG alone).
%
%   Outputs:
%     cex   {x; y} cell (column input x and its output y) on success, else {}.
%     found logical.

    cex = {}; found = false;
    if ~isa(net, 'dlnetwork')
        return;   % no autodiff available -> caller uses random sampling
    end
    if nargin < 7 || isempty(opts), opts = struct(); end
    n_restarts = getfielddef(opts, 'n_restarts', 20);
    n_steps    = getfielddef(opts, 'n_steps',    40);
    lr         = getfielddef(opts, 'lr',         0.1);
    use_fgsm   = getfielddef(opts, 'fgsm',       true);
    if isfield(opts,'seed') && ~isempty(opts.seed), rng(opts.seed); end

    lb = double(lb(:)); ub = double(ub(:));
    n  = numel(lb);
    span = ub - lb;
    fmt = resolve_format(inputFormat, inputSize);

    % Build restart seeds: box corners first (cheap, often violate), then random.
    seeds = cell(1, n_restarts);
    seeds{1} = lb; if n_restarts >= 2, seeds{2} = ub; end
    if n_restarts >= 3, seeds{3} = (lb+ub)/2; end
    for r = 4:n_restarts
        seeds{r} = lb + span.*rand(n,1);
    end

    % Targets: each HalfSpace's (G,g). A counterexample lands in ANY of them; we
    % minimize one target's worst margin per restart, but check ALL after.
    nH = numel(Hs);
    for r = 1:n_restarts
        h = mod(r-1, nH) + 1;          % cycle which unsafe set we steer toward
        [G, g] = halfspace_Gg(Hs(h));
        x = seeds{r};
        for it = 1:n_steps
            dlx = to_dlarray(x, inputSize, fmt);
            [~, grad] = dlfeval(@loss_and_grad, net, dlx, G, g);
            grad = double(extractdata(grad));
            grad = grad(:);
            if use_fgsm && it == 1
                upd = sign(grad);                    % FGSM warm-start step
            else
                gn = norm(grad); if gn < 1e-12, gn = 1; end
                upd = grad / gn;                     % normalized PGD step
            end
            x = x - lr .* span .* upd;               % descend the violation margin
            x = min(max(x, lb), ub);                 % project back into the box
            % Check ALL unsafe sets for a concrete hit.
            y = forward_eval(net, x, inputSize, fmt);
            for hh = 1:nH
                if Hs(hh).contains(double(y(:)))
                    cex = {x; double(y(:))}; found = true; return;
                end
            end
        end
    end
end

% ---- helpers ----

function v = getfielddef(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function fmt = resolve_format(inputFormat, inputSize)
    if nargin >= 1 && ~isempty(inputFormat) && ~strcmp(inputFormat, "default")
        fmt = char(inputFormat);
        if ~endsWith(fmt, 'B'), fmt = [fmt 'B']; end   % ensure a batch dim
        return;
    end
    % default: flat feature input -> CB; image (3D) -> SSCB.
    if isscalar(inputSize) || (numel(inputSize) <= 3 && nnz(inputSize > 1) <= 1)
        fmt = 'CB';
    else
        fmt = 'SSCB';
    end
end

function dlx = to_dlarray(x, inputSize, fmt)
    if isscalar(inputSize)
        xr = reshape(x, [inputSize, 1]);
    else
        xr = reshape(x, inputSize);
    end
    dlx = dlarray(single(xr), fmt);
end

function y = forward_eval(net, x, inputSize, fmt)
    dlx = to_dlarray(x, inputSize, fmt);
    y = extractdata(predict(net, dlx));
end

function [loss, grad] = loss_and_grad(net, dlx, G, g)
    y = predict(net, dlx);
    y = y(:);
    margins = G * y - g(:);          % unsafe set is {y : G*y <= g}; want all <= 0
    loss = max(margins);             % minimize the worst-violated constraint
    grad = dlgradient(loss, dlx);
end

function [G, g] = halfspace_Gg(H)
    % NNV HalfSpace stores the constraint as G*x <= g.
    G = H.G; g = H.g;
end
