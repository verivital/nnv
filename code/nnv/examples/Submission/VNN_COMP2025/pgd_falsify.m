function [cex, found] = pgd_falsify(net, lb, ub, Hs, inputSize, inputFormat, needReshape, opts)
%PGD_FALSIFY  Gradient-based (FGSM warm-start + PGD) adversarial falsification.
%   Searches for an input in the box [lb,ub] whose network output lands in one of
%   the unsafe HalfSpaces in Hs -- a concrete COUNTEREXAMPLE / SAT witness -- using
%   projected gradient descent on a violation-margin loss, with the gradient from
%   dlnetwork autodiff. Gradient-directed complement to NNV's random sampling
%   (which found ~3x fewer SAT than the field in 2025); run FIRST, fall back to
%   random. See VNNCOMP2026_STRATEGY.md (Pillar 1).
%
%   The optimization runs in the NETWORK-INPUT layout: the box and seeds are mapped
%   flat-ONNX -> net-input with the SAME reshape+permute the runner uses
%   (create_random_examples; selected by needReshape), so the forward pass is
%   correct for image/permuted inputs; the returned witness is mapped back to flat
%   ONNX order (matching lb/ub) so write_counterexample / validate_witness consume
%   it directly.
%
%   SOUNDNESS: the returned point is a concrete evaluation, but the caller MUST
%   still validate it with validate_witness() before emitting `sat`.
%
%   Inputs:
%     net          dlnetwork (required for autodiff). Non-dlnetwork -> found=false.
%     lb, ub       input bounds, flat ONNX-order column vectors.
%     Hs           array of HalfSpace (unsafe region; CE satisfies some Hs(h)).
%     inputSize    network input size (e.g. [1 784] or [32 32 3]).
%     inputFormat  'default' or a dlarray dims label ("CB","SSCB",...).
%     needReshape  0 (flat/feature) | 1,2,3 (image permutes; mirrors the runner).
%     opts         struct (optional): n_restarts(20), n_steps(40), lr(0.1),
%                  fgsm(true), seed([]).
%
%   Outputs: cex = {x; y} (flat x, output y) on success else {}; found logical.

    cex = {}; found = false;
    if ~isa(net, 'dlnetwork'); return; end
    if nargin < 7 || isempty(needReshape), needReshape = 0; end
    if nargin < 8 || isempty(opts), opts = struct(); end
    n_restarts = getfielddef(opts, 'n_restarts', 20);
    n_steps    = getfielddef(opts, 'n_steps',    40);
    lr         = getfielddef(opts, 'lr',         0.1);
    use_fgsm   = getfielddef(opts, 'fgsm',       true);
    if isfield(opts,'seed') && ~isempty(opts.seed), rng(opts.seed); end

    lb = double(lb(:)); ub = double(ub(:));
    nIn = numel(lb);
    nsz = net_input_shape(inputSize, needReshape);   % net-input array shape (unpadded)
    fmt = net_format(net, inputFormat, nsz);         % from the net's input layer

    % Box in the network-input layout (element-wise, since flat->net is a permutation).
    LBn = flat_to_net_arr(lb, inputSize, needReshape);
    UBn = flat_to_net_arr(ub, inputSize, needReshape);
    SPn = UBn - LBn;

    % Restart seeds (flat): corners, midpoint, then random.
    seeds = cell(1, n_restarts);
    seeds{1} = lb; if n_restarts >= 2, seeds{2} = ub; end
    if n_restarts >= 3, seeds{3} = (lb+ub)/2; end
    for r = 4:n_restarts, seeds{r} = lb + (ub-lb).*rand(nIn,1); end

    nH = numel(Hs);
    for r = 1:n_restarts
        h = mod(r-1, nH) + 1;                 % steer toward one unsafe set per restart
        [G, g] = halfspace_Gg(Hs(h));
        xn = flat_to_net_arr(seeds{r}, inputSize, needReshape);   % optimize in net space
        for it = 1:n_steps
            dlx = to_dlarray(xn, fmt);
            [~, grad] = dlfeval(@loss_and_grad, net, dlx, G, g);
            gnArr = reshape(double(extractdata(grad)), nsz);      % drop padded batch dim
            if use_fgsm && it == 1
                upd = sign(gnArr);                                % FGSM warm-start
            else
                gn = norm(gnArr(:)); if gn < 1e-12, gn = 1; end
                upd = gnArr / gn;                                 % normalized PGD step
            end
            xn = xn - lr .* SPn .* upd;                           % descend the margin
            xn = min(max(xn, LBn), UBn);                          % project into the box
            y = extractdata(predict(net, to_dlarray(xn, fmt)));
            for hh = 1:nH
                if Hs(hh).contains(double(y(:)))
                    cex = {net_arr_to_flat(xn, inputSize, needReshape); double(y(:))};
                    found = true; return;
                end
            end
        end
    end
end

% ---- helpers ----

function v = getfielddef(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function fmt = net_format(net, inputFormat, nsz)
    if ~isempty(inputFormat) && ~strcmp(inputFormat, "default")
        fmt = char(inputFormat); return;     % trust the caller's label as-is
    end
    L = net.Layers(1);                       % infer from the actual input layer
    if isa(L, 'nnet.cnn.layer.ImageInputLayer')
        fmt = 'SSCB';
    elseif isa(L, 'nnet.cnn.layer.Image3DInputLayer')
        fmt = 'SSSCB';
    elseif isa(L, 'nnet.cnn.layer.FeatureInputLayer') || isa(L, 'nnet.onnx.layer.FeatureInputLayer')
        fmt = 'CB';
    elseif numel(nsz) >= 3
        fmt = 'SSCB';
    else
        fmt = 'CB';
    end
end

function dlx = to_dlarray(arr, fmt)
    % Pad trailing singleton dims so the array rank matches the format label rank.
    if ndims(arr) < numel(fmt)
        arr = reshape(arr, [size(arr), ones(1, numel(fmt) - ndims(arr))]);
    end
    dlx = dlarray(single(arr), fmt);
end

function [loss, grad] = loss_and_grad(net, dlx, G, g)
    y = predict(net, dlx); y = y(:);
    margins = G * y - g(:);          % unsafe set {y : G*y <= g}; drive all <= 0
    loss = max(margins);
    grad = dlgradient(loss, dlx);
end

function [G, g] = halfspace_Gg(H)
    G = H.G; g = H.g;               % NNV HalfSpace: G*x <= g
end

% ---- flat-ONNX <-> network-input layout (mirrors run_vnncomp_instance) ----

function ns = net_input_shape(inputSize, needReshape)
    switch needReshape
        case 1, ns = [inputSize(2) inputSize(1) inputSize(3:end)];
        case 3, ns = [inputSize(3) inputSize(2) inputSize(1)];
        otherwise                                  % 0 and 2 -> inputSize
            if isscalar(inputSize), ns = [inputSize 1]; else, ns = inputSize; end
    end
end

function arr = flat_to_net_arr(v, inputSize, needReshape)
    v = double(v(:));
    switch needReshape
        case 2
            arr = permute(reshape(v, [inputSize(2) inputSize(1) inputSize(3:end)]), [2 1 3]);
        case 3
            arr = permute(reshape(v, inputSize), [3 2 1]);
        case 1
            arr = permute(reshape(v, inputSize), [2 1 3]);
        otherwise
            if isscalar(inputSize), arr = reshape(v, [inputSize 1]); else, arr = reshape(v, inputSize); end
    end
end

function v = net_arr_to_flat(arr, inputSize, needReshape)
    arr = reshape(arr, net_input_shape(inputSize, needReshape));   % drop padded dims
    switch needReshape
        case 2, arr = permute(arr, [2 1 3]);
        case 3, arr = permute(arr, [3 2 1]);
        case 1, arr = permute(arr, [2 1 3]);
    end
    v = reshape(arr, [], 1);
end
