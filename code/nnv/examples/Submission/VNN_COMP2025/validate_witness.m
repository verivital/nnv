function ok = validate_witness(net, x, lb, ub, Hs, inputSize, inputFormat, needReshape, in_tol, out_tol)
%VALIDATE_WITNESS  Confirm a SAT counterexample is REAL before emitting `sat`.
%   VNN-COMP scores an incorrect/non-replayable verdict at -150 (16x the +10 for a
%   correct one); NNV leaked 19 such "incorrect/missing-CE" points in 2025. This
%   re-evaluates a candidate on the network -- mapping the FLAT witness to the
%   network-input layout with the SAME reshape+permute the runner uses (needReshape)
%   -- and confirms it (1) respects the input box and (2) lands in the unsafe output
%   region. If it fails, the caller must emit `unknown`, never `sat`. Strategy Pillar 2.
%
%   Inputs:
%     net          dlnetwork OR an NNV NN.
%     x            candidate input, FLAT ONNX-order column vector (lb/ub order).
%     lb, ub       input box bounds (flat vectors).
%     Hs           array of HalfSpace (unsafe region).
%     inputSize    network input size.
%     inputFormat  'default' or a dlarray dims label.
%     needReshape  0 | 1,2,3 (mirrors the runner's flat<->net layout).
%     in_tol       input-constraint slack (default 1e-6).
%     out_tol      output-constraint slack (default 1e-4).
%
%   Output: ok -- true only for a validated counterexample.

    if nargin < 8  || isempty(needReshape), needReshape = 0;    end
    if nargin < 9  || isempty(in_tol),      in_tol  = 1e-6;     end
    if nargin < 10 || isempty(out_tol),     out_tol = 1e-4;     end
    ok = false;

    x  = double(x(:));
    lb = double(lb(:)); ub = double(ub(:));

    % (1) Input must lie within the box (near-zero tolerance per the rules).
    if numel(x) ~= numel(lb) || any(x < lb - in_tol) || any(x > ub + in_tol)
        return;
    end

    % (2) Re-evaluate the network on the witness (in the network-input layout).
    try
        y = eval_net(net, x, inputSize, inputFormat, needReshape);
    catch
        return;   % cannot validate -> invalid (caller -> unknown)
    end
    y = double(y(:));
    if isempty(y) || any(~isfinite(y)), return; end

    % (3) Confirm the output lands in some unsafe HalfSpace (G*y <= g).
    for h = 1:numel(Hs)
        m = Hs(h).G * y - Hs(h).g(:);
        if all(m <= out_tol), ok = true; return; end
    end
end

% ---- helpers ----

function y = eval_net(net, x, inputSize, inputFormat, needReshape)
    arr = flat_to_net_arr(x, inputSize, needReshape);   % flat ONNX -> net-input layout
    if isa(net, 'NN')                                   % NNV native net (manifest path)
        y = net.evaluate(arr);
        return;
    end
    fmt = net_format(net, inputFormat, size(arr));
    if ndims(arr) < numel(fmt)
        arr = reshape(arr, [size(arr), ones(1, numel(fmt) - ndims(arr))]);
    end
    y = extractdata(predict(net, dlarray(single(arr), fmt)));
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

function arr = flat_to_net_arr(v, inputSize, needReshape)
    % flat ONNX-order vector -> network-input array (mirrors create_random_examples).
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
