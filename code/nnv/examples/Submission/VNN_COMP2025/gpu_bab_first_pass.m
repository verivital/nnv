function res = gpu_bab_first_pass(category, onnx, vnnlib, opts)
% GPU_BAB_FIRST_PASS  Offline, READ-ONLY GPU-BaB first-pass on one vnncomp instance.
%
%   res = GPU_BAB_FIRST_PASS(category, onnx, vnnlib, opts) loads the instance with the
%   SAME loaders run_vnncomp_instance uses (load_vnncomp_network + load_vnnlib), builds the
%   net-order input box exactly as create_input_set would (replicating the needReshape
%   permute), takes target = argmax(net(box center)) (the nominal class for an argmax-
%   robustness benchmark), and runs the sound gpu_bab_try_verify pre-check. It returns the
%   verdict ONLY -- it changes no production result and is meant to be COMPARED against gold
%   (Star) to measure GPU-BaB coverage + agreement. A gpu='robust' that disagrees with a
%   gold violation is a RED FLAG to investigate (extraction/orientation), never a verdict.
%
%   res fields: verdict ('robust'|'unsafe'|'unknown'|'skip'), target, nClasses, guardErr,
%               nodes, reason. opts forwarded to gpu_bab_try_verify (maxNodes, precision...).

    if nargin < 4, opts = struct(); end
    res = struct('verdict','error','target',NaN,'nClasses',NaN,'guardErr',NaN,'nodes',0,'reason','');
    try
        [~, nnvnet, needReshape, ~, inputSize, ~, ~, ~] = load_vnncomp_network(category, onnx, vnnlib);
        property = load_vnnlib(vnnlib);
        lb = property.lb; ub = property.ub;
        if iscell(lb), lb = lb{1}; ub = ub{1}; end           % first input set (first-pass)

        [blb, bub] = i_netorder_box(lb, ub, inputSize, needReshape);
        inShape = i_inshape(nnvnet, numel(blb));
        c  = (double(blb) + double(bub)) / 2;
        yc = nnvnet.evaluate(reshape(c, inShape)); yc = yc(:);
        [~, tgt] = max(yc);

        [verdict, info] = gpu_bab_try_verify(nnvnet, blb, bub, tgt, opts);
        res.verdict  = char(verdict);
        res.target   = tgt;
        res.nClasses = numel(yc);
        res.guardErr = info.guardErr;
        res.nodes    = info.nodes;
        res.reason   = char(info.reason);
    catch ME
        res.verdict = 'error'; res.reason = ME.message;
    end
end

function [blb, bub] = i_netorder_box(lb, ub, inputSize, needReshape)
% Reproduce create_input_set's box -> net-order [H W C] mapping, then flatten column-major.
    if iscell(inputSize), inputSize = inputSize{1}; end
    is_feature = isscalar(inputSize) || (numel(inputSize) <= 3 && nnz(inputSize > 1) <= 1);
    if is_feature
        blb = double(lb(:)); bub = double(ub(:)); return;    % flat feature input -- no permute
    end
    L = reshape(double(lb), inputSize); U = reshape(double(ub), inputSize);
    if needReshape == 1
        L = permute(L, [2 1 3]);            U = permute(U, [2 1 3]);
    elseif needReshape == 2
        ns = [inputSize(2) inputSize(1) inputSize(3:end)];
        L = permute(reshape(double(lb), ns), [2 1 3 4]);
        U = permute(reshape(double(ub), ns), [2 1 3 4]);
    elseif needReshape == 3
        L = permute(L, [3 2 1]);            U = permute(U, [3 2 1]);
    end
    blb = L(:); bub = U(:);
end

function s = i_inshape(nnvnet, n)
    s = [];
    if ~isempty(nnvnet.Layers)
        L = nnvnet.Layers{1};
        if isprop(L, 'InputSize') && ~isempty(L.InputSize)
            is = double(L.InputSize(:)');
            if numel(is) == 1, s = [is 1]; else, s = is; end
        end
    end
    if isempty(s), s = [n 1]; end
end
