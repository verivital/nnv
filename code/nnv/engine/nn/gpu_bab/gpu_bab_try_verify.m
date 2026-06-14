function [verdict, info] = gpu_bab_try_verify(net, lb, ub, target, opts)
% GPU_BAB_TRY_VERIFY  Sound, ADDITIVE GPU-BaB pre-check for argmax-robustness.
%
%   [verdict, info] = GPU_BAB_TRY_VERIFY(net, lb, ub, target, opts) attempts to decide
%   whether every input in the box [lb, ub] keeps argmax(net(x)) == target, using the
%   LP-free GPU-BaB engine (CROWN/IBP + ReLU-split branch-and-bound). It is designed to be
%   wired in FRONT of NNV's Star reachability as a fast pre-check that can ONLY help:
%
%     verdict : 'robust'  -- PROVEN: a sound CROWN/BaB lower bound certifies all margins > 0
%               'unsafe'  -- a concrete counterexample, RE-CONFIRMED against net.evaluate
%               'unknown' -- engine could not decide within budget -> caller runs Star
%               'skip'    -- architecture unsupported OR the soundness guard failed
%                            -> caller runs Star (GPU-BaB contributes nothing, never harm)
%
%   SOUNDNESS (why this can never produce a wrong verdict / -150):
%     (1) nn_to_ops REFUSES any layer it cannot bound soundly (errors) -> 'skip'.
%     (2) ORIENTATION GUARD: the op-list is bounded over a flat [H W C] column-major box;
%         reach may feed a permuted box (needReshape). We require the op-list evaluation at
%         the box CENTER to equal net.evaluate(center) to tolerance, i.e. the op-list bounds
%         the SAME function the verdict is about. A transposed/mis-ordered box fails this for
%         a generic center -> 'skip'. So 'robust' is a bound on the right function.
%     (3) gpu_bab_relu_split is sound-or-unknown (robust = every BaB leaf certified).
%     (4) 'unsafe' is taken ONLY when a concrete witness ALSO misclassifies under
%         net.evaluate (double-checked here) -> a real counterexample for the real net.
%
%   opts (all optional): .precision 'single'(def)|'double'  .maxNodes 300  .guardTol 1e-4
%         .nSample 16  .cexEvery 25  (forwarded to gpu_bab_relu_split; intermediate forced 'tight')

    info = struct('reason', '', 'nodes', 0, 'cex', [], 'guardErr', NaN);
    if nargin < 5, opts = struct(); end
    guardTol = i_optget(opts, 'guardTol', 1e-4);

    % (1) extract ops; refuse-on-unsupported -> skip (caller's Star path)
    try
        ops = nn_to_ops(net);
    catch ME
        verdict = 'skip'; info.reason = ['nn_to_ops refused: ' ME.message]; return;
    end

    inShape = i_input_shape(net);
    lb = double(lb(:)); ub = double(ub(:));
    if isempty(inShape) || prod(inShape) ~= numel(lb)
        verdict = 'skip'; info.reason = 'unknown/mismatched input shape'; return;
    end

    % (2) orientation soundness guard: op-list eval == net.evaluate at the box center
    c  = (lb + ub) / 2;
    yo = gpu_bab_ibp(ops, c, c, 'double'); yo = yo(:);
    yn = net.evaluate(reshape(c, inShape));  yn = yn(:);
    if numel(yo) ~= numel(yn)
        verdict = 'skip'; info.reason = 'output-size mismatch'; return;
    end
    info.guardErr = max(abs(yo - yn));
    if info.guardErr > guardTol * max(1, max(abs(yn)))
        verdict = 'skip'; info.reason = sprintf('orientation guard failed (%.2e)', info.guardErr); return;
    end
    nClasses = numel(yn);
    if target < 1 || target > nClasses
        verdict = 'skip'; info.reason = 'target out of range'; return;
    end

    % (3) sound ReLU-split BaB; force the 'tight' intermediate path (the only one sound for
    %     conv/bn/pool -- relu_split also self-forces this, set here for clarity).
    bopts = opts; bopts.intermediate = 'tight';
    if ~isfield(bopts, 'precision'), bopts.precision = 'single'; end
    if ~isfield(bopts, 'maxNodes'),  bopts.maxNodes  = 300;      end
    [status, binfo] = gpu_bab_relu_split(ops, lb, ub, target, nClasses, bopts);
    info.nodes = binfo.nodes;

    % (4) map verdict; re-confirm any counterexample against net.evaluate
    switch status
        case 'robust'
            verdict = 'robust';
        case 'unsafe'
            yc = net.evaluate(reshape(double(binfo.cex), inShape)); yc = yc(:);
            [~, pred] = max(yc);
            if pred ~= target
                verdict = 'unsafe'; info.cex = binfo.cex;        % confirmed real witness
            else
                verdict = 'unknown'; info.reason = 'BaB cex not confirmed by net.evaluate';
            end
        otherwise
            verdict = 'unknown';
    end
end

function s = i_input_shape(net)
% Input [H W C] for an image net, or [n 1] for a feature/flat net, from the first layer.
    s = [];
    if isempty(net.Layers), return; end
    L = net.Layers{1};
    if isprop(L, 'InputSize') && ~isempty(L.InputSize)
        is = double(L.InputSize(:)');
        if numel(is) == 1, s = [is 1]; else, s = is; end
    end
end

function v = i_optget(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
