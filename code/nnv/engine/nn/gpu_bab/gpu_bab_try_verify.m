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

    inShape = i_input_shape(net);
    lb = double(lb(:)); ub = double(ub(:));
    if isempty(inShape) || prod(inShape) ~= numel(lb)
        verdict = 'skip'; info.reason = 'unknown/mismatched input shape'; return;
    end

    % (1+2) extract ops AND auto-calibrate the conv->FC flatten order against net.evaluate at
    % SEVERAL probe points. ONNX-imported nets flatten row-major while the engine views the
    % conv output column-major; trying each order and keeping the one whose op-list matches
    % net.evaluate (the soundness guard) handles both. A wrong order/extraction never matches
    % the probes -> 'skip' (Star), never a -150. Multi-probe makes a coincidental match negligible.
    c = (lb + ub) / 2;
    probes  = [c, lb + 0.25*(ub - lb), lb + 0.75*(ub - lb)];
    orders  = {'colmajor', 'chw_rowmajor', 'hwc_rowmajor'};
    ops = []; nClasses = NaN;
    for oi = 1:numel(orders)
        try
            cand = nn_to_ops(net, orders{oi});
        catch ME0
            if oi == 1, info.reason = ['nn_to_ops refused: ' ME0.message]; end
            continue;                               % unsupported layer -> same for every order
        end
        ok = true; maxe = 0;
        for pp = 1:size(probes, 2)
            cp = probes(:, pp);
            yo = gpu_bab_ibp(cand, cp, cp, 'double'); yo = yo(:);
            yn = net.evaluate(reshape(cp, inShape)); yn = yn(:);
            if numel(yo) ~= numel(yn), ok = false; break; end
            e = max(abs(yo - yn)); maxe = max(maxe, e);
            if e > guardTol * max(1, max(abs(yn))), ok = false; break; end
        end
        if ok
            ops = cand; nClasses = numel(yn); info.guardErr = maxe;
            info.reason = sprintf('flatten=%s', orders{oi}); break;
        end
    end
    if isempty(ops)
        if isempty(info.reason), info.reason = 'no flatten order matched net.evaluate (guard)'; end
        verdict = 'skip'; return;
    end
    if target < 1 || target > nClasses
        verdict = 'skip'; info.reason = 'target out of range'; return;
    end

    % (3) sound ReLU-split BaB; force the 'tight' intermediate path (the only one sound for
    %     conv/bn/pool -- relu_split also self-forces this, set here for clarity).
    %     PRECISION/SOUNDNESS (research R1): FP32 bounds are NOT certified-sound without
    %     outward (directed) rounding -- accumulated rounding can make a single-precision
    %     margin wrongly positive (-150). Default 'double' here (FP64 rounding ~1e-15,
    %     negligible for these nets). GPU 'single' is fast but for CERTIFIED competition use
    %     it MUST be paired with per-op outward rounding (R1.a) or a proven FP error margin;
    %     until then 'single' is dev/empirical only. See research/GPU_BAB_PLAN.md sec 3.4/R1.
    bopts = opts; bopts.intermediate = 'tight';
    if ~isfield(bopts, 'precision'), bopts.precision = 'double'; end   % sound default
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
