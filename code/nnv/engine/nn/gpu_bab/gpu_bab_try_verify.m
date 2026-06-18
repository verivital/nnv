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
        % i_input_shape could not determine a matching shape (e.g. a leading SequenceInput/
        % PlaceholderLayer from a 'BCT' ONNX import of a flat FC net -- relusplitter/safenlp).
        % Fall back to a FLAT [n 1] input, valid for feature/FC nets. The orientation guard
        % below REQUIRES net.evaluate(reshape(c, inShape)) to match the op-list at several
        % probes, so a wrong shape (e.g. a genuine image net) fails the guard -> 'skip'
        % (sound: a mis-shaped query never yields a verdict).
        inShape = [numel(lb) 1];
    end

    % (1+2) extract ops AND auto-calibrate the conv->FC flatten order against net.evaluate at
    % SEVERAL probe points. ONNX-imported nets flatten row-major while the engine views the
    % conv output column-major; trying each order and keeping the one whose op-list matches
    % net.evaluate (the soundness guard) handles both. A wrong order/extraction never matches
    % the probes -> 'skip' (Star), never a -150. Multi-probe makes a coincidental match negligible.
    c = (lb + ub) / 2;
    % SOUNDNESS (audit finding): the guard must DISTINGUISH a wrong conv->FC flatten
    % permutation. COLLINEAR probes lb+s*(ub-lb) are BLIND to column permutations on a uniform
    % L-inf box -- for any permutation P, (P-I)*ones=0, so the s-dependent error cancels and
    % the standard L-inf query (ub-lb uniform) admits a wrong flatten that passes the guard
    % => false robust. Use NON-uniform, distinct-per-coordinate probe points (low-discrepancy
    % fractional ramps): the probe images are spatially varying, so the conv output differs at
    % permuted positions and a wrong order cannot match net.evaluate. Deterministic (no rng).
    n = numel(lb); idx = (1:n)';
    f1 = mod(idx * 0.6180339887498949, 1);
    f2 = mod(idx * 1.3247179572447460 + 0.37, 1);
    f3 = mod(idx * 0.7548776662466927 + 0.11, 1);
    probes  = [c, lb + f1.*(ub-lb), lb + f2.*(ub-lb), lb + f3.*(ub-lb), lb + (1-f1).*(ub-lb)];
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
    % `target` may be a vnnlib PROPERTY (cell/struct with .Hg) instead of a class index. In
    % that case derive the target from the net's OWN center prediction and VERIFY the spec is
    % argmax-robustness for it (every unsafe halfspace a positive multiple of e_target - e_j,
    % g<=0). Then a 'robust' verdict (argmax==target on the whole box) PROVABLY implies the
    % vnnlib unsat (the output avoids every unsafe halfspace). A spec that does not match this
    % exact pattern -> 'skip' (sound: we never emit a verdict about a property we did not prove).
    if isstruct(target) || iscell(target)
        % FAIL CLOSED: this is an additive pre-check, so any failure deriving/validating the
        % target (a net.evaluate that throws, a malformed/empty property) must return 'skip'
        % (-> caller runs Star), never error or emit a verdict.
        try
            yc = net.evaluate(reshape(c, inShape)); yc = yc(:);
            [~, tIdx] = max(yc);
            okSpec = i_is_argmax_spec(target, tIdx, nClasses);
        catch ME
            verdict = 'skip'; info.reason = ['target derivation failed: ' ME.message]; return;
        end
        if ~okSpec
            verdict = 'skip'; info.reason = 'spec not argmax-robustness for the center prediction'; return;
        end
        target = tIdx;
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
    % R1 SOUNDNESS (audit finding): a caller-supplied precision='single' is NOT sound (FP32
    % rounding can make a margin spuriously positive => false robust). FORCE 'double' unless
    % the caller EXPLICITLY opts into unsound single (dev / GPU speed; a SOUND single path
    % needs per-op outward rounding, R1.a, not yet implemented).
    if isfield(opts, 'allowUnsoundSingle') && isequal(opts.allowUnsoundSingle, true)
        if ~isfield(bopts, 'precision'), bopts.precision = 'single'; end
    else
        bopts.precision = 'double';
    end
    % a negative leaf margin would certify nodes whose true class does NOT strictly dominate
    % (false robust); clamp to >= 0 (sound-or-unknown).
    if isfield(bopts, 'margin') && bopts.margin < 0, bopts.margin = 0; end
    if ~isfield(bopts, 'maxNodes'),  bopts.maxNodes  = 300;      end
    % per-node alpha+beta CROWN tightening (env override for dev/tuning; sound for any value --
    % alpha in [0,1], beta>=0 are valid relaxations). Lets the BaB cross the convex barrier the
    % fixed-slope bound can't, on the hardest deep nodes.
    % only override when the parsed env is finite + in range (str2double of unset/non-numeric is NaN;
    % a NaN alphaIter/betaIter would propagate through max(alphaIter,betaIter) and silently disable the knob)
    if ~isfield(bopts, 'alphaIter'), v = str2double(getenv('NNV_BAB_ALPHA_ITERS')); if isfinite(v) && v >= 0, bopts.alphaIter = v; end; end
    if ~isfield(bopts, 'betaIter'),  v = str2double(getenv('NNV_BAB_BETA_ITERS'));  if isfinite(v) && v >= 0, bopts.betaIter  = v; end; end
    if ~isfield(bopts, 'alphaLr'),   v = str2double(getenv('NNV_BAB_ALPHA_LR'));    if isfinite(v) && v >  0, bopts.alphaLr   = v; end; end
    % DEVICE (opts.device 'cpu' default | 'gpu'): run the whole BaB on the GPU by moving the ops'
    % weight tensors + the input box to gpuArray ONCE here. The IBP/CROWN passes are gpuArray-
    % overloaded pure linear algebra, so all heavy compute + the per-node intermediate bounds stay
    % RESIDENT on device; only tiny scalars (the certified test, the split-neuron index, a witness)
    % cross back -- no per-iteration bound transfers. GPU implies single precision (T4 FP64 is
    % ~1/16 rate), so it requires the explicit allowUnsoundSingle opt-in (FP32 is dev/empirical,
    % not certified-sound without outward rounding -- R1). CPU path is unchanged (the default).
    if isfield(opts, 'device') && strcmp(opts.device, 'gpu')
        if gpuDeviceCount < 1
            verdict = 'skip'; info.reason = 'device=gpu requested but no GPU available'; return;
        end
        if ~(isfield(opts, 'allowUnsoundSingle') && isequal(opts.allowUnsoundSingle, true))
            verdict = 'skip'; info.reason = 'device=gpu needs allowUnsoundSingle (FP32 not certified-sound)'; return;
        end
        bopts.precision = 'single';
        ops = i_ops_to_gpu(ops);
        lb = gpuArray(single(lb)); ub = gpuArray(single(ub));
        info.device = 'gpu';
    end
    % ---- M1: sound-FP32 root pre-check EMIT --------------------------------------------------
    % When sound-FP32 is enabled (NNV_SOUND_FP32_TIGHT set) on the single-precision path,
    % gpu_bab_crown_tight returns a PROVABLY SOUND lower bound on the spec margins (every CROWN
    % bound outward-widened by a per-op running-error radius; validated G1 0/99, PR #385). So
    % min(margins) > 0 at the ROOT certifies argmax-robustness with NO branch-and-bound and NO
    % FP64 confirm -- a sound ~9x emit for the root-certifiable tier. In this mode we do NOT fall
    % through to the (not-yet-sound-FP32) batched BaB for a verdict (its per-node spec_dag is not
    % yet hardened -- milestone M3b): a non-certifying root returns 'unknown' so the caller runs
    % the sound fallback. FAIL-CLOSED: any crown_tight/vmag error -> 'unknown', never a verdict.
    info.soundFP32 = false;
    soundFP32 = strcmp(bopts.precision, 'single') && ~isempty(getenv('NNV_SOUND_FP32_TIGHT'));
    if soundFP32
        rootCert = false;
        try
            Cspec = i_argmax_spec_C(target, nClasses, lb);          % rows e_target - e_j (== relu_split.m)
            mrg = gpu_bab_crown_tight(ops, lb, ub, Cspec, 'single');% sound lower bound on C*f(x)
            mrg = gather(double(mrg(:)));
            rootCert = ~isempty(mrg) && all(isfinite(mrg)) && all(mrg > 0);
        catch ME1
            info.reason = [info.reason '; sound-FP32 root pre-check error: ' ME1.message];
        end
        if rootCert
            verdict = 'robust'; info.soundFP32 = true; info.nodes = 0;
            info.reason = sprintf('%s; sound-FP32 root pre-check robust (min margin %.4g)', info.reason, min(mrg));
        else
            verdict = 'unknown';
            info.reason = [info.reason '; sound-FP32 root did not certify (deep BaB pending M3b)'];
        end
        return;
    end

    % Engine: opts.engine='batched' uses gpu_bab_relu_split_batched (batched-DFS, processes a
    % whole BaB-node frontier per kernel -- GPU-saturating, for FC + sequential conv); default
    % is the serial gpu_bab_relu_split. The batched bounding is bound-for-bound identical
    % (sound), so the verdict is the same; only throughput differs. (batched refuses 'add'/
    % maxpool DAG ops -> it errors -> caught below -> 'skip', never an unsound verdict.)
    useBatched = isfield(opts, 'engine') && strcmp(opts.engine, 'batched');
    try
        if useBatched
            [status, binfo] = gpu_bab_relu_split_batched(ops, lb, ub, target, nClasses, bopts);
        else
            [status, binfo] = gpu_bab_relu_split(ops, lb, ub, target, nClasses, bopts);
        end
    catch ME
        % a batched-engine refusal (e.g. residual/maxpool DAG) or any engine error -> skip,
        % so the caller falls back to Star (additive: GPU-BaB never harms).
        verdict = 'skip'; info.reason = ['engine error: ' ME.message]; return;
    end
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

function C = i_argmax_spec_C(target, K, like_arr)
% Argmax-robustness spec matrix, IDENTICAL to gpu_bab_relu_split.m's construction: row j
% (j ~= target) = e_target - e_j, so C*f(x) = f_target - f_j and 'robust' iff every margin > 0.
% Built 'like' the input box so it inherits gpuArray/single on the GPU path.
    C = -eye(K, 'like', like_arr);
    C(:, target) = C(:, target) + 1;
    C(target, :) = [];
end

function ok = i_is_argmax_spec(prop, target, K)
% SOUND check: is the vnnlib property EXACTLY argmax-robustness for class `target`? True iff
% every unsafe halfspace {y: G*y <= g} is a positive multiple of e_target - e_j (g <= 0) for
% some j ~= target. Then any output in the unsafe region has y_target <= y_j (+ a non-positive
% margin) for some j => argmax(y) ~= target. So proving argmax==target on the whole box
% ('robust') implies the output avoids EVERY unsafe halfspace => the vnnlib instance is UNSAT.
% (The gpu-bab proves target dominates ALL classes, >= what the spec's listed j's require, so
% it is sound-but-conservative if the spec lists a subset.) Anything else -> false -> skip.
    ok = false;
    if iscell(prop)
        if isempty(prop), return; end          % empty cell -> not a spec -> skip
        prop = prop{1};
    end
    if ~isstruct(prop) || ~isfield(prop, 'Hg'), return; end   % malformed -> fail closed
    Hg = prop.Hg;
    if isempty(Hg), return; end
    for i = 1:numel(Hg)
        G = double(Hg(i).G); g = double(Hg(i).g);
        if size(G,1) ~= 1 || size(G,2) ~= K, return; end     % single row over the K outputs
        if any(g > 1e-9), return; end                        % margin must be non-positive
        nz = find(abs(G) > 1e-12);
        if numel(nz) ~= 2 || ~ismember(target, nz), return; end   % exactly target vs one other
        other = nz(nz ~= target);
        if G(target) <= 0 || G(other) >= 0, return; end      % target +, the other -
        if abs(abs(G(target)) - abs(G(other))) > 1e-6 * max(1, abs(G(target))), return; end  % ~ e_target - e_j
    end
    ok = true;
end

function ops = i_ops_to_gpu(ops)
% Move each op's numeric DATA tensors (weights/bias/scale/shift) to gpuArray(single) so the
% bound passes run on device. Leaves dlconv hyperparameters (pool/stride/pad/dil) and metadata
% (type/src/inputs/shapes) on the host -- they are scalars/small index vectors dlconv expects on
% the host, and moving them would not help. Done ONCE per query (not per BaB node).
    for k = 1:numel(ops)
        o = ops{k};
        for f = {'W', 'b', 'scale', 'shift'}
            if isfield(o, f{1}) && isnumeric(o.(f{1}))
                o.(f{1}) = gpuArray(single(o.(f{1})));
            end
        end
        ops{k} = o;
    end
end
