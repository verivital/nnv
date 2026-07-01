function status = i_vit_sound_unsat(category, onnx, lb, ub, prop)
%I_VIT_SOUND_UNSAT  Sound ViT 'unsat' (verified-robust) via the native ViTCrown
%   linear-bound verifier. Returns 1 ONLY when robustness is PROVEN against the exact
%   vnnlib box + output halfspaces; returns 2 (defer) on anything else (not the 32x32
%   CIFAR ViT spec / weight bundle missing / shape mismatch / not robust / exception),
%   so the caller falls through to the existing falsify->reach({})->unknown path with
%   ZERO change in behaviour.
%
%   Verdict logic (sound, matches verify_specification): the vnnlib output property is
%   an OR of "unsafe" halfspaces {G_k*y <= g_k}. The instance is robust iff the reach
%   set misses EVERY unsafe halfspace. A halfspace k is provably missed if SOME row r
%   is always violated over the input box, i.e. min_box(G_kr*y) > g_kr. ViTCrown's
%   optimizeAlpha returns a SOUND lower bound of (C*logits) over the box, so with
%   C = G_all this gives min_box(G_kr*y); subtract g and require, per halfspace, that
%   its best row exceeds 0. (No argmax/label assumption -- driven off the authoritative
%   prop halfspaces.) NO LP solver, NO external verifier (see ViTCrown.m provenance).
%
%   Inputs: category (string), onnx (path), lb/ub (vnnlib input box, already the
%   normalized c-major 3*32*32 layout ViTCrown expects -- fed AS-IS, no re-normalize),
%   prop (cell/array of HalfSpace specs from load_vnnlib).
    status = 2;                                   % default: defer / unknown
    try
        if ~contains(category, "vit"), return; end
        if iscell(lb) || iscell(ub), return; end  % single normalized box only
        lb = lb(:); ub = ub(:);
        if numel(lb) ~= 3*32*32 || numel(ub) ~= numel(lb), return; end

        bundle = [char(onnx) '.vitbundle.mat'];    % weights-only .mat (built in prepare)
        if exist(bundle, 'file') ~= 2, return; end % no bundle -> sound fall-through
        M   = ViTReach.load(bundle);
        ops = ViTCrown.toOps(M);
        if ops{1}.dim ~= numel(lb), return; end    % layout mismatch -> defer

        % Match verify_specification EXACTLY: the spec is property{1}.Hg, an array of
        % DISJUNCT unsafe halfspaces (OR). Robust iff the reach set misses EVERY
        % disjunct. Disjunct cp = {y : G_cp*y <= g_cp} (G_cp possibly multi-row = AND).
        if iscell(prop), Hg = prop{1}.Hg; else, Hg = prop.Hg; end
        np = numel(Hg);
        if np < 1, return; end
        Gall = []; gall = []; disjOf = [];          % stack rows, tag with disjunct index
        for cp = 1:np
            Gc = double(Hg(cp).G); gc = double(Hg(cp).g(:));
            if size(Gc,2) ~= 10, return; end        % spec not over the 10 logits -> defer
            Gall = [Gall; Gc]; gall = [gall; gc];               %#ok<AGROW>
            disjOf = [disjOf; cp*ones(size(Gc,1),1)];           %#ok<AGROW>
        end

        nIter = i_envnum('NNV_VIT_ALPHA_ITERS', 100);
        lr    = i_envnum('NNV_VIT_ALPHA_LR', 0.03);
        [cl, cu] = ViTCrown.refineBounds(ops, lb, ub, 2);          % CROWN intermediate bounds
        rowLB = ViTCrown.optimizeAlpha(ops, lb, ub, cl, cu, Gall, ...
                    struct('nIter', nIter, 'lr', lr)) - gall;       % sound min_box(G*y) - g

        % CRITICAL SOUNDNESS GUARD: a non-finite bound (NaN/Inf, e.g. CROWN overflow on
        % some nets) is NOT a proof of anything -> defer to unknown. Without this, the
        % loop below is unsound: max(NaN) <= tol evaluates to false, so 'robust' would
        % never be cleared and a FALSE 'unsat' would be emitted (matches run_vit_crown's
        % safe all(mg>0) polarity, where NaN -> false -> unknown).
        if ~all(isfinite(rowLB)), return; end

        % disjunct cp is MISSED iff some row is always violated: max_r(rowLB) > 0.
        % Robust iff EVERY disjunct missed (min over disjuncts > tol). For the standard
        % argmax spec (single-row disjuncts), this is all(margins > 0) = the binding margin.
        tol = 1e-7;                                % FP guard above bound-computation noise
        robust = true;
        for cp = 1:np
            if ~(max(rowLB(disjOf == cp)) > tol), robust = false; break; end   % require positive proof
        end
        if robust, status = 1; end                 % SOUND unsat (verified robust)
    catch ME
        fprintf('ViT sound-unsat skipped (%s) -> fall through to falsify/unknown\n', ME.message);
        status = 2;
    end
end

function v = i_envnum(name, dflt)
    v = str2double(getenv(name));
    if isnan(v) || v <= 0, v = dflt; end
end
