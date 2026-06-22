function [margins, preL, preU, scoreCell, soundFP32] = gpu_bab_crown_spec_dag(ops, x_lb, x_ub, C, precision, fixings, rootBounds, alphaCell, vmag, cSel)
% BINDING-SPEC BRANCHING (10th arg cSel, optional): when non-empty, a per-NODE single
% binding-spec coefficient (nSpec_eff x nOut x B, typically nSpec_eff=1 -- the least-avoided
% disjunct's row for each node) REPLACES the full C as the backward seed. The 4th output
% scoreCell then reflects the BaBSR sensitivity of ONLY the binding disjunct (the spec that is
% actually blocking certification), so i_pick_splits branches on the neuron that helps the
% binding spec instead of one helping an already-avoided spec.
% ⚠ WHEN cSel IS PASSED, the returned `margins` (and preL/preU) ARE the SELECTED-row bound, NOT the
%   full-C spec margins -- the backward seed becomes cSel (see below), so they correspond to cSel's
%   rows only. The caller MUST use ONLY the 4th output (scoreCell) from a cSel call and DISCARD its
%   margins; CERTIFICATION must always use the full-C margins from a SEPARATE cSel-absent call. (The
%   sole caller, gpu_bab_relu_split_batched, does exactly this: `[~,~,~,scB] = ...spec_dag(...,cSel)`
%   for the score, and certifies on the full-C `margins` from i_bound_batch.) So this is "bound-
%   invariant" at the SYSTEM level (split CHOICE only; never the certified bound) but NOT a no-op on
%   this function's `margins` output -- do not read cSel-path margins as full-spec.
% Absent (default cSel={}) -> byte-identical to the prior full-spec path.
% SOUND-FP32 (M3b): optional `vmag` (cell(nOps+1,1) of DOUBLE per-op output value-magnitude
% majorants from gpu_bab_ibp(...,'double'), computed ONCE at the BaB root since the input box is
% fixed across nodes). When supplied on the single-precision path, the batched backward accumulates
% a transient `derr` (nSpec x B) bounding the FP32 rounding error of every backward op, and the
% final margins are widened OUTWARD by (rad + derr) so `margins_single <= margins_double <= true`
% (monotone-sound -> can EMIT, no FP64 confirm). 5th output `soundFP32` = (single && vmag applied),
% the fail-closed gate the caller MUST check before trusting a single-precision 'robust'. Absent vmag
% (or double) -> byte-identical to the prior unsound screen (never emits). See FP32_SOUND_RECIPE §M3b.
% NOTE: i_gamma_sd / i_outward_rad_sd below MUST stay identical to gpu_bab_crown_tight.m's
% i_gamma / i_outward_rad (M0 will factor them to one shared file; until then, diff char-for-char).
% AMORTIZED alpha-CROWN: optional alphaCell (cell(nOps,1) of per-relu lower-slope vectors,
% dim_k x 1, in [0,1]) lets the caller bound a whole BaB frontier with FIXED root-optimized
% slopes -- no autodiff, no per-node gradient tape -> large frontier. Unset -> min-area (default,
% bound-for-bound unchanged). Any alpha in [0,1] is a sound lower ReLU slope, so this is sound.
% GPU_BAB_CROWN_SPEC_DAG  Sound CROWN lower bound on a linear output spec C*f(x), batched
%   over B node columns, for feedforward conv nets INCLUDING residual DAGs
%   (affine/conv/normaffine/avgpool/relu/add). The generalisation of gpu_bab_crown_spec from
%   FC to conv: conv/BN/avgpool/add are LINEAR (exact interval forward + exact adjoint
%   backward), only ReLU is relaxed. The backward coefficient A is nSpec-by-dim-by-B (the SPEC
%   rows, nSpec ~ nClasses-1, NOT the layer width), so the memory is nSpec*HWC*B -- feasible
%   for conv, unlike batching the tight intermediate bounds' eye(nk) seed.
%
%   FULL DAG: the forward IBP caches each op's output bounds (so a residual 'add' fetches BOTH
%   its inputs and any op consuming a non-previous op bounds correctly) and the backward CROWN
%   routes each op's coefficient to op.src -- an 'add' routes it UNCHANGED to BOTH inputs. A
%   pure chain (every src==k-1, no 'add') reduces to the rolling forward + sequential backward
%   exactly (bound-for-bound unchanged). 'maxpool' is NOT handled here (per-node window
%   relaxation is not yet batched) and ERRORS -> sound-by-refusal (caller runs the tight path).
%
%   [margins, preL, preU] = GPU_BAB_CROWN_SPEC_DAG(ops, x_lb, x_ub, C, precision, fixings, rootBounds)
%     ops        : op list (nn_to_ops); affine/conv/normaffine/avgpool/relu/add (DAG via op.src,
%                  'add' via op.inputs=[a b]); NO 'maxpool'
%     x_lb,x_ub  : n-by-B input box columns (B = batch/node dim)
%     C          : nSpec-by-nOut spec
%     fixings    : optional cell(nOps,1) of per-relu dim_k-by-B node clamps (-1/0/+1)
%     rootBounds : optional struct with fields .preL,.preU (cell(nOps,1), dim_k-by-1 TIGHT
%                  pre-activation bounds computed ONCE at the root by gpu_bab_crown_tight).
%                  When given, the loose per-node IBP forward is SKIPPED and each node's
%                  pre-activation bounds are the (broadcast) root bounds clamped per node by
%                  the fixings -- tight bounds at batched speed (the #1 tightness lever).
%     margins    : nSpec-by-B lower bound on C*f over each node's clamped box
%     preL,preU  : per-relu (clamped) pre-activation bounds, dim_k-by-B
%
%   SOUNDNESS: interval-conv forward (Wp*lb+Wn*ub, the tightest linear interval) is sound; the
%   conv/avgpool/normaffine/add backward adjoints are EXACT (linear, no relaxation -- 'add'
%   routes the coefficient unchanged to both summands); the ReLU relaxation is the standard
%   sign-aware lower/upper line; the per-node clamps only tighten pre-activation bounds. With
%   rootBounds, the reused root bounds hold over the FULL input box -- a superset of every
%   node's sub-region -- and the per-neuron clamp (active: l>=0, inactive: u<=0) is the split's
%   domain restriction, so the bounds stay sound and are tighter than the per-node IBP. The DAG
%   backward is sound by induction: ops are topologically ordered, so every consumer of op k is
%   processed before op k, hence skipA{k} (the total coefficient on op k's output) is complete
%   when op k is bounded. For every x in node k's box, C*f(x) >= margins(:,k).

    if nargin < 5 || isempty(precision), precision = 'single'; end
    if nargin < 6, fixings = {}; end
    if nargin < 7, rootBounds = []; end
    if nargin < 8, alphaCell = {}; end
    if nargin < 9, vmag = {}; end
    if nargin < 10, cSel = {}; end
    B = size(x_lb, 2); nSpec = size(C, 1); nOps = numel(ops); n = size(x_lb, 1);
    % BINDING-SPEC: a per-node seed (nSpec_eff x nOut x B) overrides the full C as the backward
    % seed (nSpec_eff rows, typically 1). Only the score uses it; bound-invariant for branching.
    useCSel = ~isempty(cSel);
    if useCSel, nSpec = size(cSel, 1); end
    preL = cell(nOps,1); preU = cell(nOps,1);
    % SOUND-FP32 (M3b) state. doErr=true only on the single path WITH a (double) vmag majorant ->
    % the backward widens OUTWARD by derr (per-op FP32 roundoff, nSpec x B) + a final-contraction rad.
    % vmag absent / double -> doErr=false -> NO widening -> byte-identical to the prior screen.
    % SOUND-FP32 (M3b Step 3): doErr -> accumulate a sound bound `derr` (nSpec x B) on the FP32
    % backward roundoff (mirroring the validated gpu_bab_crown_tight derr, batched) + a final
    % contraction `rad`, and widen the margins OUTWARD (down) by (rad+derr) so
    % `margins_single <= margins_double <= true`. soundFP32=doErr is the EMIT gate the caller checks.
    % SOUNDNESS NB: the derr GEMMs are computed in DOUBLE (pagemtimes(double(|A|),vmag); vmag is the
    % double IBP majorant) -- pagemtimes(single,double) errors AND a single vmag is too coarse. Each
    % term is cast back to `precision` for the (single) derr accumulator (the i_gamma 2x inflation
    % absorbs that cast). vmag absent / double -> doErr=false -> NO widening -> byte-identical screen.
    % FAIL-CLOSED (adversarial review HOLE 3): the sound widening requires the ROOT-tight bounds,
    % which crown_tight already outward-widened (its preL/preU via i_backward+vmag). The
    % isempty(rootBounds) IBP-FORWARD branch below computes cl/cu in bare single WITHOUT widening,
    % so the relu relaxation it builds is NOT a sound over-approx -> never emit from it. Require
    % rootBounds for doErr; the emit path (cifar/tiny) always supplies them (rootTight).
    doErr = strcmp(precision, 'single') && ~isempty(vmag) && ~isempty(rootBounds);
    soundFP32 = doErr;
    us = single(eps('single') / 2);            % single unit roundoff (the d=d+... add-chain rounds here)
    if doErr
        derr = zeros(nSpec, B, precision);     % sound bound on the FP32 backward roundoff (accumulated)
        dmag = zeros(nSpec, B, precision);     % running sum of |d-terms| (the cross-op d=d+... chain)
    end

    if isempty(rootBounds)
        % ---- forward IBP (batched, FULL DAG), pre-activation bounds at each ReLU, with node
        % clamps. cl/cu cache each op's output bounds (index k+1 = op k; index 1 = op 0 = input)
        % so a residual 'add' fetches BOTH its inputs (out[a]+out[b], LINEAR -> exact) and any op
        % consuming a non-previous op (skip branch) bounds correctly. Sequential nets reduce to
        % the rolling bounds exactly. ----
        cl = cell(nOps+1,1); cu = cell(nOps+1,1);
        cl{1} = cast(x_lb, precision); cu{1} = cast(x_ub, precision);
        for k = 1:nOps
            op = ops{k};
            if strcmp(op.type, 'add')
                a = op.inputs(1)+1; b = op.inputs(2)+1;
                cl{k+1} = cl{a} + cl{b}; cu{k+1} = cu{a} + cu{b};
                continue;
            end
            if strcmp(op.type, 'concat')
                ic = op.inputs + 1;                       % stack input bounds (LINEAR -> exact)
                cl{k+1} = vertcat(cl{ic}); cu{k+1} = vertcat(cu{ic});
                continue;
            end
            if strcmp(op.type, 'product')                 % bilinear: store stacked input bounds + corner interval
                ia = op.inputs(1)+1; ib = op.inputs(2)+1;
                la = cl{ia}; ua = cu{ia}; lyy = cl{ib}; uyy = cu{ib};
                preL{k} = [la; lyy]; preU{k} = [ua; uyy];  % feed the McCormick backward
                c1=la.*lyy; c2=la.*uyy; c3=ua.*lyy; c4=ua.*uyy;
                cl{k+1} = min(min(c1,c2),min(c3,c4)); cu{k+1} = max(max(c1,c2),max(c3,c4));
                continue;
            end
            s = op.src + 1; lb = cl{s}; ub = cu{s};
            switch op.type
                case 'affine'
                    W = cast(op.W, precision); bb = cast(op.b(:), precision);
                    Wp = max(W,0); Wn = min(W,0);
                    cl{k+1} = Wp*lb + Wn*ub + bb; cu{k+1} = Wp*ub + Wn*lb + bb;
                case 'conv'
                    [cl{k+1}, cu{k+1}] = i_conv_ibp(op, lb, ub, precision);
                case 'normaffine'
                    sf = i_bcast_flat(op.scale, op.shape, precision);
                    tf = i_bcast_flat(op.shift, op.shape, precision);
                    pos = sf >= 0;
                    cl{k+1} = (sf.*lb).*pos + (sf.*ub).*(~pos) + tf;
                    cu{k+1} = (sf.*ub).*pos + (sf.*lb).*(~pos) + tf;
                case 'avgpool'
                    [cl{k+1}, cu{k+1}] = i_avgpool_ibp(op, lb, ub, precision);
                case 'maxpool'
                    preL{k} = lb; preU{k} = ub;            % window input bounds for the backward relaxation
                    [cl{k+1}, cu{k+1}] = i_maxpool_ibp(op, lb, ub, precision);
                case 'relu'
                    if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                        fx = fixings{k};
                        lb(fx == 1)  = max(lb(fx == 1),  0);
                        ub(fx == -1) = min(ub(fx == -1), 0);
                    end
                    preL{k} = lb; preU{k} = ub;
                    cl{k+1} = max(lb,0); cu{k+1} = max(ub,0);
                otherwise
                    error('gpu_bab_crown_spec_dag:op', ...
                        'Unsupported op "%s" (affine/conv/normaffine/avgpool/relu/add only).', op.type);
            end
        end
    else
        % ---- ROOT-TIGHT REUSE: skip the loose per-node IBP forward; build each ReLU's
        % pre-activation bounds from the TIGHT root bounds (broadcast to all B nodes), then
        % clamp per node by the fixings. The backward pass below only reads preL/preU at ReLU
        % ops, so no forward propagation is needed. SOUND (root bounds hold over the full box
        % >= each node's region; the own-neuron clamp is the split's domain restriction) and
        % much tighter than IBP. Non-relu ops (incl. 'add') carry no relu bounds; validate support.
        if ~isstruct(rootBounds) || ~all(isfield(rootBounds, {'preL','preU'}))
            error('gpu_bab_crown_spec_dag:rootBounds', ...
                'rootBounds must be a struct with fields preL and preU (per-op cell arrays).');
        end
        tmpl = cast(x_lb(1), precision); % scalar template for device+precision of the broadcasts (no n-by-B copy)
        for k = 1:nOps
            tk = ops{k}.type;
            if strcmp(tk, 'relu')
                l = repmat(cast(rootBounds.preL{k}, 'like', tmpl), 1, B);
                u = repmat(cast(rootBounds.preU{k}, 'like', tmpl), 1, B);
                if ~isempty(fixings) && numel(fixings) >= k && ~isempty(fixings{k})
                    fx = fixings{k};
                    l(fx == 1)  = max(l(fx == 1),  0);
                    u(fx == -1) = min(u(fx == -1), 0);
                end
                preL{k} = l; preU{k} = u;
            elseif strcmp(tk, 'product')
                % bilinear: reuse the root stacked input bounds (broadcast). No per-node clamp --
                % products carry no ReLU fixings; root bounds hold over the full box >= node region
                % (sound, possibly looser than re-concretizing). Backward reads preL{k}=[xb;yb].
                preL{k} = repmat(cast(rootBounds.preL{k}, 'like', tmpl), 1, B);
                preU{k} = repmat(cast(rootBounds.preU{k}, 'like', tmpl), 1, B);
            elseif strcmp(tk, 'maxpool')
                % maxpool: reuse the root window input bounds (broadcast) for the backward relaxation.
                % No per-node ReLU clamp (maxpool carries no fixings); root bounds hold over the full
                % box >= the node region (sound, possibly looser). Backward reads preL{k}/preU{k}.
                preL{k} = repmat(cast(rootBounds.preL{k}, 'like', tmpl), 1, B);
                preU{k} = repmat(cast(rootBounds.preU{k}, 'like', tmpl), 1, B);
            elseif ~any(strcmp(tk, {'affine','conv','normaffine','avgpool','add','concat'}))
                error('gpu_bab_crown_spec_dag:op', ...
                    'Unsupported op "%s" (affine/conv/normaffine/avgpool/relu/maxpool/add/concat/product only).', tk);
            end
        end
    end

    % ---- backward CROWN (batched over B, FULL DAG): lower bound on C*f ----
    % skipA{k} = accumulated nSpec-by-width_k-by-B coefficient on op k's OUTPUT. A single-input
    % op routes its coefficient (after its linear/relax backward) to op.src; an 'add' routes it
    % UNCHANGED to BOTH inputs (LINEAR -> exact). inputSkipA = total coefficient on op 0 (the
    % engine input). Topological order => every consumer of op k is visited before op k, so
    % skipA{k} is complete when op k is processed (sound by induction). For a pure chain this is
    % bound-for-bound the rolling pass (skipA{k-1} == the old running A).
    skipA = cell(nOps, 1);
    if useCSel
        skipA{nOps} = cast(cSel, precision);             % per-node binding-spec seed (nSpec_eff x nOut x B)
    else
        skipA{nOps} = repmat(cast(C, precision), 1, 1, B);   % spec sits on the network output (op nOps)
    end
    d = zeros(nSpec, B, precision);
    inputSkipA = [];
    % BaBSR sensitivity score (computed only if the caller asks for the 4th output). At each ReLU,
    % neuron i's relaxation LOWERS the bound by exactly Aneg(s,i,b)*bu(i,b) for spec s; the total
    % slack a split there can recover is sum_s |Aneg(s,i,b)|*bu(i,b). Branching on the largest such
    % score (output-sensitivity x gap) closes the bound in far fewer nodes than largest-gap alone.
    wantScore = nargout >= 4;
    scoreCell = cell(nOps, 1);
    for k = nOps:-1:1
        A = skipA{k};
        if isempty(A), continue; end
        skipA{k} = [];                                   % free: keep only live branches resident
        op = ops{k};
        if doErr && ~isempty(getenv('NNV_DEBUG_DERR'))   % per-op derr trace (cumulative BEFORE op k)
            fprintf('[derrtrace] before op %2d %-10s cumderr_min=%.5g\n', k, op.type, gather(min(derr(:))));
        end
        if strcmp(op.type, 'add')
            for ii = 1:numel(op.inputs)
                s = op.inputs(ii);
                if doErr && (s == 0 || ~isempty(skipA{s}))   % a MERGE (accumulating +) -> coeff add-chain rounding
                    if s == 0, vm = double(vmag{1}); else, vm = double(vmag{s + 1}); end
                    derr = derr + i_dterm_sd(A, reshape(vm, [], 1, 1), nOps, nSpec, B);
                end
                if s == 0
                    if isempty(inputSkipA), inputSkipA = A; else, inputSkipA = inputSkipA + A; end
                elseif isempty(skipA{s}), skipA{s} = A;
                else, skipA{s} = skipA{s} + A;
                end
            end
            continue;
        end
        if strcmp(op.type, 'concat')
            % out = [in_1; in_2; ...] -> SLICE the coefficient's dim axis (axis 2) back to each
            % input's block (transpose of stacking, LINEAR -> exact). A is nSpec x outDim x B.
            off = 0;
            for ii = 1:numel(op.inputs)
                sz = op.sizes(ii); s = op.inputs(ii);
                Ablk = A(:, off+(1:sz), :); off = off + sz;
                if doErr && ((s == 0 && ~isempty(inputSkipA)) || (s ~= 0 && ~isempty(skipA{s})))  % MERGE (HOLE 1)
                    if s == 0, vm = double(vmag{1}); else, vm = double(vmag{s + 1}); end
                    derr = derr + i_dterm_sd(Ablk, reshape(vm, [], 1, 1), nOps, nSpec, B);
                end
                if s == 0
                    if isempty(inputSkipA), inputSkipA = Ablk; else, inputSkipA = inputSkipA + Ablk; end
                elseif isempty(skipA{s}), skipA{s} = Ablk;
                else, skipA{s} = skipA{s} + Ablk;
                end
            end
            continue;
        end
        if strcmp(op.type, 'product')
            % out = in[a].*in[b], bilinear: McCormick (gpu_bab_mul_relax), sign-aware, route a
            % DISTINCT coefficient to each input (like 'add' but Ax/Ay differ). A is nSpec x wa x B;
            % preL{k}=[xb;yb] is 2*wa x B (stacked input bounds). Sound: +coeff uses the LOWER
            % plane for a spec lower bound; -coeff the UPPER plane (the al/au rule).
            wa = op.sizes(1); lz = preL{k}; uz = preU{k};
            la = lz(1:wa,:); ua = uz(1:wa,:); lyy = lz(wa+1:end,:); uyy = uz(wa+1:end,:);
            [aL,bL,cL,aU,bU,cU] = gpu_bab_mul_relax(la, ua, lyy, uyy, [], [], precision);  % each wa x B
            if doErr   % McCormick slope-update error scales as the PRODUCT va*vb (the slopes are
                       % UNBOUNDED, unlike relu's [0,1]) + the cL/cU intercept (per-node) + its d-add
                va = double(vmag{op.inputs(1)+1}); vb = double(vmag{op.inputs(2)+1}); vab = va .* vb;  % wa x 1
                ccl = double(abs(cL) + abs(cU));                                          % wa x B (per-node)
                derr = derr + i_dterm_sd(A, reshape(2*vab, wa, 1, 1), 2,  nSpec, B) ...
                            + i_dterm_sd(A, reshape(ccl,   wa, 1, B), wa, nSpec, B);
                dmag = dmag + single(reshape(pagemtimes(double(abs(A)), reshape(ccl, wa, 1, B)), nSpec, B));
                derr = derr + us * dmag;
            end
            % spec_dag computes a LOWER bound on C*f (like its ReLU relax): +coeff -> LOWER plane,
            % -coeff -> UPPER plane.
            Apos = max(A,0); Aneg = min(A,0);
            aLr = reshape(aL,1,wa,B); aUr = reshape(aU,1,wa,B);
            bLr = reshape(bL,1,wa,B); bUr = reshape(bU,1,wa,B);
            Ax = Apos.*aLr + Aneg.*aUr; Ay = Apos.*bLr + Aneg.*bUr;
            d = d + reshape(pagemtimes(Apos, reshape(cL,wa,1,B)), nSpec, B) ...
                  + reshape(pagemtimes(Aneg, reshape(cU,wa,1,B)), nSpec, B);
            sa = op.inputs(1);
            if doErr && ((sa == 0 && ~isempty(inputSkipA)) || (sa ~= 0 && ~isempty(skipA{sa})))  % MERGE (HOLE 2)
                if sa == 0, vm = double(vmag{1}); else, vm = double(vmag{sa + 1}); end
                derr = derr + i_dterm_sd(Ax, reshape(vm, [], 1, 1), nOps, nSpec, B);
            end
            if sa == 0
                if isempty(inputSkipA), inputSkipA = Ax; else, inputSkipA = inputSkipA + Ax; end
            elseif isempty(skipA{sa}), skipA{sa} = Ax;
            else, skipA{sa} = skipA{sa} + Ax;
            end
            sb = op.inputs(2);
            if doErr && ((sb == 0 && ~isempty(inputSkipA)) || (sb ~= 0 && ~isempty(skipA{sb})))  % MERGE (HOLE 2)
                if sb == 0, vm = double(vmag{1}); else, vm = double(vmag{sb + 1}); end
                derr = derr + i_dterm_sd(Ay, reshape(vm, [], 1, 1), nOps, nSpec, B);
            end
            if sb == 0
                if isempty(inputSkipA), inputSkipA = Ay; else, inputSkipA = inputSkipA + Ay; end
            elseif isempty(skipA{sb}), skipA{sb} = Ay;
            else, skipA{sb} = skipA{sb} + Ay;
            end
            continue;
        end
        if doErr, vin = double(vmag{op.src + 1}); end     % input value-magnitude majorant (op 0 = input)
        switch op.type
            case 'affine'
                W = cast(op.W, precision); bb = cast(op.b(:), precision);
                if doErr   % A*W contracts over out-width = size(A,2); output value-mag (no cancel) = |W|*vin+|b|
                    omg = double(abs(W)) * double(vin) + double(abs(bb));            % out-width x 1
                    derr = derr + i_dterm_sd(A, reshape(omg, [], 1, 1), size(A,2), nSpec, B);
                    dmag = dmag + i_dterm_raw_sd(A, reshape(double(abs(bb)), [], 1, 1), nSpec, B);  % d=d+A*b chain
                    derr = derr + us * dmag;
                end
                d = d + reshape(pagemtimes(A, bb), nSpec, B);
                A = pagemtimes(A, W);
            case 'conv'
                if doErr   % magnitude adjoint transpconv(|A|,|W|) [len kh*kw*outCh] + bias [len prod(outShape)]
                    opm = op; opm.W = abs(op.W); opm.b = abs(op.b);
                    [Amag, dmc] = i_conv_backward(abs(A), zeros(nSpec, B, precision), opm, precision); % nSpec x inDim x B, nSpec x B
                    mC = size(op.W,1) * size(op.W,2) * size(op.W,4); nO = prod(op.outShape);
                    % Amag/dmc are SINGLE-precision reductions (lengths ~mC / nO), so they can round DOWN
                    % by up to gamma_raw -- the only derr operands NOT computed in double here. The
                    % the (1+2^-7) gamma inflation covers the cast + accumulation, NOT this operand
                    % shortfall (a single reduction can round down by gamma_nO), so inflate it explicitly.
                    % Inflate each by (1+2*gamma_raw) >= 1/(1-gamma_raw) -> a sound over-estimate of the
                    % true magnitude (adversarial-review finding 2026-06-18; ~0.8% for nO=65536, so the
                    % derr stays ~halved while remaining a rigorous bound).
                    Amag = Amag .* (1 + 2*i_gamma_raw_sd(mC, precision));
                    dmc  = dmc  .* (1 + 2*i_gamma_raw_sd(nO, precision));
                    derr = derr + i_dterm_sd(Amag, reshape(double(vin), [], 1, 1), mC, nSpec, B) ...
                                + i_gamma_sd(nO, precision) * dmc;                    % bias reduction length
                    dmag = dmag + dmc; derr = derr + us * dmag;                       % d=d+bias add-chain
                end
                [A, d] = i_conv_backward(A, d, op, precision);
            case 'normaffine'
                sf = i_bcast_flat(op.scale, op.shape, precision);
                tf = i_bcast_flat(op.shift, op.shape, precision);
                if doErr   % TWO FP32 ops: slope multiply A.*sf' AND intercept matmul d=d+A*tf
                    sfm = double(i_bcast_flat(abs(op.scale), op.shape, precision));   % |slope| flat (dim x 1)
                    tfm = double(i_bcast_flat(abs(op.shift), op.shape, precision));   % |shift| flat (dim x 1)
                    Asf = abs(A) .* reshape(single(sfm), 1, [], 1);                   % |A| .* |sf|' (nSpec x dim x B)
                    derr = derr + i_dterm_sd(Asf, reshape(double(vin), [], 1, 1), 2,          nSpec, B) ...
                                + i_dterm_sd(A,   reshape(tfm,        [], 1, 1), numel(tfm),  nSpec, B);
                    dmag = dmag + i_dterm_raw_sd(A, reshape(tfm, [], 1, 1), nSpec, B);% d=d+A*tf add-chain
                    derr = derr + us * dmag;
                end
                d = d + reshape(pagemtimes(A, tf), nSpec, B);
                A = A .* reshape(sf, 1, [], 1);
            case 'avgpool'
                [A, d] = i_avgpool_backward(A, d, op, precision);
                if doErr   % A is now on the INPUT; no bias -> no d-add
                    derr = derr + i_dterm_sd(A, reshape(double(vin), [], 1, 1), prod(op.pool), nSpec, B);
                end
            case 'maxpool'
                % Sound CROWN LOWER relaxation through maxpool (ported batched from gpu_bab_crown_tight
                % i_maxpool_backward; spec_dag is lower-direction only). soundFP32 emit is NOT supported
                % for maxpool -> refuse (the trust-FP32 screen runs soundFP32=false, so it is unaffected).
                if doErr
                    error('gpu_bab_crown_spec_dag:maxpoolSoundFP32', ...
                        'maxpool sound-FP32 emit not supported -- refused (use the FP64 confirm).');
                end
                [A, d] = i_maxpool_backward_batched(A, d, op, preL{k}, preU{k}, precision);
            case 'relu'
                l = preL{k}; u = preU{k};                  % dim x B
                ak = []; if ~isempty(alphaCell) && numel(alphaCell) >= k, ak = alphaCell{k}; end
                [au, bu, al] = i_relu_relax(l, u, precision, ak);
                dim = size(l, 1);
                Apos = max(A, 0); Aneg = min(A, 0);
                if doErr   % slope mult (|A_in|<=|A|, no cancel) + the bu intercept (PER-NODE dim x B, R1a)
                           % |A| covers the d-pass; +4u for the au=u/(u-l)/bu=-au*l derivation rounding
                    Abu  = i_dterm_raw_sd(A, reshape(double(abs(bu)), dim, 1, B), nSpec, B);   % |A|*|bu| (nSpec x B)
                    derr = derr + i_dterm_sd(A, reshape(double(vin), dim, 1, 1), 8, nSpec, B) ...    % m=8: slope multiply (~u) + the relu RELAXATION-DERIVATION rounding (~3u*vin at z=u; scales with vin, NOT bu) -- review 2026-06-18
                                + single(i_gamma_sd(dim,'double') + double(4)*double(us)) * Abu;
                    dmag = dmag + Abu; derr = derr + us * dmag;        % d=d+Aneg*bu add-chain
                end
                d = d + reshape(pagemtimes(Aneg, reshape(bu, dim, 1, B)), nSpec, B);
                if wantScore
                    % sum_s |Aneg(s,i,b)| * bu(i,b) -> dim x B (the BaBSR split score for this layer)
                    scoreCell{k} = reshape(sum(abs(Aneg), 1), dim, B) .* bu;
                end
                A = Apos .* reshape(al, 1, dim, B) + Aneg .* reshape(au, 1, dim, B);
            otherwise
                error('gpu_bab_crown_spec_dag:op', ...
                    'Unsupported op "%s" in backward (affine/conv/normaffine/avgpool/relu/add only).', op.type);
        end
        s = op.src;                                       % route coefficient to op.src
        if doErr && (s == 0 || ~isempty(skipA{s}))        % a MERGE (accumulating +) -> coeff add-chain rounding
            if s == 0, vm = double(vmag{1}); else, vm = double(vmag{s + 1}); end
            derr = derr + i_dterm_sd(A, reshape(vm, [], 1, 1), nOps, nSpec, B);
        end
        if s == 0
            if isempty(inputSkipA), inputSkipA = A; else, inputSkipA = inputSkipA + A; end
        elseif isempty(skipA{s}), skipA{s} = A;
        else, skipA{s} = skipA{s} + A;
        end
    end

    if isempty(inputSkipA)
        if doErr, derr = derr + us * abs(d); margins = d - derr;   % constant output: only the d add-chain rounds
        else,     margins = d; end                                 % output independent of input (degenerate)
        return;
    end
    A = inputSkipA;
    Apos = max(A, 0); Aneg = min(A, 0);
    lbcol = reshape(cast(x_lb, precision), n, 1, B);
    ubcol = reshape(cast(x_ub, precision), n, 1, B);
    aposx = reshape(pagemtimes(Apos, lbcol), nSpec, B);
    anegx = reshape(pagemtimes(Aneg, ubcol), nSpec, B);
    margins = aposx + anegx + d;
    if doErr   % widen the LOWER bound OUTWARD (down) by the final-contraction rad + the accumulated
               % per-op derr, so margins_single <= margins_double <= true (monotone-sound, can EMIT)
        rad  = i_outward_rad_sd(A, x_lb, x_ub, precision);
        derr = derr + us * abs(d) + cast(3,precision)*us * (abs(aposx) + abs(anegx));  % '+d' u*|d| + the two final adds fl(fl(aposx+anegx)+d) ~2u*P (3*us*P w/ margin; reduced rad k no longer covers them) -- re-review CHECK 3, 2026-06-18
        margins = margins - rad - derr;
        if ~isempty(getenv('NNV_DEBUG_DERR'))             % M3b tightening diagnosis (no behaviour change)
            mraw = margins + rad + derr;                  % the un-widened margin
            [mw, iw] = min(margins(:));                   % binding spec (widened)
            fprintf('[derr] raw_min=%.6g widened_min=%.6g | binding: raw=%.6g rad=%.6g derr=%.6g\n', ...
                gather(min(mraw(:))), gather(mw), gather(mraw(iw)), gather(rad(iw)), gather(derr(iw)));
            vm = cellfun(@(v) max(abs(double(v(:)))), vmag);
            fprintf('[derr] vmag(IBP) max-per-op: %s\n', mat2str(gather(vm(:))', 4));
        end
    end
end

% ---------------------------------------------------------------------------------------
function [olb, oub] = i_conv_ibp(op, lb, ub, precision)
% Interval conv (batched over B): the Wp/Wn split of the affine IBP with the matmul replaced
% by dlconv. Sound (tightest interval for the linear conv map). lb/ub are prod(inShape)-by-B.
    ish = op.inShape; osh = op.outShape; B = size(lb,2);
    W = cast(op.W, precision); Wp = max(W,0); Wn = min(W,0);
    bb = reshape(cast(op.b(:), precision), [1 1 osh(3)]);
    L4 = dlarray(reshape(lb, [ish(1) ish(2) ish(3) B]), 'SSCB');
    U4 = dlarray(reshape(ub, [ish(1) ish(2) ish(3) B]), 'SSCB');
    pad2 = [op.pad(1) op.pad(3); op.pad(2) op.pad(4)];      % [t l; b r]
    args = {'Stride', op.stride, 'Padding', pad2, 'DilationFactor', op.dil};
    Lo = dlconv(L4, Wp, bb, args{:}) + dlconv(U4, Wn, 0, args{:});
    Hi = dlconv(U4, Wp, bb, args{:}) + dlconv(L4, Wn, 0, args{:});
    olb = reshape(extractdata(Lo), [prod(osh) B]);
    oub = reshape(extractdata(Hi), [prod(osh) B]);
end

function [A2, d2] = i_conv_backward(A, d, op, precision)
% Exact CROWN backward through a conv (linear), batched over B: fold B into the
% dltranspconv batch dim (nSpec*B). A: nSpec x prod(outShape) x B.
    nSpec = size(A,1); B = size(A,3);
    osh = op.outShape; ish = op.inShape; W = cast(op.W, precision);
    Aperm = permute(A, [2 1 3]);                            % dim x nSpec x B
    A4 = dlarray(reshape(Aperm, [osh(1) osh(2) osh(3) nSpec*B]), 'SSCB');
    Afull = extractdata(dltranspconv(A4, W, 0, 'Stride', op.stride, 'Cropping', 0, 'DilationFactor', op.dil));
    pt = op.pad(1); pl = op.pad(3);
    Ain = zeros([ish(1) ish(2) ish(3) nSpec*B], 'like', Afull);
    hi = min(ish(1), size(Afull,1)-pt); wi = min(ish(2), size(Afull,2)-pl);
    if hi>0 && wi>0, Ain(1:hi,1:wi,:,:) = Afull(pt+(1:hi), pl+(1:wi), :, :); end
    A2 = reshape(Ain, [prod(ish) nSpec B]);                 % inDim x nSpec x B
    A2 = permute(A2, [2 1 3]);                              % nSpec x inDim x B
    bc = reshape(cast(op.b(:), precision), [1 1 osh(3)]);
    A4d = reshape(extractdata(A4), [osh(1) osh(2) osh(3) nSpec B]);
    dinc = reshape(sum(A4d .* bc, [1 2 3]), [nSpec B]);
    d2 = d + dinc;
end

function [olb, oub] = i_avgpool_ibp(op, lb, ub, precision)
% Interval avgpool (batched): mean over each non-overlapping window -> monotone, so the
% exact interval is the avgpool of the bounds. lb/ub prod(inShape)-by-B.
    ish = op.inShape; osh = op.outShape; B = size(lb,2);
    L4 = reshape(cast(lb,precision), [ish(1) ish(2) ish(3) B]);
    U4 = reshape(cast(ub,precision), [ish(1) ish(2) ish(3) B]);
    olb = reshape(i_pool_mean(L4, op), [prod(osh) B]);
    oub = reshape(i_pool_mean(U4, op), [prod(osh) B]);
end

function Y = i_pool_mean(X, op)
% Non-overlapping mean pool (stride==pool, unpadded) of X [H W C B].
    osh = op.outShape; kh = op.pool(1); kw = op.pool(2); B = size(X,4);
    Y = zeros([osh(1) osh(2) osh(3) B], 'like', X);
    for oh = 1:osh(1)
        for ow = 1:osh(2)
            rh = (oh-1)*op.stride(1) + (1:kh); rw = (ow-1)*op.stride(2) + (1:kw);
            Y(oh,ow,:,:) = mean(mean(X(rh,rw,:,:),1),2);
        end
    end
end

function [A2, d2] = i_avgpool_backward(A, d, op, precision)
% Exact CROWN backward through non-overlapping avgpool (batched over B): distribute
% A_out/(kh*kw) uniformly to each window's input cells. A: nSpec x prod(outShape) x B.
    nSpec = size(A,1); B = size(A,3); osh = op.outShape; ish = op.inShape;
    kh = op.pool(1); kw = op.pool(2);
    Aperm = permute(A, [2 1 3]);                            % dim x nSpec x B
    A4 = reshape(cast(Aperm,precision), [osh(1) osh(2) osh(3) nSpec*B]);
    Aup = repelem(A4, kh, kw, 1, 1) / (kh*kw);
    Ain = zeros([ish(1) ish(2) ish(3) nSpec*B], precision);
    hi = min(ish(1), size(Aup,1)); wi = min(ish(2), size(Aup,2));
    Ain(1:hi,1:wi,:,:) = Aup(1:hi,1:wi,:,:);
    A2 = permute(reshape(Ain, [prod(ish) nSpec B]), [2 1 3]);
    d2 = d;
end

function [olb, oub] = i_maxpool_ibp(op, lb, ub, precision)
% Interval maxpool (batched): max over each window is monotone, so the exact interval is the
% maxpool of the bounds. lb/ub prod(inShape)-by-B.
    ish = op.inShape; osh = op.outShape; B = size(lb,2);
    L4 = reshape(cast(lb,precision), [ish(1) ish(2) ish(3) B]);
    U4 = reshape(cast(ub,precision), [ish(1) ish(2) ish(3) B]);
    olb = reshape(i_pool_max(L4, op), [prod(osh) B]);
    oub = reshape(i_pool_max(U4, op), [prod(osh) B]);
end

function Y = i_pool_max(X, op)
% Max pool of X [H W C B] (stride/overlap allowed, unpadded -- enforced at extraction).
    osh = op.outShape; kh = op.pool(1); kw = op.pool(2); B = size(X,4);
    Y = -inf([osh(1) osh(2) osh(3) B], 'like', X);
    for oh = 1:osh(1)
        for ow = 1:osh(2)
            rh = (oh-1)*op.stride(1) + (1:kh); rw = (ow-1)*op.stride(2) + (1:kw);
            Y(oh,ow,:,:) = max(max(X(rh,rw,:,:),[],1),[],2);
        end
    end
end

function [A2, d] = i_maxpool_backward_batched(A, d, op, l, u, precision)
% Batched (over the B node columns) sound CROWN LOWER backward through maxpool, ported from the
% serial gpu_bab_crown_tight i_maxpool_backward (LOWER direction only -- spec_dag computes lower
% margins). For each output cell o = max over its window:
%   lower  y_o >= x_m       (m = argmax of the window's LOWER bounds; max >= any element)
%   upper  y_o <= x_m       if m DOMINATES (l_m >= u_i for every other window i) -- exact
%          y_o <= max_i u_i otherwise (a constant; sound). The constant upper is robust-by-default
%          (an i_maxpool_relax 'decided' miss only loosens -- never unsound).
% A: nSpec x n_out x B ; l,u: n_in x B ; d: nSpec x B. The window relaxation differs per node
% column b, so build per-b selection and scatter (overlapping windows accumulate in the matmul).
% PERF (known, sound): this gathers to host + forms CPU sparse [n_out x n_in] per b, so for large
% HWC maxpools / large B it dominates and negates batching. The window->flat index maps are static
% per op, so a GPU-resident scatter (precompute maps once, argmax/umax on-device, no dense matrices)
% is the planned speedup -- tracked as future work; the current path is correct, just CPU-bound.
    nSpec = size(A,1); B = size(A,3);
    n_in = prod(op.inShape); n_out = prod(op.outShape);
    A2 = zeros(nSpec, n_in, B, 'like', A);
    for b = 1:B
        [mIdx, decided, umax] = i_maxpool_relax(op, gather(l(:,b)), gather(u(:,b)));
        Lmat  = sparse(1:n_out, mIdx, 1, n_out, n_in);             % y_o >= x_{m_o}
        dec   = find(decided);
        Umat  = sparse(dec, mIdx(dec), 1, n_out, n_in);            % decided: y_o <= x_{m_o}
        Ubias = zeros(n_out, 1); Ubias(~decided) = umax(~decided); % undecided: y_o <= umax
        Ab = double(gather(A(:,:,b))); Apos = max(Ab,0); Aneg = min(Ab,0);
        A2(:,:,b) = cast(Apos*Lmat + Aneg*Umat, precision);
        d(:,b)    = d(:,b) + cast(Aneg*Ubias, precision);
    end
end

function [mIdx, decided, umax] = i_maxpool_relax(op, l, u)
% Per output cell: m = flat input index of the window's argmax-LOWER element; decided = that
% element dominates (its lower bound >= every other window upper => max is exactly it); umax =
% max window upper. Unpadded windows (enforced at extraction); overlap allowed. (Verbatim port
% of gpu_bab_crown_tight i_maxpool_relax.)
    ish = op.inShape; osh = op.outShape;
    H = ish(1); W = ish(2); C = ish(3); Ho = osh(1); Wo = osh(2);
    kh = op.pool(1); kw = op.pool(2); sh = op.stride(1); sw = op.stride(2);
    L = reshape(double(l), [H W C]); U = reshape(double(u), [H W C]);
    n_out = prod(osh);
    mIdx = ones(n_out, 1); decided = false(n_out, 1); umax = zeros(n_out, 1);
    for c = 1:C
        for ow = 1:Wo
            rw = (ow-1)*sw + (1:kw);
            for oh = 1:Ho
                rh = (oh-1)*sh + (1:kh);
                lw = L(rh, rw, c); uw = U(rh, rw, c);
                [maxl, p] = max(lw(:));
                [pr, pc] = ind2sub([kh kw], p);
                o = sub2ind([Ho Wo C], oh, ow, c);
                mIdx(o) = sub2ind([H W C], rh(pr), rw(pc), c);
                uw2 = uw; uw2(p) = -inf;
                decided(o) = maxl >= max(uw2(:));               % argmax-lower dominates others
                umax(o) = max(uw(:));
            end
        end
    end
end

function v = i_bcast_flat(x, sh, precision)
% Broadcast x ([1 1 C] / [H W C] / scalar) to a flat [prod(sh) x 1] column (column-major).
    v = reshape(zeros([sh(1) sh(2) sh(3)], precision) + cast(x, precision), [], 1);
end

function [au, bu, al] = i_relu_relax(l, u, precision, alpha)
% Per-neuron ReLU relaxation over [l,u] (elementwise, dim x B): stable-on (l>=0) identity,
% stable-off (u<=0) zero, unstable upper line au*z+bu / lower line al*z. al = min-area {0,1}
% by default; if alpha (dim x 1, the AMORTIZED root slopes) is given, al(unstable) = alpha
% (clamped [0,1] -> sound), broadcast over the B node columns.
    if nargin < 4, alpha = []; end
    au = zeros(size(l), precision); bu = zeros(size(l), precision); al = zeros(size(l), precision);
    act = (l >= 0); au(act) = 1; al(act) = 1;
    unst = (l < 0) & (u > 0); dn = u(unst) - l(unst);
    au(unst) = u(unst) ./ dn; bu(unst) = -au(unst) .* l(unst);
    if isempty(alpha)
        al(unst) = cast(u(unst) >= -l(unst), precision);     % min-area binary {0,1}
    else
        ab = repmat(cast(alpha(:), precision), 1, size(l,2));  % dim x 1 -> dim x B (root slopes)
        al(unst) = min(max(ab(unst), 0), 1);                   % clamp [0,1] (sound)
    end
end

% ===== SOUND-FP32 (M3b) helpers. i_gamma_sd / i_outward_rad_sd MUST stay identical to =========
% gpu_bab_crown_tight.m i_gamma / i_outward_rad (M0 will factor to one shared file). ===========
function t = i_dterm_sd(A, OP, m, nSpec, B)
% One sound-FP32 derr term for the batched backward: gamma(m) * (|A| * OP), where the contraction
% is done in DOUBLE (pagemtimes(single,double) ERRORS; a single vmag is too coarse), then cast back
% to single for the (single) derr accumulator (the 2x inflation in i_gamma_sd absorbs that cast).
%   A  : nSpec x w x B coefficient (single)
%   OP : the DOUBLE value-magnitude operand, w x 1 x 1 (shared, broadcast over B) OR w x 1 x B (per-node)
%   m  : the FP contraction length for this op (deterministic Higham gamma_m)
    t = single(i_gamma_sd(m, 'double') * reshape(pagemtimes(double(abs(A)), OP), nSpec, B));
end

function t = i_dterm_raw_sd(A, OP, nSpec, B)
% The raw magnitude |A| * OP (no gamma), for the dmag add-chain accumulator. Same DOUBLE GEMM as
% i_dterm_sd. A: nSpec x w x B (single); OP: double, w x 1 x 1 (shared) or w x 1 x B (per-node).
    t = single(reshape(pagemtimes(double(abs(A)), OP), nSpec, B));
end

function g = i_gamma_sd(m, precision)
% Deterministic Higham factor gamma_m = m*u/(1-m*u) for a length-m FP32 contraction; fail-closed to
% realmax when m*u>=0.5. INFLATION (M3b-T): spec_dag computes each derr term in DOUBLE
% (i_dterm_sd: single(gamma * pagemtimes(double|A|, double vmag))), so the only single-precision
% rounding in the derr value is the final single() cast (<= u) PLUS the single derr-accumulator sum
% (round-down ~ gamma_N for N derr-adds). A (1 + 2^-7) inflation covers BOTH the cast AND the
% accumulation for N <= 2^17 derr-adds (any net), provably -> a sound over-estimate of the Higham
% bound. (crown_tight.m's i_gamma now matches at (1+2^-7) -- same double-GEMM derr. The old 2x made
% derr ~2x larger than necessary, blocking cifar emit; adversarial review 2026-06-18 confirmed the
% relu-slope (m=8), final-add (3*us*|a|), and accumulator coverages this 2^-7 + the conv inflation.)
    u = single(eps('single') / 2);
    m = single(m);
    den = 1 - m * u;
    if den <= single(0.5)
        g = cast(realmax('single'), precision); return;
    end
    g = cast(1 + 2^-7, precision) * (m * u) / den;   % (1+2^-7): covers the final cast AND the single derr-accumulator sum for N<=2^17 derr-adds (any net). Was (1+2^-10) which covered only N<=16384 -- adversarial review 2026-06-18.
end

function g = i_gamma_raw_sd(m, precision)
% UN-inflated Higham factor gamma_m = m*u/(1-m*u): the relative down-rounding bound of a length-m
% SINGLE reduction. Used to over-estimate single-computed derr OPERANDS (conv Amag/dmc) via the
% factor (1+2*gamma_raw). That factor is a sound over-estimate -- (1+2g) >= 1/(1-g) <=> g(1-2g) >= 0
% -- ONLY for g <= 0.5, i.e. m*u <= 1/3 (den >= 2/3). Beyond g=0.5 the (1+2g) inflation UNDER-covers
% 1/(1-g), so FAIL-CLOSE to realmax (operand -> +Inf -> derr -> +Inf -> margin -Inf -> never emits)
% at den <= 2/3. This is UNREACHABLE for real convolutions (needs reduction length m > ~5.6e6; CIFAR
% nO=prod(outShape)=65536 gives gamma_raw ~ 4e-3), but it makes the over-estimate a PROVABLE invariant
% rather than a reachability argument (adversarial re-review 2026-06-18).
    u = single(eps('single') / 2);
    m = single(m);
    den = 1 - m * u;
    if den <= single(2/3)
        g = cast(realmax('single'), precision); return;
    end
    g = cast((m * u) / den, precision);
end

function rad = i_outward_rad_sd(A, x_lb, x_ub, precision)
% Batched final-contraction outward radius (mirrors gpu_bab_crown_tight i_outward_rad, computed in
% DOUBLE). A: nSpec x n x B (input-space coeff), x_lb/x_ub: n x B -> rad: nSpec x B (single).
    nSpec = size(A,1); B = size(A,3); n = size(x_lb,1);
    if ~strcmp(precision, 'single'), rad = zeros(nSpec, B, 'like', A); return; end
    u = single(eps('single') / 2); k = single(1 + 2^-7);    % M3b-T: double GEMM -> covers cast + accumulation (was 2x); final-add u*|a| covered at the margin block -- review 2026-06-18
    den = 1 - single(n) * u;
    if den <= single(0.5)
        rad = realmax('single') * ones(nSpec, B, 'like', A); return;
    end
    gbar = double(k) * (double(n) * double(u)) / double(den);
    xmag = max(abs(double(x_lb)), abs(double(x_ub)));            % n x B (double)
    rad = single(gbar * reshape(pagemtimes(double(abs(A)), reshape(xmag, n, 1, B)), nSpec, B) ...
                 + double(n) * realmin('single'));
end
