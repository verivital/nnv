function [verdict, info] = gpu_bab_halfspace_verify(net, lb, ub, prop, opts)
% GPU_BAB_HALFSPACE_VERIFY  Sound, ADDITIVE GPU-BaB pre-check for GENERAL-halfspace safety
%   specs (the control-style benchmarks: cersyve, lsnc_relu, linearizenn). Generalises the
%   argmax pre-check (gpu_bab_try_verify) to an arbitrary vnnlib unsafe region.
%
%   prop.Hg is an array of HalfSpace objects; each HalfSpace {y : G*y <= g} is ONE unsafe
%   DISJUNCT (polytope). The vnnlib unsafe region is their UNION. The property HOLDS (unsat /
%   safe) iff the reachable output avoids EVERY disjunct. A disjunct {G*y <= g} is AVOIDED iff
%   SOME row i has G_i*y - g_i > 0 for the WHOLE reachable set (R lies entirely outside row i's
%   halfspace -> outside the polytope). So:
%       'robust' (safe)  <=>  for every disjunct d,  max_i ( lb(G_d,i * y) - g_d,i ) > guardTol
%   where lb(.) is a sound CROWN lower bound over the input box. (argmax is the special case:
%   each disjunct is a single row e_target - e_j, and "all disjuncts avoided" == "all margins>0".)
%
%   verdict : 'robust'  -- PROVEN safe (every disjunct provably avoided)  -> caller emits unsat
%             'unknown' -- root CROWN did not avoid some disjunct          -> caller runs Star
%             'skip'    -- unsupported net / orientation guard failed       -> caller runs Star
%
%   SOUNDNESS (never a wrong unsat / -150):
%     (1) nn_to_ops REFUSES any layer it cannot bound -> 'skip'.
%     (2) ORIENTATION GUARD: the op-list must reproduce net.evaluate at several non-uniform
%         probe points (the bounds describe the SAME function the spec is about) -> else 'skip'.
%     (3) gpu_bab_crown_tight is a sound CROWN LOWER bound on C*y; G_d,i*y >= margins -> the
%         avoidance test margins-g > tol is a sufficient (sound) condition. DOUBLE precision
%         (FP64 ~1e-15) so rounding cannot make an avoidance spuriously hold.
%     (4) This pre-check only ever returns 'robust' on the PROVEN side; anything it cannot prove
%         is 'unknown'/'skip' -> Star. It can only ADD sound unsats, never a wrong verdict.
%
%   This first cut uses the ROOT tight CROWN bound (no BaB yet); benchmarks it cannot decide at
%   the root fall through to Star. ReLU-split BaB over the disjunctive predicate is the follow-on.

    info = struct('reason', '', 'nDisjuncts', 0, 'guardErr', NaN);
    if nargin < 5, opts = struct(); end
    guardTol = i_optget(opts, 'guardTol', 1e-4);
    precision = 'double';                              % certified path: FP64 (sound)

    lb = double(lb(:)); ub = double(ub(:));
    inShape = [numel(lb) 1];                           % control nets are flat FC feature inputs

    % (1+2) extract ops + orientation guard (multi-probe, non-uniform; deterministic, no rng).
    c = (lb + ub) / 2;
    n = numel(lb); idx = (1:n)';
    f1 = mod(idx * 0.6180339887498949, 1);
    f2 = mod(idx * 1.3247179572447460 + 0.37, 1);
    f3 = mod(idx * 0.7548776662466927 + 0.11, 1);
    probes = [c, lb + f1.*(ub-lb), lb + f2.*(ub-lb), lb + f3.*(ub-lb)];
    orders = {'colmajor', 'chw_rowmajor', 'hwc_rowmajor'};
    ops = []; nOut = NaN;
    for oi = 1:numel(orders)
        try
            cand = nn_to_ops(net, orders{oi});
        catch ME0
            if oi == 1, info.reason = ['nn_to_ops refused: ' ME0.message]; end
            continue;
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
        if ok, ops = cand; nOut = numel(yn); info.guardErr = maxe; info.reason = sprintf('flatten=%s', orders{oi}); break; end
    end
    if isempty(ops)
        if isempty(info.reason), info.reason = 'no flatten order matched net.evaluate (guard)'; end
        verdict = 'skip'; return;
    end

    % parse prop.Hg disjuncts -> per-disjunct (G,g); stack rows into one spec C for ONE bound pass
    [Gd, gd, rows, ok] = i_parse_disjuncts(prop, nOut);
    if ~ok
        verdict = 'skip'; info.reason = 'unsupported/empty Hg (could not parse disjuncts)'; return;
    end
    info.nDisjuncts = numel(Gd);
    C = cat(1, Gd{:});                                 % all rows of all disjuncts (nRows x nOut)
    gAll = cat(1, gd{:});                              % matching offsets

    % (3) sound disjunctive ReLU-split BaB (FP64). opts.spec routes the disjunct structure into the
    % batched BaB, which certifies SAFE iff EVERY disjunct is avoided (some row separated) -- and
    % SPLITS unstable ReLUs to tighten beyond the (typically too-loose) root CROWN. margin=guardTol
    % is the FP64 soundness slack. 'robust' is a sound UNSAT proof; anything else -> Star.
    spec = struct();
    spec.C = C; spec.g = gAll(:); spec.rowGroups = rows;
    bopts = struct('precision', precision, 'spec', spec, 'margin', guardTol, 'rootTight', true, ...
                   'maxNodes', i_optget(opts, 'maxNodes', 5000), 'maxFrontier', i_optget(opts, 'maxFrontier', 256), ...
                   'alphaIter', i_optget(opts, 'alphaIter', 20), 'betaIter', i_optget(opts, 'betaIter', 20));
    ev = getenv('NNV_HS_ALPHA'); if ~isempty(ev), bopts.alphaIter = str2double(ev); end
    ev = getenv('NNV_HS_BETA');  if ~isempty(ev), bopts.betaIter  = str2double(ev); end
    ev = getenv('NNV_HS_MAXNODES'); if ~isempty(ev), bopts.maxNodes = str2double(ev); end
    try
        [st, binfo] = gpu_bab_relu_split_batched(ops, lb, ub, 1, nOut, bopts);
    catch ME
        verdict = 'skip'; info.reason = ['halfspace BaB errored: ' ME.message]; return;
    end
    info.nodes = binfo.nodes;
    if strcmp(st, 'robust')
        verdict = 'robust';                            % every unsafe disjunct provably avoided over the box
    else
        verdict = 'unknown'; info.reason = sprintf('BaB %s (nodes=%d, rounds=%d)', st, binfo.nodes, binfo.rounds);
    end
end

% ---------------------------------------------------------------------------------------------
function [Gd, gd, rows, ok] = i_parse_disjuncts(prop, nOut)
% prop.Hg: array of HalfSpace objects, each one unsafe disjunct {y : G*y <= g}. Return per-disjunct
% G (k_d x nOut) and g (k_d x 1), plus rows{d} = the indices of disjunct d in the stacked spec.
    Gd = {}; gd = {}; rows = {}; ok = false;
    if ~isfield(prop, 'Hg') && ~isprop(prop, 'Hg'), return; end
    Hg = prop.Hg;
    if isempty(Hg), return; end
    if ~iscell(Hg), Hg = num2cell(Hg); end             % HalfSpace array -> cell of disjuncts
    r0 = 0;
    for d = 1:numel(Hg)
        hd = Hg{d};
        if ~(isprop(hd,'G') && isprop(hd,'g')), return; end
        G = double(hd.G); g = double(hd.g(:));
        if isempty(G) || size(G,2) ~= nOut || size(G,1) ~= numel(g), return; end
        Gd{end+1} = G; gd{end+1} = g; %#ok<AGROW>
        rows{end+1} = (r0 + 1) : (r0 + size(G,1)); %#ok<AGROW>
        r0 = r0 + size(G,1);
    end
    ok = ~isempty(Gd);
end

function v = i_optget(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
