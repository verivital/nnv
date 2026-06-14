function [status, counterEx] = verify_multinet(nnvnet, property, reachOptions)
    % [status, counterEx] = verify_multinet(nnvnet, property, reachOptions)
    %
    % Verify a VNN-LIB 2.0 multi-network `equal-to` property (property.multinet
    % as produced by load_vnnlib2, Phase 2) on ONE shared network: `equal-to`
    % means both declared networks are the SAME model, evaluated on two inputs
    % coupled by property.multinet.jointC/jointd.
    %
    % Inputs:
    %   nnvnet       - the (single, shared) NN object
    %   property     - load_vnnlib2 output with property.multinet present
    %   reachOptions - NNV reach options struct; defaults to
    %                  struct('reachMethod','exact-star') (plan section 3.1:
    %                  approx-star decorrelates Y_f/Y_g and typically yields
    %                  unknown; exact-star keeps the joint correlation)
    %
    % Outputs (VNN-COMP convention, mirrors run_vnncomp_instance):
    %   status    = 0 -> sat: counterexample found AND validated by evaluating
    %                    the ORIGINAL network on both halves + replaying ALL
    %                    parsed input constraints on the concrete pair
    %             = 1 -> unsat: the joint product-net reach set is provably
    %                    disjoint from the unsafe cross-network region
    %             = 2 -> unknown: anything else (non-whitelisted layer, shape
    %                    mismatch, any thrown error, any doubt)
    %   counterEx = {x_f, x_g} (flat column vectors, parser row-major ONNX
    %               order) when status == 0; [] otherwise.
    %
    % SOUNDNESS (-150 rule): every path defaults to unknown. The whole body is
    % wrapped in try/catch; a verdict is only emitted from the two defensible
    % paths above. Reach returning "intersection nonempty" is NOT promoted to
    % sat (even under exact-star) because that would require extracting and
    % replay-validating a concrete witness from the intersecting star.

    status = 2; counterEx = [];
    if nargin < 3 || isempty(reachOptions) || ~isstruct(reachOptions)
        reachOptions = struct('reachMethod', 'exact-star');
    end
    if ~isfield(reachOptions, 'reachMethod') || isempty(reachOptions.reachMethod)
        % NN.validate_reach_options reads reachMethod unconditionally
        reachOptions.reachMethod = 'exact-star';
    end
    try
        [status, counterEx] = verify_multinet_impl(nnvnet, property, reachOptions);
    catch
        % ANY error -> unknown (never let an exception escape into the runner)
        status = 2; counterEx = [];
    end
end

% =========================================================================

function [status, counterEx] = verify_multinet_impl(nnvnet, property, reachOptions)
    status = 2; counterEx = [];

    % ---- 1) property validation (anything off -> unknown) ------------------
    % FALSIFICATION-FIRST ordering: the sat path only needs nnvnet.evaluate (which
    % handles the FULL imported net -- ImageInput/ElementwiseAffine/Flatten wrappers
    % and all), so validate the property and try to falsify BEFORE assembling the
    % FC+ReLU-only product net. The product net (and its whitelist) is built later and
    % is needed ONLY for the unsat/reach direction. This recovers easily-falsifiable
    % equal-to instances (measured on monotonic_acasxu: ~50% of feasible samples
    % violate) that the old whitelist-first order discarded as unknown.
    if ~isstruct(property) || ~isfield(property, 'multinet'), return; end
    if isfield(property, 'unsupported') && property.unsupported, return; end
    mn = property.multinet;
    req = {'names', 'inShapes', 'outShapes', 'equivKind', 'jointLb', 'jointUb', ...
           'jointC', 'jointd', 'crossProp'};
    for k = 1:numel(req)
        if ~isfield(mn, req{k}), return; end
    end
    if ~strcmp(mn.equivKind, 'equal')
        % isomorphic-to (different weights) must NEVER fall through to the
        % equal-to falsifier/product: g is a DIFFERENT network there
        return;
    end
    jlb = double(mn.jointLb(:)); jub = double(mn.jointUb(:));
    twoN = numel(jlb);
    if twoN == 0 || numel(jub) ~= twoN || mod(twoN, 2) ~= 0, return; end
    if any(~isfinite(jlb)) || any(~isfinite(jub)) || any(jlb > jub), return; end
    n = twoN/2;
    if prod(mn.inShapes{1}) ~= n || prod(mn.inShapes{2}) ~= n, return; end
    outDim = prod(mn.outShapes{1});
    if outDim < 1 || prod(mn.outShapes{2}) ~= outDim, return; end
    Cj = double(mn.jointC); dj = double(mn.jointd(:));
    if isempty(Cj)
        Cj = zeros(0, twoN); dj = zeros(0, 1);
    end
    if size(Cj, 2) ~= twoN || size(Cj, 1) ~= numel(dj), return; end
    hs = mn.crossProp;
    if isempty(hs) || ~isa(hs, 'HalfSpace'), return; end
    for h = 1:numel(hs)
        if size(hs(h).G, 2) ~= 2*outDim || size(hs(h).G, 1) ~= numel(hs(h).g)
            return;
        end
    end
    % network I/O widths must match the spec (first FC consumes the input,
    % last FC produces the output -- whitelist guarantees both exist)
    fcIdx = find(cellfun(@(L) isa(L, 'FullyConnectedLayer'), nnvnet.Layers));
    if isempty(fcIdx), return; end
    if nnvnet.Layers{fcIdx(1)}.InputSize ~= n, return; end
    if nnvnet.Layers{fcIdx(end)}.OutputSize ~= outDim, return; end

    % ---- 3) falsification FIRST (sound, validated sat witness) -------------
    % rejection-sample the joint box; jointC/jointd constrain the CONCRETE
    % stacked x = [x_f; x_g] DIRECTLY (parser convention -- see
    % apply_mn_input_expr in load_vnnlib2.m), so feasibility is checked on x
    % itself, NOT on star predicates. `==` couplings are measure-zero under
    % box sampling, so detected equality pairs are REPAIRED algebraically
    % before the feasibility check (the repaired x satisfies them exactly).
    nSamples = 1000;
    if isfield(reachOptions, 'falsifySamples') && ~isempty(reachOptions.falsifySamples)
        nSamples = reachOptions.falsifySamples;
    end
    [eqP, eqQ, eqCP, eqCQ, eqD] = detect_equality_pairs(Cj, dj);
    for s = 1:nSamples
        x = jlb + rand(twoN, 1).*(jub - jlb);
        % repair two-term equality couplings: a_p*x_p + a_q*x_q = d
        %   -> x_q = (d - a_p*x_p)/a_q  (exact for the parser's +/-1 rows)
        for k = 1:numel(eqP)
            x(eqQ(k)) = (eqD(k) - eqCP(k)*x(eqP(k)))/eqCQ(k);
        end
        if any(x < jlb) || any(x > jub), continue; end      % repair left the box
        if ~isempty(Cj) && any(Cj*x > dj), continue; end    % infeasible -> reject
        xf = x(1:n); xg = x(n+1:end);
        yf = double(nnvnet.evaluate(xf)); yf = yf(:);
        yg = double(nnvnet.evaluate(xg)); yg = yg(:);
        if numel(yf) ~= outDim || numel(yg) ~= outDim
            return;   % network output does not match the declared shape -> unknown
        end
        y = [yf; yg];
        for h = 1:numel(hs)
            % STRICT interior check (G*y < g componentwise): crossProp stores
            % the <=-relaxation of possibly-STRICT asserts, so a boundary point
            % (G*y == g) could fail the original strict assert on the official
            % witness replay. Requiring strict interiority is sound (at most we
            % miss boundary witnesses -> fall through to reach/unknown).
            if all(double(hs(h).G)*y - double(hs(h).g) < 0)
                status = 0;
                counterEx = {xf, xg};
                return;
            end
        end
    end

    % ---- product net (FC+ReLU whitelist) -- needed ONLY for unsat/reach ----
    % every layer must be FullyConnectedLayer or ReluLayer; anything else cannot be
    % soundly stacked into [f(x_f); f(x_g)] -> no unsat proof -> unknown (the sat
    % path above already ran on the original net regardless of this).
    [prodNet, okp] = multinet_product_net(nnvnet);
    if ~okp
        return;
    end

    % ---- 4) product-net consistency cross-check ----------------------------
    % the unsat path trusts h([xf;xg]) == [f(xf); f(xg)]; verify it CONCRETELY
    % on a few box samples (no feasibility needed -- evaluation equivalence
    % must hold pointwise everywhere). Guards layer-order/Connections quirks.
    for t = 1:3
        x = jlb + rand(twoN, 1).*(jub - jlb);
        y1 = double(nnvnet.evaluate(x(1:n)));
        y2 = double(nnvnet.evaluate(x(n+1:end)));
        yref = [y1(:); y2(:)];
        yprod = double(prodNet.evaluate(x)); yprod = yprod(:);
        if numel(yprod) ~= numel(yref) || norm(yprod - yref) > 1e-4*max(1, norm(yref))
            return;   % product does not replicate f twice -> unknown
        end
    end

    % ---- 5) joint input Star ------------------------------------------------
    % Star(lb, ub) (via Box.toStar) builds  x = c + G*alpha  with
    %     c = (lb+ub)/2 = V(:,1),
    %     G = V(:,2:end) = diag((ub-lb)/2) with all-zero (PINNED, lb==ub)
    %         columns REMOVED,
    %     alpha in [-1, 1]^m  (m = number of non-degenerate dims).
    % So alpha is NOT the raw input offset: it is scaled by the half-width and
    % pinned dims have no predicate at all. The parser's jointC/jointd
    % constrain the CONCRETE x; substituting x = c + G*alpha gives the
    % predicate-space constraints to append:
    %     jointC*x <= jointd  <=>  (jointC*G)*alpha <= jointd - jointC*c
    % Composing through S0.V (instead of assuming x = c + I*alpha) keeps the
    % algebra correct regardless of which columns Box dropped.
    S0 = Star(jlb, jub);
    c0 = double(S0.V(:, 1));
    G0 = double(S0.V(:, 2:end));
    Ca = Cj*G0;
    da = dj - Cj*c0;
    % rows whose variables are all pinned map to all-zero alpha rows
    % (0*alpha <= da(r)):
    %   da(r) >= 0 -> trivially true, drop the row;
    %   da(r) <  0 -> the input region is empty; technically a vacuous unsat
    %                 but we refuse to claim a verdict from an empty region
    zr = ~any(Ca, 2);
    if any(da(zr) < 0)
        return;
    end
    Ca = Ca(~zr, :); da = da(~zr);
    S = Star(double(S0.V), [double(S0.C); Ca], [double(S0.d); da], ...
             double(S0.predicate_lb), double(S0.predicate_ub));

    % ---- 6) joint reach + spec check ----------------------------------------
    R = prodNet.reach(S, reachOptions);
    spec = struct();
    spec.Hg = hs;   % field assignment (struct('Hg', hs) would explode an hs array)
    res = verify_specification(R, {spec});
    if res == 1
        % reach set (exact or over-approximate -- both sound for this
        % direction) does not intersect the unsafe cross region -> unsat
        status = 1;
    else
        % res == 2 is NOT promoted to sat even under exact-star: that would
        % require extracting a concrete witness from the intersecting star and
        % replay-validating it. Stay unknown.
        status = 2;
    end
end

% =========================================================================

function [eqP, eqQ, eqCP, eqCQ, eqD] = detect_equality_pairs(Cj, dj)
    % find two-row inequality pairs encoding equalities: rows (i, j) with
    % Cj(j,:) == -Cj(i,:) and dj(j) == -dj(i) represent  Cj(i,:)*x == dj(i).
    % Only two-term rows are returned (the parser emits exactly those); the
    % falsifier uses them to repair sampled points onto the equality manifold.
    eqP = zeros(0, 1); eqQ = zeros(0, 1);
    eqCP = zeros(0, 1); eqCQ = zeros(0, 1); eqD = zeros(0, 1);
    m = size(Cj, 1);
    used = false(m, 1);
    for i = 1:m
        if used(i), continue; end
        nzi = find(Cj(i, :));
        if numel(nzi) ~= 2, continue; end
        for j = i+1:m
            if used(j), continue; end
            if isequal(Cj(j, :), -Cj(i, :)) && dj(j) == -dj(i)
                used(i) = true; used(j) = true;
                eqP(end+1, 1) = nzi(1);          %#ok<AGROW>
                eqQ(end+1, 1) = nzi(2);          %#ok<AGROW>
                eqCP(end+1, 1) = Cj(i, nzi(1));  %#ok<AGROW>
                eqCQ(end+1, 1) = Cj(i, nzi(2));  %#ok<AGROW>
                eqD(end+1, 1) = dj(i);           %#ok<AGROW>
                break;
            end
        end
    end
end
