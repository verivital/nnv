function run_acas_rl_tll(varargin)
%RUN_ACAS_RL_TLL ToolComparison FC-net half: ACAS Xu + RL controllers.
%
%   Refresh of the NNV 2.0 / CAV'23 head-to-head on:
%       acas_p3   (45 ACAS Xu networks, property 3)
%       acas_p4   (45 ACAS Xu networks, property 4)
%       rl        (50 reinforcement-learning VNNLIB properties, fixed random subset)
%
%   (TLLverify was dropped from the active comparison because both tools'
%   fast methods give V=0 across the full grid; bundled CAV'23 rows are
%   preserved at acas_rl_tll/legacy/results_tllverify.{mat,csv}.)
%
%   Tools:
%       'nnv'         -> NNV 3.0 reachability
%       'mw_estimate' -> estimateNetworkOutputBounds + manual bound-check
%                        (requires AI Verification Library / AIVL)
%       'mw_abc'      -> verifyNetworkRobustness(net, vnnlib, Algorithm="alpha-beta-crown")
%                        (ACAS only; requires R2026a bridge)
%
%   NNV algorithm grid (per benchmark, mirrors NNV 2.0 / CAV'23 paper):
%       acas_p3 / acas_p4: full grid -- approx-star,
%                          relax-star-range-{25,50,75,100}, exact-star
%       rl (and other VNNLIB-style benchmarks): approx-star only
%
%   Override per call via 'algorithms' option. The legacy alias
%   'relax-star-50' maps to 'relax-star-range-50' for back-compat with
%   persisted result rows from earlier runs.
%
%   The mw_deeppoly variant of verifyNetworkRobustness only accepts argmax-style
%   robustness specs, not VNNLIB half-space properties; it's exercised in
%   mnist_resnet/run_mnist_resnet.m instead.
%
%   Results append incrementally to results/results_<benchmark>.mat with the
%   canonical schema from tool_utils. Instances already present in a result
%   file are skipped (resumable).
%
%   Examples:
%     run_acas_rl_tll()                                    % all benchmarks, default tools
%     run_acas_rl_tll('benchmarks',{'acas_p3','rl'})
%     run_acas_rl_tll('tools', {'nnv','mw_estimate'})
%     run_acas_rl_tll('algorithms', {'approx-star','relax-star-range-50'})
%     run_acas_rl_tll('timeout', 60, 'numNets', 5)         % smoke

    p = inputParser;
    addParameter(p, 'benchmarks', {'acas_p3','acas_p4','rl','oval21','collins_rul'});
    % mw_abc (alpha-beta-CROWN bridge) is intentionally absent from defaults:
    % the bridge requires R2026a. Override via 'tools' to include 'mw_abc'.
    addParameter(p, 'tools',      {'nnv','mw_estimate'});
    addParameter(p, 'algorithms', {});    % empty = use default algorithms_for(tool) per tool
    addParameter(p, 'timeout',    300);   % ACAS/RL global cap. CAV'23 had no cap (one outlier at 10479s);
                                          % we cap at 300s -- instances that need >300s also needed >900s
                                          % in prior runs, so paper-faithful T-counts are preserved while
                                          % wall-clock drops ~3x. VNNCOMP-style benchmarks (oval21,
                                          % collins_rul) use per-instance timeouts from instances.csv
                                          % directly (this default does NOT apply to them).
    addParameter(p, 'resultsDir', fullfile(fileparts(mfilename('fullpath')), 'results'));
    addParameter(p, 'rlSeed',     1);     % CAV'23 used rng(1) then randperm(200,50)
    addParameter(p, 'rerun',      {});    % statuses to re-run (e.g., {'error','timeout'})
    addParameter(p, 'numNets',    Inf);   % cap on # ACAS networks per benchmark (smoke: e.g. 5)
    parse(p, varargin{:});
    opts = p.Results;
    if ~isfolder(opts.resultsDir), mkdir(opts.resultsDir); end

    u = tool_utils();
    fprintf("acas_rl_tll: benchmarks = {%s}\n", strjoin(opts.benchmarks, ", "));
    fprintf("acas_rl_tll: tools      = {%s}\n", strjoin(opts.tools, ", "));
    fprintf("acas_rl_tll: timeout    = %g s\n", opts.timeout);
    fprintf("acas_rl_tll: resultsDir = %s\n\n", opts.resultsDir);

    for b = 1:numel(opts.benchmarks)
        bench   = opts.benchmarks{b};
        matFile = fullfile(opts.resultsDir, sprintf("results_%s.mat", bench));
        fprintf("=== %s -> %s ===\n", bench, matFile);
        switch bench
            case {'acas_p3','acas_p4'}
                run_acas(bench, matFile, opts, u);
            case 'rl'
                run_rl(matFile, opts, u);
            case 'oval21'
                run_image_vnnlib_csv('oval21', ...
                    fullfile(nnv_root(),'examples','NNV2.0','Submission','CAV2023','NNV_vs_MATLAB','oval21'), ...
                    matFile, opts, u);
            case 'collins_rul'
                run_image_vnnlib_csv('collins_rul', ...
                    fullfile(fileparts(mfilename('fullpath')), 'collins_rul_cnn_2022'), ...
                    matFile, opts, u);
            otherwise
                warning("Unknown benchmark: %s -- skipping", bench);
        end
    end

    fprintf("\nacas_rl_tll: done. Regenerate tables with:  make_acas_rl_tll_table\n");
end

% =========================================================================
% PER-BENCHMARK DRIVERS
% =========================================================================

function run_acas(bench, matFile, opts, u)
    if ~isempty(opts.rerun)
        n = u.purge_status(matFile, opts.rerun);
        if n > 0, fprintf("  purged %d rows with status in {%s}\n", n, strjoin(opts.rerun, ', ')); end
    end
    assetDir   = find_cav23_subdir('acas');
    acasOnnx   = fullfile(assetDir, 'onnx');
    propVnnlib = fullfile(assetDir, 'vnnlib', sprintf('prop_%s.vnnlib', last_char(bench)));
    [XLower, XUpper, H, g] = acas_property(bench);
    netFiles   = dir(fullfile(acasOnnx, '*.onnx'));
    if isfield(opts, 'numNets') && isfinite(opts.numNets)
        netFiles = netFiles(1:min(numel(netFiles), opts.numNets));
    end

    for i = 1:numel(netFiles)
        instance_id = string(netFiles(i).name);
        onnx        = fullfile(acasOnnx, char(instance_id));
        for t = 1:numel(opts.tools)
            tool = opts.tools{t};
            if strcmp(tool, 'mw_deeppoly'), continue; end     % N/A for ACAS half-spaces
            algs = algorithms_for(tool, bench, opts.algorithms);
            for k = 1:numel(algs)
                alg = algs{k};
                if u.has_instance(matFile, tool, bench, instance_id, alg), continue; end
                fprintf("  [%s %-12s %-22s] %s ... ", bench, tool, alg, instance_id);
                [status, tsec] = run_one( @() verifyAcas(onnx, propVnnlib, XLower, XUpper, H, g, tool, alg), opts.timeout);
                row = u.new_row(tool, bench, instance_id, status, tsec, alg, opts.timeout);
                u.append_to_mat(matFile, row);
                fprintf("%s (%s s)\n", status, u.format_time(tsec, opts.timeout));
            end
        end
    end
end

function run_rl(matFile, opts, u)
    if ~isempty(opts.rerun)
        n = u.purge_status(matFile, opts.rerun);
        if n > 0, fprintf("  purged %d rows with status in {%s}\n", n, strjoin(opts.rerun, ', ')); end
    end
    assetDir  = find_cav23_subdir('rl_benchmarks');
    csvFile   = fullfile(assetDir, 'instances.csv');
    T         = readtable(csvFile, 'Delimiter', ',', 'ReadVariableNames', false);
    T.Properties.VariableNames(1:3) = {'onnx_rel','vnnlib_rel','timeout'};
    rng(opts.rlSeed);
    inst_idxs = randperm(200, 50);

    for kk = 1:numel(inst_idxs)
        i           = inst_idxs(kk);
        onnx        = fullfile(assetDir, string(T.onnx_rel{i}));
        vnnlib      = fullfile(assetDir, string(T.vnnlib_rel{i}));
        instance_id = sprintf("%s|%s", T.onnx_rel{i}, T.vnnlib_rel{i});
        for t = 1:numel(opts.tools)
            tool = opts.tools{t};
            if ismember(tool, {'mw_deeppoly','mw_abc'}), continue; end
            algs = algorithms_for(tool, 'rl', opts.algorithms);
            for k = 1:numel(algs)
                alg = algs{k};
                if u.has_instance(matFile, tool, 'rl', instance_id, alg), continue; end
                fprintf("  [rl %-12s %-22s] %s ... ", tool, alg, instance_id);
                [status, tsec] = run_one( @() verifyRL(onnx, vnnlib, tool, alg), opts.timeout);
                row = u.new_row(tool, 'rl', instance_id, status, tsec, alg, opts.timeout);
                u.append_to_mat(matFile, row);
                fprintf("%s (%s s)\n", status, u.format_time(tsec, opts.timeout));
            end
        end
    end
end

function run_image_vnnlib_csv(bench, assetDir, matFile, opts, u)
% Generic driver for VNNCOMP-style image-input + VNNLIB-output benchmarks.
% Reads instances.csv (cols: onnx_rel, vnnlib_rel, timeout), iterates each
% instance through the (tool, algorithm) grid, and persists rows.
    if ~isempty(opts.rerun)
        n = u.purge_status(matFile, opts.rerun);
        if n > 0, fprintf("  purged %d rows with status in {%s}\n", n, strjoin(opts.rerun, ', ')); end
    end
    csvFile = fullfile(assetDir, 'instances.csv');
    T       = readtable(csvFile, 'Delimiter', ',', 'ReadVariableNames', false);
    T.Properties.VariableNames(1:3) = {'onnx_rel','vnnlib_rel','timeout'};

    for i = 1:height(T)
        onnxRel  = strip_dot_slash(string(T.onnx_rel{i}));
        vnnRel   = strip_dot_slash(string(T.vnnlib_rel{i}));
        onnx     = char(fullfile(assetDir, onnxRel));
        vnnlib   = char(fullfile(assetDir, vnnRel));
        instance_id = sprintf("%s|%s", onnxRel, vnnRel);
        % VNNCOMP semantics: each instances.csv row carries its own timeout
        % (30 s on easy props, up to 1800 s on hard ones). We honor it
        % directly -- opts.timeout governs ACAS/RL only, not VNNCOMP-style
        % benchmarks where the suite ships authoritative per-instance caps.
        instTimeout = double(T.timeout(i));
        if ~isfinite(instTimeout) || instTimeout <= 0, instTimeout = opts.timeout; end
        effTimeout = instTimeout;
        % Restart parpool between instances to avoid worker-state corruption.
        % We saw AIVL `verifyNetworkRobustness` hang on Collins RUL full_window
        % nets after NNV exact-star ran on the same parpool worker.
        delete(gcp('nocreate'));
        for t = 1:numel(opts.tools)
            tool = opts.tools{t};
            if ismember(tool, {'mw_deeppoly','mw_abc'}), continue; end
            algs = algorithms_for(tool, bench, opts.algorithms);
            for k = 1:numel(algs)
                alg = algs{k};
                if u.has_instance(matFile, tool, bench, instance_id, alg), continue; end
                fprintf("  [%s %-12s %-22s] %s ... ", bench, tool, alg, instance_id);
                [status, tsec] = run_one( @() verifyImageVnnlib(onnx, vnnlib, tool, alg), effTimeout);
                row = u.new_row(tool, bench, instance_id, status, tsec, alg, effTimeout);
                u.append_to_mat(matFile, row);
                fprintf("%s (%s s)\n", status, u.format_time(tsec, effTimeout));
            end
        end
    end
end

function s = strip_dot_slash(s)
    s = string(s);
    if startsWith(s, "./"), s = extractAfter(s, 2); end
end

function ensure_aivl_on_path()
% Make sure the manually-extracted AIVL Support Package wins over the
% matlabroot stub. parfeval workers don't auto-run startup.m, so we do
% it idempotently at each call. No-op if AIVL is already on path or not
% installed.
    persistent done
    if ~isempty(done) && done, return; end
    sps = dir(fullfile(userpath(), 'SupportPackages', 'R*', 'toolbox', 'nnet', 'supportpackages', 'aivnv'));
    for k = 1:numel(sps)
        addpath(fullfile(sps(k).folder, sps(k).name), '-begin');
    end
    done = true;
end

% =========================================================================
% VERIFY FUNCTIONS
% =========================================================================

function [status, tsec] = verifyAcas(onnx, vnnlibFile, XLower, XUpper, H, g, tool, alg)
% Verify a single ACAS Xu (network, property) with tool/algorithm.
    status = "error"; tsec = NaN; %#ok<NASGU>
    switch tool
        case 'nnv'
            netNNV  = matlab2nnv(load_mw_network(onnx, 'BCSS'));
            IS      = ImageStar(XLower, XUpper);
            reachOpt = reach_opt_for(alg);
            t = tic;
            R = netNNV.reach(IS, reachOpt);
            Rstar = [];
            for k = 1:numel(R), Rstar = [Rstar R(k).toStar]; end %#ok<AGROW>
            status = nnv_halfspace_status(Rstar, H, g);
            tsec   = toc(t);

        case 'mw_estimate'
            netMW = rebuild_for_aivl(load_mw_network(onnx, 'BCSS'));
            XL = dlarray(XLower, "CB"); XU = dlarray(XUpper, "CB");
            t = tic;
            [YLower, YUpper] = estimateNetworkOutputBounds(netMW, XL, XU);
            status = acas_mw_estimate_status(YLower, YUpper);
            tsec   = toc(t);

        case 'mw_abc'
            netMW = rebuild_for_aivl(load_mw_network(onnx, 'BCSS'));
            t = tic;
            r = verifyNetworkRobustness(netMW, vnnlibFile, Algorithm="alpha-beta-crown");
            status = mw_verifyresult_to_status(r);
            tsec   = toc(t);

        otherwise
            error("verifyAcas: tool %s not supported", tool);
    end
end

function [status, tsec] = verifyRL(onnx, vnnlibFile, tool, alg)
% Verify a single RL benchmark instance.
    status = "error"; tsec = NaN; %#ok<NASGU>
    switch tool
        case 'nnv'
            nn       = matlab2nnv(load_mw_network(onnx, 'BC'));
            property = load_vnnlib(vnnlibFile);
            IS       = Star(property.lb, property.ub);
            reachOpt = reach_opt_for(alg);
            t = tic;
            R = nn.reach(IS, reachOpt);
            status = nnv_rl_status(R, property.prop);
            tsec   = toc(t);

        case 'mw_estimate'
            netMW = rebuild_for_aivl(load_mw_network(onnx, 'BC'));
            [XLower, XUpper, output] = load_vnnlib_matlab(vnnlibFile);
            XL = dlarray(XLower, "CB"); XU = dlarray(XUpper, "CB");
            t = tic;
            [lb, ub] = estimateNetworkOutputBounds(netMW, XL, XU);
            status = mw_evalbounds_status(lb, ub, output);
            tsec   = toc(t);

        otherwise
            error("verifyRL: tool %s not supported", tool);
    end
end

function [status, tsec] = verifyImageVnnlib(onnx, vnnlibFile, tool, alg)
% Verify a single image-input + VNNLIB-output instance (oval21, collins_rul, ...).
%   Mirrors verifyRL but uses BCSS input + ImageStar (reshaped from VNNLIB lb/ub
%   to the network's image input dims). Same nnv_rl_status / mw_evalbounds_status
%   reuse the OR-of-half-spaces logic from CAV'23.
    status = "error"; tsec = NaN; %#ok<NASGU>
    switch tool
        case 'nnv'
            nn       = matlab2nnv(load_mw_network(onnx, 'BCSS'));
            property = load_vnnlib(vnnlibFile);
            inSize   = nn.Layers{1}.InputSize;
            lb       = reshape(property.lb, inSize);
            ub       = reshape(property.ub, inSize);
            IS       = ImageStar(lb, ub);
            reachOpt = reach_opt_for(alg);
            t = tic;
            R = nn.reach(IS, reachOpt);
            status = nnv_vnnlib_status(R, property.prop);
            tsec   = toc(t);

        case 'mw_estimate'
            % Use load_vnnlib (NNV-side parser) to get lb/ub + prop struct.
            % NNV's load_vnnlib_matlab parser has a paren-counting bug on
            % (or (and (<= Y_i Y_j))) outputs (oval21, collins_rul); this
            % path avoids it by evaluating the disjunctive halfspaces
            % manually against the AIVL output bounds.
            ensure_aivl_on_path();
            netMW    = rebuild_for_aivl(load_mw_network(onnx, 'BCSS'));
            property = load_vnnlib(vnnlibFile);
            inSize   = netMW.Layers(1).InputSize;
            XL = dlarray(reshape(property.lb, inSize), "SSCB");
            XU = dlarray(reshape(property.ub, inSize), "SSCB");
            t = tic;
            [yL, yU] = estimateNetworkOutputBounds(netMW, XL, XU);
            yL = extractdata(yL); yU = extractdata(yU);
            yL = yL(:); yU = yU(:);
            status = mw_vnnlib_status(yL, yU, property.prop);
            tsec   = toc(t);

        otherwise
            error("verifyImageVnnlib: tool %s not supported", tool);
    end
end

% =========================================================================
% STATUS NORMALIZATION
% =========================================================================

function status = nnv_halfspace_status(Set, H, g)
% Port of CAV'23 verifyNNV() from verifyP3.m. H*y <= g describes the UNSAFE region.
    status = "verified";
    for k = 1:numel(Set)
        S = Set(k).intersectHalfSpace(H, g);
        if isempty(S)
            continue;
        elseif isempty(Set(k).intersectHalfSpace(-H, -g))
            status = "violated"; return;
        else
            status = "unknown"; return;
        end
    end
end

function status = nnv_rl_status(R, prop)
% Port of verifyNNV() from verify_rl_nnv.m (supports OR-of-half-spaces).
    status = "unknown"; %#ok<NASGU>
    prop = prop{1};
    np   = numel(prop);
    nr   = numel(R);
    if np == 1
        status = "verified";
        for k = 1:nr
            Set = R(k);
            if isa(Set, 'ImageStar'), Set = Set.toStar; end
            S = Set.intersectHalfSpace(prop.Hg.G, prop.Hg.g);
            if isempty(S)
                status = "verified";
            elseif isempty(Set.intersectHalfSpace(-prop.Hg.G, -prop.Hg.g))
                status = "violated"; return;
            else
                status = "unknown";  return;
            end
        end
        return;
    end
    status = "verified";
    for cp = 1:np
        for k = 1:nr
            Set = R(k);
            if isa(Set, 'ImageStar'), Set = Set.toStar; end
            S = Set.intersectHalfSpace(prop(cp).Hg.G, prop(cp).Hg.g);
            if isempty(S), continue; end
            if isempty(Set.intersectHalfSpace(-prop(cp).Hg.G, -prop(cp).Hg.g))
                status = "violated"; return;
            else
                status = "unknown";
            end
        end
    end
end

function status = nnv_vnnlib_status(R, prop)
% VNNLIB-spec verifier supporting both shapes that load_vnnlib returns:
%   (a) RL-style: prop is a cell wrapping an ARRAY of structs, each with .Hg
%       a single halfspace struct (so np = numel(prop), prop(k).Hg.G).
%   (b) oval21-style: prop is a cell wrapping a SINGLE struct whose .Hg is
%       a struct array of halfspaces (so np = numel(prop.Hg), prop.Hg(k).G).
% Both encode an OR of half-spaces; we flatten to a list of (G, g) pairs.
    if iscell(prop), prop = prop{1}; end
    halves = struct('G',{},'g',{});
    if isstruct(prop) && isscalar(prop) && isfield(prop,'Hg')
        % oval21-style
        Hg = prop.Hg;
        for k = 1:numel(Hg), halves(end+1) = struct('G', Hg(k).G, 'g', Hg(k).g); end %#ok<AGROW>
    elseif isstruct(prop)
        % RL-style: array of structs, each has .Hg
        for k = 1:numel(prop)
            Hg = prop(k).Hg;
            for j = 1:numel(Hg), halves(end+1) = struct('G', Hg(j).G, 'g', Hg(j).g); end %#ok<AGROW>
        end
    else
        status = "unknown"; return;
    end

    np = numel(halves);
    nr = numel(R);
    Sets = cell(nr, 1);
    for k = 1:nr
        S = R(k);
        if isa(S, 'ImageStar'), S = S.toStar; end
        Sets{k} = S;
    end

    if np == 1
        h = halves(1);
        status = "verified";
        for k = 1:nr
            Inter = Sets{k}.intersectHalfSpace(h.G, h.g);
            if isempty(Inter), continue; end
            if isempty(Sets{k}.intersectHalfSpace(-h.G, -h.g))
                status = "violated"; return;
            end
            status = "unknown"; return;
        end
        return;
    end

    status = "verified";
    for cp = 1:np
        h = halves(cp);
        for k = 1:nr
            Inter = Sets{k}.intersectHalfSpace(h.G, h.g);
            if isempty(Inter), continue; end
            if isempty(Sets{k}.intersectHalfSpace(-h.G, -h.g))
                status = "violated"; return;
            end
            status = "unknown";
        end
    end
end

function status = mw_vnnlib_status(yL, yU, prop)
% Evaluate an OR-of-halfspaces VNNLIB property against AIVL output bounds.
% Used as a load_vnnlib_matlab-free path for benchmarks whose VNNLIB output
% specs trip up that parser (oval21, collins_rul). The property encodes
% counterexample conditions: status="violated" iff at least one half-space
% (G*y <= g) is forced over [yL,yU]; status="verified" iff ALL halfspaces
% are infeasible over [yL,yU]; otherwise unknown.
%
% Half-space feasibility test using interval arithmetic on G*y:
%   max(G*y - g) <= 0  ->  forced (definite counterexample)
%   min(G*y - g) >  0  ->  infeasible (no counterexample for this disjunct)
%
% Reuses nnv_vnnlib_status's halfspace flattening so RL-style and
% oval21-style prop encodings both work.
    if iscell(prop), prop = prop{1}; end
    halves = struct('G',{},'g',{});
    if isstruct(prop) && isscalar(prop) && isfield(prop, 'Hg')
        Hg = prop.Hg;
        for k = 1:numel(Hg), halves(end+1) = struct('G', Hg(k).G, 'g', Hg(k).g); end %#ok<AGROW>
    elseif isstruct(prop)
        for k = 1:numel(prop)
            Hg = prop(k).Hg;
            for j = 1:numel(Hg), halves(end+1) = struct('G', Hg(j).G, 'g', Hg(j).g); end %#ok<AGROW>
        end
    else
        status = "unknown"; return;
    end

    forced = false; allInfeasible = true;
    for k = 1:numel(halves)
        h = halves(k);
        Gpos = max(h.G, 0); Gneg = min(h.G, 0);
        Gy_max = Gpos * yU + Gneg * yL;
        Gy_min = Gpos * yL + Gneg * yU;
        gv = h.g(:);
        if all(Gy_max - gv <= 0)
            forced = true; break;
        end
        if ~all(Gy_min - gv > 0)
            allInfeasible = false;
        end
    end
    if forced
        status = "violated";
    elseif allInfeasible
        status = "verified";
    else
        status = "unknown";
    end
end

function status = acas_mw_estimate_status(YLower, YUpper)
% Port of verifyMAT() from CAV'23 verifyP3.m.
    YLower = extractdata(YLower); YUpper = extractdata(YUpper);
    if YLower(1) > YUpper(2) && YLower(1) > YUpper(3) && ...
       YLower(1) > YUpper(4) && YLower(1) > YUpper(5)
        status = "verified";
    elseif YUpper(1) <= YLower(2) && YUpper(1) <= YLower(3) && ...
           YUpper(1) <= YLower(4) && YUpper(1) <= YLower(5)
        status = "violated";
    else
        status = "unknown";
    end
end

function status = mw_evalbounds_status(lb, ub, output)
% Port of verifyMAT() from CAV'23 verify_rl_matlab.m / verify_tllverify_matlab.m.
    lb = extractdata(lb); ub = extractdata(ub); %#ok<NASGU>
    n = numel(output);
    r = ones(n, 1);
    for i = 1:n
        r(i) = eval(output{i}{1});
        if ~r(i)
            r(i) = eval(output{i}{2});
            if ~r(i), r(i) = 2; else,  r(i) = 0; break; end
        end
    end
    if all(r == 1),      status = "violated";
    elseif any(r == 0),  status = "verified";
    else,                status = "unknown";
    end
end

function status = mw_verifyresult_to_status(r)
    r = string(r);
    switch r
        case "verified", status = "verified";
        case "violated", status = "violated";
        otherwise,       status = "unknown";
    end
end

% =========================================================================
% MISC HELPERS
% =========================================================================

function algs = algorithms_for(tool, benchmark, override)
% Default algorithm set per (tool, benchmark). Mirrors the NNV 2.0 / CAV'23
% paper methodology: ACAS p3/p4 use the full grid (CAV'23 Tables 2-3); RL
% and other VNNLIB-style benchmarks were originally evaluated on approx-star
% only. Override by passing a non-empty 'algorithms' option to the driver.
    switch tool
        case 'nnv'
            full_grid_range = {'approx-star', ...
                               'relax-star-range-25', 'relax-star-range-50', ...
                               'relax-star-range-75', 'relax-star-range-100', ...
                               'exact-star'};
            full_grid_area  = {'approx-star', ...
                               'relax-star-area-25', 'relax-star-area-50', ...
                               'relax-star-area-75', 'relax-star-area-100'};
            switch string(benchmark)
                case {"acas_p3","acas_p4"}
                    algs = full_grid_range;
                case "oval21"
                    % CAV'23 used approx-star only; we extend with the area-relax
                    % grid for the ATVA26 paper. exact-star is documented as
                    % intractable on CIFAR Conv+ReLU networks; not in defaults.
                    algs = full_grid_area;
                case "collins_rul"
                    % Small RUL CNNs; exact-star is tractable here (~few s/inst).
                    algs = [full_grid_area, {'exact-star'}];
                otherwise
                    algs = {'approx-star'};
            end
        case 'mw_estimate', algs = {'estimate-bounds'};
        case 'mw_deeppoly', algs = {'deep-poly'};
        case 'mw_abc',      algs = {'alpha-beta-crown'};
        otherwise, error("Unknown tool: %s", tool);
    end
    if nargin >= 3 && ~isempty(override)
        algs = algs(ismember(algs, override));
    end
end

function reachOpt = reach_opt_for(alg)
% NNV reachability options keyed by algorithm string. All variants use
% reach-method 'relax-star-range' (per CAV'23) for FC nets. The legacy
% alias 'relax-star-50' maps to the canonical 'relax-star-range-50' so
% existing persisted result rows keep matching.
    reachOpt = struct;
    reachOpt.reachOption = "single";
    reachOpt.numCores    = 1;
    reachOpt.relaxFactor = 0;
    switch alg
        case 'approx-star'
            reachOpt.reachMethod = 'approx-star';
        case {'relax-star-range-25'}
            reachOpt.reachMethod = 'relax-star-range';
            reachOpt.relaxFactor = 0.25;
        case {'relax-star-range-50','relax-star-50'}    % alias kept for back-compat
            reachOpt.reachMethod = 'relax-star-range';
            reachOpt.relaxFactor = 0.5;
        case {'relax-star-range-75'}
            reachOpt.reachMethod = 'relax-star-range';
            reachOpt.relaxFactor = 0.75;
        case {'relax-star-range-100'}
            reachOpt.reachMethod = 'relax-star-range';
            reachOpt.relaxFactor = 1.0;
        case 'relax-star-area-25'
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 0.25;
        case 'relax-star-area-50'
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 0.5;
        case 'relax-star-area-75'
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 0.75;
        case 'relax-star-area-100'
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 1.0;
        case 'exact-star'
            reachOpt.reachMethod = 'exact-star';
            reachOpt.reachOption = "parallel";
            reachOpt.numCores    = min(8, maxNumCompThreads);
        otherwise
            error("Unknown NNV algorithm: %s", alg);
    end
end

function [XLower, XUpper, H, g] = acas_property(bench)
    H = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
    g = [0; 0; 0; 0];
    switch bench
        case 'acas_p3'
            XLower = [-0.303531156; -0.009549297; 0.493380324; 0.3; 0.3];
            XUpper = [-0.298552812;  0.009549297; 0.5;        0.5; 0.5];
        case 'acas_p4'
            XLower = [-0.303531156; -0.009549297; 0.0;         0.318181818; 0.083333333];
            XUpper = [-0.298552812;  0.009549297; 0.0;         0.5;         0.166666667];
    end
end

function net = load_mw_network(onnx, inputFmt)
% R2025b-compatible ONNX loader. Prefer importNetworkFromONNX (returns
% dlnetwork with affine folds). Fall back to the manual ElementwiseAffineLayer
% fold from CAV'23 on older opsets. Accepts inputFmt in {'BCSS','BC'}.
%
% importNetworkFromONNX writes auto-generated custom layers to a `+pkgname`
% directory in the current working directory. When the cwd is read-only (as
% with a bind-mounted host source tree owned by a different UID), the import
% fails. We `cd` to tempdir for the import and restore the original cwd.
    oldDir   = pwd;
    importDir = tempname; mkdir(importDir);
    cleanupCd  = onCleanup(@() safe_cd(oldDir));
    cleanupRm  = onCleanup(@() safe_rmdir(importDir));
    cd(importDir);
    try
        net = importNetworkFromONNX(onnx, InputDataFormats=string(inputFmt));
        return;
    catch
    end
    try
        net = importONNXNetwork(onnx, InputDataFormats=string(inputFmt));
        Layers = net.Layers;
        n = numel(Layers);
        keep = [];
        for i = 1:(n-1)
            if isa(Layers(i), 'nnet.onnx.layer.ElementwiseAffineLayer') && i > 1
                Layers(i-1).Bias = Layers(i).Offset;
            else
                keep(end+1) = i; %#ok<AGROW>
            end
        end
        Layers = Layers(keep);
        net    = dlnetwork(Layers);
    catch ME2
        error("load_mw_network: both importNetworkFromONNX and importONNXNetwork failed: %s", ME2.message);
    end
end

function safe_cd(d), try, cd(d); catch, end, end %#ok<TRYNC>
function safe_rmdir(d), try, rmdir(d, 's'); catch, end, end %#ok<TRYNC>

function d = find_cav23_subdir(sub)
% Find NNV_vs_MATLAB/<sub> under the NNV examples tree.
    candidates = { ...
        fullfile(nnv_root(), 'examples', 'NNV2.0', 'Submission', 'CAV2023', 'NNV_vs_MATLAB', sub), ...
        fullfile('/home/matlab/nnv/code/nnv/examples/NNV2.0/Submission/CAV2023/NNV_vs_MATLAB', sub)};
    for i = 1:numel(candidates)
        if isfolder(candidates{i}), d = candidates{i}; return; end
    end
    error("Cannot locate CAV'23 %s directory; tried:\n  %s", sub, strjoin(candidates, [newline '  ']));
end

function r = nnv_root()
% This file lives at code/nnv/examples/NNV3.0/ToolComparison/acas_rl_tll/run_acas_rl_tll.m
% nnv_root is code/nnv/ (4 dir-levels up from this file's directory).
    here = fileparts(mfilename('fullpath'));      % .../acas_rl_tll
    r = fileparts(fileparts(fileparts(fileparts(here))));   % .../code/nnv
end

function c = last_char(s), c = s(end); end

function [status, tsec] = run_one(fn, timeout)
% Execute fn under a wall-clock timeout via parfeval (process pool).
    persistent havePar
    if isempty(havePar)
        havePar = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    end
    if havePar && timeout > 0 && isfinite(timeout)
        pool = gcp('nocreate'); if isempty(pool), pool = parpool('local'); end
        F = parfeval(pool, fn, 2);
        done = wait(F, 'finished', timeout);
        if ~done
            cancel(F);
            status = "timeout"; tsec = timeout; return;
        end
        try
            [status, tsec] = fetchOutputs(F);
        catch ME
            status = "error"; tsec = NaN;
            fprintf(2, "\n    run_one fetchOutputs error:\n%s\n", getReport(ME, 'extended'));
        end
    else
        try
            [status, tsec] = fn();
        catch ME
            status = "error"; tsec = NaN;
            fprintf(2, "\n    run_one error:\n%s\n", getReport(ME, 'extended'));
        end
    end
end
