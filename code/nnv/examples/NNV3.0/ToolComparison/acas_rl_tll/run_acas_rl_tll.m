function run_acas_rl_tll(varargin)
%RUN_ACAS_RL_TLL ToolComparison FC-net half: ACAS Xu, RL controllers, TLLverify.
%
%   Refresh of the NNV 2.0 / CAV'23 head-to-head on:
%       acas_p3   (45 ACAS Xu networks, property 3)
%       acas_p4   (45 ACAS Xu networks, property 4)
%       rl        (50 reinforcement-learning VNNLIB properties, fixed random subset)
%       tllverify (32 TLL VNNLIB properties)
%
%   Tools:
%       'nnv'         -> NNV 3.0 reachability
%       'mw_estimate' -> estimateNetworkOutputBounds + manual bound-check
%                        (requires AI Verification Toolbox)
%       'mw_abc'      -> verifyNetworkRobustness(net, vnnlib, Algorithm="alpha-beta-crown")
%                        (ACAS + TLLverify only; requires R2026a bridge)
%
%   NNV reach-method grid (CAV'23-style):
%       'approx-star'              fastest, baseline
%       'relax-star-range-25/50/75/100'   range-based relax with relaxFactor 0.25/0.5/0.75/1.0
%       'exact-star'               sound and complete; 8-core parallel
%       'relax-star-50'            legacy alias for 'relax-star-range-50' (back-compat
%                                  with persisted result rows from earlier runs)
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
    addParameter(p, 'benchmarks', {'acas_p3','acas_p4','rl','tllverify'});
    % mw_abc (alpha-beta-CROWN bridge) is intentionally absent from defaults:
    % the bridge requires R2026a. Override via 'tools' to include 'mw_abc'.
    addParameter(p, 'tools',      {'nnv','mw_estimate'});
    addParameter(p, 'algorithms', {});    % empty = use default algorithms_for(tool) per tool
    addParameter(p, 'timeout',    900);   % CAV'23 ACAS exact-star had max ~10000s -- cap outliers.
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
            case 'tllverify'
                run_tllverify(matFile, opts, u);
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
            algs = algorithms_for(tool, opts.algorithms);
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
            algs = algorithms_for(tool, opts.algorithms);
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

function run_tllverify(matFile, opts, u)
    if ~isempty(opts.rerun)
        n = u.purge_status(matFile, opts.rerun);
        if n > 0, fprintf("  purged %d rows with status in {%s}\n", n, strjoin(opts.rerun, ', ')); end
    end
    assetDir  = find_cav23_subdir('tllverify');
    csvFile   = fullfile(assetDir, 'instances.csv');
    T         = readtable(csvFile, 'Delimiter', ',', 'ReadVariableNames', false);
    T.Properties.VariableNames(1:3) = {'onnx_rel','vnnlib_rel','timeout'};

    for i = 1:height(T)
        onnx        = fullfile(assetDir, string(T.onnx_rel{i}));
        vnnlib      = fullfile(assetDir, string(T.vnnlib_rel{i}));
        instance_id = sprintf("%s|%s", T.onnx_rel{i}, T.vnnlib_rel{i});
        for t = 1:numel(opts.tools)
            tool = opts.tools{t};
            if strcmp(tool, 'mw_deeppoly'), continue; end
            algs = algorithms_for(tool, opts.algorithms);
            for k = 1:numel(algs)
                alg = algs{k};
                if u.has_instance(matFile, tool, 'tllverify', instance_id, alg), continue; end
                fprintf("  [tllverify %-8s %-22s] %s ... ", tool, alg, instance_id);
                [status, tsec] = run_one( @() verifyTLL(onnx, vnnlib, tool, alg), opts.timeout);
                row = u.new_row(tool, 'tllverify', instance_id, status, tsec, alg, opts.timeout);
                u.append_to_mat(matFile, row);
                fprintf("%s (%s s)\n", status, u.format_time(tsec, opts.timeout));
            end
        end
    end
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

function [status, tsec] = verifyTLL(onnx, vnnlibFile, tool, alg)
% Verify a single TLLverify instance.
    status = "error"; tsec = NaN; %#ok<NASGU>
    switch tool
        case 'nnv'
            nn       = matlab2nnv(rebuild_for_aivl(load_mw_network(onnx, 'BC')));
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

        case 'mw_abc'
            netMW = rebuild_for_aivl(load_mw_network(onnx, 'BC'));
            t = tic;
            r = verifyNetworkRobustness(netMW, vnnlibFile, Algorithm="alpha-beta-crown");
            status = mw_verifyresult_to_status(r);
            tsec   = toc(t);

        otherwise
            error("verifyTLL: tool %s not supported", tool);
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

function algs = algorithms_for(tool, override)
% Default algorithm set per tool (CAV'23-style grid for NNV).
    switch tool
        case 'nnv'
            algs = {'approx-star', ...
                    'relax-star-range-25', 'relax-star-range-50', ...
                    'relax-star-range-75', 'relax-star-range-100', ...
                    'exact-star'};
        case 'mw_estimate', algs = {'estimate-bounds'};
        case 'mw_deeppoly', algs = {'deep-poly'};
        case 'mw_abc',      algs = {'alpha-beta-crown'};
        otherwise, error("Unknown tool: %s", tool);
    end
    if nargin >= 2 && ~isempty(override)
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
