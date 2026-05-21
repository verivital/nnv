function run_vnnlib_half(opts, u)
%RUN_VNNLIB_HALF  driver for the five VNNLIB-style benchmarks:
%   acas_xu_p3, acas_xu_p4, rl, oval21, collins_rul
%
%   opts struct fields:
%     mode       'smoke' | 'default'
%     tools      cell, subset of {'nnv','aivl'}
%     algorithms cell, optional algorithm filter (empty = use per-bench defaults)
%
%   u: the tool_utils() handle.

    here = fileparts(mfilename('fullpath'));
    tc_root = fileparts(here);
    results_dir = fullfile(tc_root, 'results');
    if ~isfolder(results_dir), mkdir(results_dir); end

    BENCHES = {'acas_xu_p3','acas_xu_p4','rl','collins_rul','oval21'};

    % Optional benchmarks filter (passed in via opts.benchmarks).
    if isfield(opts, 'benchmarks') && ~isempty(opts.benchmarks)
        BENCHES = BENCHES(ismember(BENCHES, opts.benchmarks));
    end

    for b = 1:numel(BENCHES)
        bench = BENCHES{b};
        assetDir = fullfile(tc_root, 'benchmarks', bench);
        matFile  = fullfile(results_dir, sprintf('%s.mat', bench));

        if ~isfolder(assetDir)
            fprintf('[%s] SKIP — asset dir missing: %s\n', bench, assetDir);
            continue;
        end

        % Decompress any .gz in onnx/ and vnnlib/ that haven't been touched yet.
        decompress_gz_dir(fullfile(assetDir, 'onnx'),   '*.onnx.gz');
        decompress_gz_dir(fullfile(assetDir, 'vnnlib'), '*.vnnlib.gz');

        csvFile = fullfile(assetDir, 'instances.csv');
        T = readtable(csvFile, 'Delimiter', ',', 'ReadVariableNames', false);
        T.Properties.VariableNames(1:3) = {'onnx_rel','vnnlib_rel','timeout'};

        % Smoke mode: 1 instance per benchmark.
        if strcmp(opts.mode, 'smoke')
            T = T(1, :);
        end

        for i = 1:height(T)
            onnxRel = strip_dot_slash(string(T.onnx_rel{i}));
            vnnRel  = strip_dot_slash(string(T.vnnlib_rel{i}));
            onnxRel = regexprep(onnxRel, '\.gz$', '');
            vnnRel  = regexprep(vnnRel,  '\.gz$', '');
            onnx    = char(fullfile(assetDir, onnxRel));
            vnnlib  = char(fullfile(assetDir, vnnRel));
            instance_id = sprintf('%s|%s', onnxRel, vnnRel);

            inst_to = double(T.timeout(i));
            if ~isfinite(inst_to) || inst_to <= 0, inst_to = 300; end

            for t = 1:numel(opts.tools)
                tool = opts.tools{t};
                algs = algorithms_for(tool, bench, opts.mode, opts.algorithms);
                for k = 1:numel(algs)
                    alg = algs{k};
                    if u.has_instance(matFile, tool, bench, instance_id, alg)
                        continue;
                    end
                    fprintf('  [%s %-4s %-22s] %s ... ', bench, tool, alg, instance_id);
                    [status, tsec] = run_with_timeout( ...
                        @() verify_instance(onnx, vnnlib, bench, tool, alg), inst_to);
                    row = u.new_row(tool, bench, instance_id, status, tsec, alg, inst_to);
                    u.append_to_mat(matFile, row);
                    fprintf('%s (%s s)\n', status, u.format_time(tsec, inst_to));
                end
            end
        end
    end
end

% =========================================================================
% Algorithm grid
% =========================================================================

function algs = algorithms_for(tool, bench, mode, override)
%ALGORITHMS_FOR Per-(tool, benchmark, mode) algorithm list.
    switch tool
        case 'nnv'
            switch bench
                case 'acas_xu_p3'
                    full = {'approx-star','exact-star','relax-star-range-50'};
                case 'acas_xu_p4'
                    full = {'approx-star','exact-star','relax-star-range-50'};
                case 'oval21'
                    full = {'approx-star'};
                case 'rl'
                    full = {'approx-star','exact-star','relax-star-range-50'};
                case 'collins_rul'
                    full = {'approx-star'};
                otherwise
                    full = {'approx-star'};
            end
        case 'aivl'
            % All five use estimate-bounds. deep-poly with VNNLIB-form
            % output specs requires verifyNetworkRobustness with VNNLIB ingest,
            % which lands in R2026a. Under R2025b the AIVL baseline is
            % estimate-bounds + manual half-space check.
            full = {'estimate-bounds'};
        otherwise
            error('algorithms_for: unknown tool %s', tool);
    end

    if strcmp(mode, 'smoke')
        if strcmp(tool, 'nnv'), algs = full(1);
        else, algs = {};   % skip AIVL in smoke (avoid warmup latency)
        end
    else
        algs = full;
    end

    if nargin >= 4 && ~isempty(override)
        algs = algs(ismember(algs, override));
    end
end

% =========================================================================
% Per-instance verification
% =========================================================================

function [status, tsec] = verify_instance(onnx, vnnlib, bench, tool, alg)
%VERIFY_INSTANCE Dispatch to per-benchmark verification logic.
%   Bench config:
%     acas_xu_p3  — image (BCSS) → ImageStar (half-space VNNLIB)
%     acas_xu_p4  — image (BCSS) → ImageStar (half-space VNNLIB)
%     oval21      — image (BCSS) → ImageStar (half-space VNNLIB)
%     rl          — FC (BC)      → Star      (half-space VNNLIB)
%     collins_rul — image (BCSS) → ImageStar (half-space VNNLIB)
    status = "error"; tsec = NaN; %#ok<NASGU>
    switch bench
        case 'acas_xu_p3'
            [status, tsec] = verify_image_vnnlib(onnx, vnnlib, tool, alg, 'BCSS', '');
        case 'acas_xu_p4'
            [status, tsec] = verify_image_vnnlib(onnx, vnnlib, tool, alg, 'BCSS', '');
        case 'oval21'
            % CIFAR ResNet — BCSS input, ImageStar wrap, half-space VNNLIB.
            [status, tsec] = verify_image_vnnlib(onnx, vnnlib, tool, alg, 'BCSS', '');
        case 'rl'
            % FC input (BC), Star wrap, half-space VNNLIB output.
            [status, tsec] = verify_fc_vnnlib(onnx, vnnlib, tool, alg, '', 'Star');
        case 'collins_rul'
            % Small 1D CNN — BCSS input, ImageStar wrap, half-space VNNLIB.
            [status, tsec] = verify_image_vnnlib(onnx, vnnlib, tool, alg, 'BCSS', '');
        otherwise
            error('verify_instance: unknown benchmark %s', bench);
    end
end

function [status, tsec] = verify_fc_vnnlib(onnx, vnnlib, tool, alg, outFmt, ~)
%FC-input + VNNLIB output. rl benchmark only.
%   outFmt: '' or 'BC' — passed through to load_mw_network.
    if nargin < 5, outFmt = ''; end
    status = "error"; tsec = NaN; %#ok<NASGU>
    switch tool
        case 'nnv'
            nn       = matlab2nnv(load_mw_network(onnx, 'BC', outFmt));
            property = load_vnnlib(vnnlib);
            IS       = Star(property.lb, property.ub);
            reachOpt = reach_opt_for(alg);
            t = tic;
            R = nn.reach(IS, reachOpt);
            status = nnv_vnnlib_status_local(R, property.prop);
            tsec   = toc(t);

        case 'aivl'
            netMW = rebuild_for_aivl(load_mw_network(onnx, 'BC', outFmt));
            switch alg
                case 'estimate-bounds'
                    [XLower, XUpper, output] = load_vnnlib_matlab(vnnlib);
                    XL = dlarray(XLower, 'CB'); XU = dlarray(XUpper, 'CB');
                    t = tic;
                    [lb, ub] = estimateNetworkOutputBounds(netMW, XL, XU);
                    status = aivl_evalbounds_status_local(lb, ub, output);
                    tsec   = toc(t);
                otherwise
                    error('verify_fc_vnnlib: aivl alg %s not supported', alg);
            end
        otherwise
            error('verify_fc_vnnlib: tool %s not supported', tool);
    end
end

function [status, tsec] = verify_image_vnnlib(onnx, vnnlib, tool, alg, fmt, outFmt)
%Image-input + VNNLIB output. acas_xu_p3.
    if nargin < 6, outFmt = ''; end
    status = "error"; tsec = NaN; %#ok<NASGU>
    switch tool
        case 'nnv'
            nn       = matlab2nnv(load_mw_network(onnx, fmt, outFmt));
            property = load_vnnlib(vnnlib);
            inSize   = nn.Layers{1}.InputSize;
            lb       = reshape(property.lb, inSize);
            ub       = reshape(property.ub, inSize);
            IS       = ImageStar(lb, ub);
            reachOpt = reach_opt_for(alg);
            t = tic;
            R = nn.reach(IS, reachOpt);
            status = nnv_vnnlib_status_local(R, property.prop);
            tsec   = toc(t);

        case 'aivl'
            ensure_aivl_on_path_local();
            netMW    = rebuild_for_aivl(load_mw_network(onnx, fmt, outFmt));
            switch alg
                case 'estimate-bounds'
                    property = load_vnnlib(vnnlib);
                    firstLayer = netMW.Layers(1);
                    isFeature  = isa(firstLayer, 'nnet.cnn.layer.FeatureInputLayer') ...
                              || isa(firstLayer, 'nnet.onnx.layer.FeatureInputLayer');
                    if isFeature
                        % R2025b often imports ACAS-style nets as
                        % FeatureInputLayer (0 spatial dims) even when we
                        % requested 'BCSS'. AIVL rejects 'SSCB' for these
                        % ("Layer expects 0 spatial dims, received 2"). Use
                        % 'CB' for flat-input nets.
                        XL = dlarray(property.lb(:), 'CB');
                        XU = dlarray(property.ub(:), 'CB');
                    else
                        % Image input layer: reshape lb/ub to the layer's
                        % [H W C] InputSize and use 'SSCB'.
                        inSize = firstLayer.InputSize;
                        if isscalar(inSize),         inSize = [1, 1, inSize];      end
                        if numel(inSize) == 2,       inSize = [inSize, 1];          end
                        XL = dlarray(reshape(property.lb, inSize), 'SSCB');
                        XU = dlarray(reshape(property.ub, inSize), 'SSCB');
                    end
                    t = tic;
                    [yL, yU] = estimateNetworkOutputBounds(netMW, XL, XU);
                    yL = extractdata(yL); yU = extractdata(yU);
                    yL = yL(:); yU = yU(:);
                    status = aivl_vnnlib_status_local(yL, yU, property.prop);
                    tsec   = toc(t);
                otherwise
                    error('verify_image_vnnlib: aivl alg %s not supported', alg);
            end
        otherwise
            error('verify_image_vnnlib: tool %s not supported', tool);
    end
end

% =========================================================================
% Status helpers — kept local so this driver stays self-contained while
% the heavier tool_utils + rebuild_for_aivl live in utils/.
% =========================================================================

function status = nnv_vnnlib_status_local(R, prop)
    if iscell(prop), prop = prop{1}; end
    halves = struct('G',{},'g',{});
    if isstruct(prop) && isscalar(prop) && isfield(prop,'Hg')
        Hg = prop.Hg;
        for k = 1:numel(Hg)
            halves(end+1) = struct('G', Hg(k).G, 'g', Hg(k).g); %#ok<AGROW>
        end
    elseif isstruct(prop)
        for k = 1:numel(prop)
            Hg = prop(k).Hg;
            for j = 1:numel(Hg)
                halves(end+1) = struct('G', Hg(j).G, 'g', Hg(j).g); %#ok<AGROW>
            end
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

function status = aivl_evalbounds_status_local(lb, ub, output)
% Status from AIVL estimateNetworkOutputBounds against a load_vnnlib_matlab
% output spec. Two output formats are possible depending on the benchmark:
%   (a) cell-of-string-pairs: output{i}={"all(lb >= ...)", "all(ub <= ...)"}
%       (rl, acas, tllverify — CAV'23 convention; eval'd against lb/ub).
%   (b) struct with .G and .g: half-space form G*y <= g
%       (collins_rul, oval21 — VNNLIB direct).
    lb = extractdata(lb); ub = extractdata(ub); %#ok<NASGU>
    if iscell(output)
        % v1 path (mirrors aivl_evalbounds_status in
        % ToolComparison/acas_rl_tll/run_acas_rl_tll.m:671).
        n = numel(output);
        r = ones(n, 1);
        for i = 1:n
            r(i) = eval(output{i}{1});
            if ~r(i)
                r(i) = eval(output{i}{2});
                if ~r(i), r(i) = 2; else, r(i) = 0; break; end
            end
        end
        if all(r == 1),      status = "violated";
        elseif any(r == 0),  status = "verified";
        else,                status = "unknown";
        end
    elseif isstruct(output) && isfield(output, 'G') && isfield(output, 'g')
        % Half-space form: verified iff G*y <= g is impossible over [lb,ub].
        lb = lb(:); ub = ub(:);
        G = output.G; g = output.g;
        nrows = size(G,1);
        verified = true;
        violated_all = true;
        for k = 1:nrows
            gk = G(k,:);
            maxv = sum(max(gk.*lb, gk.*ub));
            minv = sum(min(gk.*lb, gk.*ub));
            if minv <= g(k), verified = false;     end
            if maxv >  g(k), violated_all = false; end
        end
        if verified,         status = "verified";
        elseif violated_all, status = "violated";
        else,                status = "unknown";
        end
    else
        status = "unknown";
    end
end

function status = aivl_vnnlib_status_local(yL, yU, prop)
% Evaluate an OR-of-halfspaces VNNLIB property against AIVL output bounds.
% The property encodes counterexample (unsafe) conditions:
%   status="violated" iff at least one half-space (G*y <= g) is forced
%                     over [yL,yU] (definite counterexample)
%   status="verified" iff EVERY halfspace is infeasible over [yL,yU]
%   otherwise         unknown
%
% Half-space feasibility test using interval arithmetic on G*y with the
% Gpos/Gneg split (standard interval bound propagation):
%   max(G*y - g) <= 0  ->  forced  (definite counterexample)
%   min(G*y - g) >  0  ->  infeasible (no counterexample for this disjunct)
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

    if isempty(halves), status = "unknown"; return; end

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

function ensure_aivl_on_path_local()
    persistent done
    if ~isempty(done) && done, return; end
    matches = dir(fullfile(userpath(), 'SupportPackages', 'R*', 'toolbox', 'nnet', 'supportpackages', 'aivnv'));
    if ~isempty(matches)
        addpath(matches(1).folder, '-begin');
    end
    done = true;
end

% =========================================================================
% Utilities
% =========================================================================

function s = strip_dot_slash(s)
    s = string(s);
    if startsWith(s, "./"), s = extractAfter(s, 2); end
end

function decompress_gz_dir(d, pat)
    if ~isfolder(d), return; end
    gzs = dir(fullfile(d, pat));
    for k = 1:numel(gzs)
        gzPath = fullfile(gzs(k).folder, gzs(k).name);
        plain  = regexprep(gzPath, '\.gz$', '');
        if ~isfile(plain), gunzip(gzPath); end
    end
end

function [status, tsec] = run_with_timeout(fn, timeout_s)
%RUN_WITH_TIMEOUT Process-pool parfeval with wall-clock timeout enforcement.
    persistent havePar
    if isempty(havePar)
        havePar = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    end

    if havePar && timeout_s > 0 && isfinite(timeout_s)
        pool = gcp('nocreate');
        if isempty(pool)
            % evalc captures parpool's "Starting parallel pool..." /
            % "Connected to parallel pool with N workers." chatter so it
            % doesn't flood the smoke log.
            evalc('parpool(''local'');');
            pool = gcp('nocreate');
        end
        F  = parfeval(pool, fn, 2);
        ok = wait(F, 'finished', timeout_s);
        if ~ok
            cancel(F);
            status = "timeout";
            tsec   = timeout_s;
            return;
        end
        try
            [status, tsec] = fetchOutputs(F);
        catch ME
            status = "error";
            tsec   = NaN;
            fprintf(2, '    run_with_timeout fetchOutputs error:\n%s\n', getReport(ME, 'extended'));
        end
    else
        t0 = tic;
        try
            [status, tsec] = fn();
        catch ME
            status = "error";
            tsec   = toc(t0);
            fprintf(2, '    run_with_timeout inline error:\n%s\n', getReport(ME, 'extended'));
        end
    end
end
