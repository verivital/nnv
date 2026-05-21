function run_argmax_half(opts, u)
%RUN_ARGMAX_HALF  driver for mnist_resnet8 (argmax L∞ robustness).
%
%   Single benchmark, single NNV algorithm (relax-star-area-50), single
%   AIVL algorithm (deep-poly via verifyNetworkRobustness). Argmax-form,
%   so we read ytrue from the bundled testset (not from a VNNLIB file).

    here = fileparts(mfilename('fullpath'));
    tc_root = fileparts(here);
    results_dir = fullfile(tc_root, 'results');
    if ~isfolder(results_dir), mkdir(results_dir); end

    bench = 'mnist_resnet8';
    % Optional benchmarks filter: if a non-empty filter was provided and
    % does not include this benchmark, skip silently.
    if isfield(opts, 'benchmarks') && ~isempty(opts.benchmarks) ...
            && ~ismember(bench, opts.benchmarks)
        return;
    end
    assetDir = fullfile(tc_root, 'benchmarks', bench);
    matFile  = fullfile(results_dir, sprintf('%s.mat', bench));
    netFile  = fullfile(assetDir, 'models', 'mnist_resnet8.mat');
    tsFile   = fullfile(assetDir, 'models', 'mnist_resnet8_testset.mat');

    if ~isfile(netFile) || ~isfile(tsFile)
        fprintf('[%s] SKIP — missing bundled model or testset\n', bench);
        return;
    end

    [net_mw, net_nnv, X, y] = load_mnist_resnet8(netFile, tsFile);

    % Mode-dependent scope.
    if strcmp(opts.mode, 'smoke')
        nImages   = 1;
        epsilons  = 1/255;
    else
        % Reduced from 50 imgs to 25 imgs; preserves 4-eps scaling story.
        % 25 × 4 = 100 instances × 2 algos = 200 calls; ~20-30 min total.
        nImages   = min(25, size(X, 4));
        epsilons  = [1/255, 2/255, 4/255, 8/255];
    end

    per_instance_timeout = 300;       % 5-min cap per (image, eps) call

    for iEps = 1:numel(epsilons)
        eps = epsilons(iEps);
        for iImg = 1:nImages
            instance_id = sprintf('img%d|eps=%s', iImg, num2str(eps,'%.5g'));
            % No parpool restart between argmax instances: mnist_resnet8's
            % relax-star-area-50 and AIVL deep-poly are both single-threaded;
            % the per-instance restart added ~30 s overhead per instance with
            % no isolation benefit.

            for t = 1:numel(opts.tools)
                tool = opts.tools{t};
                algs = algs_for(tool, opts.mode, opts.algorithms);
                for k = 1:numel(algs)
                    alg = algs{k};
                    if u.has_instance(matFile, tool, bench, instance_id, alg)
                        continue;
                    end
                    fprintf('  [%s %-4s %-22s] %s ... ', bench, tool, alg, instance_id);
                    [status, tsec] = run_with_timeout( ...
                        @() verify_resnet8(net_mw, net_nnv, X(:,:,:,iImg), y(iImg), eps, tool, alg), ...
                        per_instance_timeout);
                    row = u.new_row(tool, bench, instance_id, status, tsec, alg, per_instance_timeout);
                    u.append_to_mat(matFile, row);
                    fprintf('%s (%s s)\n', status, u.format_time(tsec, per_instance_timeout));
                end
            end
        end
    end
end

% =========================================================================

function algs = algs_for(tool, mode, override)
    switch tool
        case 'nnv',  full = {'relax-star-area-50'};
        case 'aivl', full = {'deep-poly'};
        otherwise,   error('algs_for: unknown tool %s', tool);
    end
    % Smoke and default both run their full algorithm list (only 1 alg per tool
    % on mnist_resnet8 either way). AIVL inclusion is gated upstream by
    % run_toolcomparison.m's probe; if AIVL is missing from path, opts.tools is
    % filtered there and this function isn't invoked for 'aivl'.
    algs = full;
    if nargin >= 3 && ~isempty(override)
        algs = algs(ismember(algs, override));
    end
end

function [net_mw, net_nnv, X, y] = load_mnist_resnet8(netFile, tsFile)
    S = load(netFile);
    if ~isfield(S, 'net')
        error('mnist_resnet8: %s has no variable ''net''', netFile);
    end
    net_mw = S.net;
    % Strip SoftmaxLayer — AIVL verifyNetworkRobustness rejects it.
    smIdx = find(arrayfun(@(L) isa(L,'nnet.cnn.layer.SoftmaxLayer'), net_mw.Layers), 1);
    if ~isempty(smIdx)
        net_mw = removeLayers(net_mw, net_mw.Layers(smIdx).Name);
        net_mw = initialize(net_mw);
    end
    net_nnv = matlab2nnv(net_mw);

    T = load(tsFile);
    X = T.Xtest;
    y = T.Ytest;
end

function [status, tsec] = verify_resnet8(net_mw, net_nnv, x0, ytrue, eps, tool, alg)
    status = "error"; tsec = NaN; %#ok<NASGU>
    x0d = single(x0);
    if ndims(x0d) == 2 %#ok<ISMAT>
        x0d = reshape(x0d, size(x0d,1), size(x0d,2), 1);
    end
    XLow = max(x0d - single(eps), 0);
    XUp  = min(x0d + single(eps), 255);

    switch tool
        case 'nnv'
            t = tic;
            IS = ImageStar(XLow, XUp);
            reachOpt = reach_opt_for(alg);
            R = net_nnv.reach(IS, reachOpt);
            tsec = toc(t);
            status = nnv_robust_status_local(R, ytrue);
        case 'aivl'
            matches = dir(fullfile(userpath(), 'SupportPackages', 'R*', 'toolbox', 'nnet', 'supportpackages', 'aivnv'));
            if ~isempty(matches)
                addpath(matches(1).folder, '-begin');
            end
            switch alg
                case 'deep-poly'
                    XL = dlarray(XLow, "SSCB");
                    XU = dlarray(XUp,  "SSCB");
                    ytrueIdx = double(ytrue);
                    t = tic;
                    r = verifyNetworkRobustness(net_mw, XL, XU, ytrueIdx);
                    tsec = toc(t);
                    status = aivl_categorical_to_status(r);
                otherwise
                    error('verify_resnet8: aivl alg %s not supported', alg);
            end
        otherwise
            error('verify_resnet8: unknown tool %s', tool);
    end
end

function status = nnv_robust_status_local(R, ytrue)
    yIdx = double(ytrue);
    status = "verified";
    for k = 1:numel(R)
        Set = R(k);
        if isa(Set, 'ImageStar'), Set = Set.toStar; end
        nOut = Set.dim;
        for j = 1:nOut
            if j == yIdx, continue; end
            H = zeros(1, nOut); H(yIdx) = 1; H(j) = -1;
            S = Set.intersectHalfSpace(H, 0);
            if ~isempty(S)
                if isempty(Set.intersectHalfSpace(-H, 0))
                    status = "violated"; return;
                else
                    status = "unknown"; return;
                end
            end
        end
    end
end

function status = aivl_categorical_to_status(r)
    s = string(r);
    switch s
        case "verified", status = "verified";
        case "violated", status = "violated";
        otherwise,       status = "unknown";
    end
end

function [status, tsec] = run_with_timeout(fn, timeout_s) %#ok<INUSD>
%RUN_WITH_TIMEOUT Run fn() inline (no parfeval). For the argmax half
%(mnist_resnet8) the per-instance work is single-threaded and short
%(<15 s typical), so the parfeval pool overhead (~30 s spawn per
%instance) dominated wall time. Running inline trades per-instance
%timeout enforcement for ~10x faster end-to-end wall.
    t0 = tic;
    try
        [status, tsec] = fn();
    catch ME
        status = "error";
        tsec   = toc(t0);
        fprintf(2, '    run_with_timeout inline error:\n%s\n', getReport(ME, 'extended'));
    end
end
