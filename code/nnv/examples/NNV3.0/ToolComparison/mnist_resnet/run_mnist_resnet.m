function run_mnist_resnet(varargin)
%RUN_MNIST_RESNET ToolComparison ResNet half: argmax robustness on MNIST-ResNet-8.
%
%   First head-to-head against AIVL's verifyNetworkRobustness on a residual
%   network -- enabled by R2024b's `additionLayer` support.
%
%   Models:
%       mnist_resnet8          ~8 conv + 3 residual blocks (additionLayer)
%
%   Tools:
%       'nnv'          -> NNV 3.0 reachability
%       'mw_deeppoly'  -> verifyNetworkRobustness (DeepPoly, R2024b+)
%
%   NNV reach-method grid (CAV'23-style, area-based for image classification):
%       'approx-star'              fastest
%       'relax-star-area-25/50/75/100'    relaxFactor 0.25/0.5/0.75/1.0
%       'exact-star'               (intractable on ResNet; not in defaults)
%       'relax-star-50'            legacy alias for 'relax-star-area-50'
%
%   Perturbation: L_inf box, eps in {1/255, 2/255, 4/255}, ~50 test points
%   per model.
%
%   Examples:
%     run_mnist_resnet()
%     run_mnist_resnet('models',{'mnist_resnet8'})
%     run_mnist_resnet('epsilons', [1/255, 4/255])
%     run_mnist_resnet('numPoints', 25)
%     run_mnist_resnet('algorithms', {'approx-star','relax-star-area-50'})

    p = inputParser;
    addParameter(p, 'models',     {'mnist_resnet8'});
    addParameter(p, 'tools',      {'nnv','mw_deeppoly'});
    addParameter(p, 'algorithms', {});
    addParameter(p, 'epsilons',   [1/255, 2/255, 4/255]);
    addParameter(p, 'numPoints',  50);
    addParameter(p, 'timeout',    600);
    addParameter(p, 'rerun',      {});
    addParameter(p, 'resultsDir', fullfile(fileparts(mfilename('fullpath')), 'results'));
    parse(p, varargin{:});
    opts = p.Results;
    if ~isfolder(opts.resultsDir), mkdir(opts.resultsDir); end

    u = tool_utils();
    fprintf("mnist_resnet: models    = {%s}\n", strjoin(opts.models, ", "));
    fprintf("mnist_resnet: epsilons  = [%s]\n", num2str(opts.epsilons));
    fprintf("mnist_resnet: points    = %d per (model, eps)\n", opts.numPoints);
    fprintf("mnist_resnet: timeout   = %g s\n\n", opts.timeout);

    for m = 1:numel(opts.models)
        model   = opts.models{m};
        matFile = fullfile(opts.resultsDir, sprintf("expC_%s.mat", model));
        if ~isempty(opts.rerun)
            n = u.purge_status(matFile, opts.rerun);
            if n > 0, fprintf("  purged %d rows with status in {%s}\n", n, strjoin(opts.rerun, ', ')); end
        end
        [net_mw, net_nnv, X, y] = load_model_and_data(model, opts.numPoints);
        fprintf("=== %s : %d test points -> %s ===\n", model, numel(y), matFile);

        for ei = 1:numel(opts.epsilons)
            eps = opts.epsilons(ei);
            for pi = 1:numel(y)
                instance_id = sprintf("%s|eps=%.4f|pt=%d", model, eps, pi);
                x0   = X(:,:,:,pi);
                ytrue = y(pi);

                for t = 1:numel(opts.tools)
                    tool = opts.tools{t};
                    algs = algorithms_for(tool, model, opts.algorithms);
                    for k = 1:numel(algs)
                        alg = algs{k};
                        if u.has_instance(matFile, tool, model, instance_id, alg), continue; end
                        fprintf("  [%s %-13s %-22s] %s ... ", model, tool, alg, instance_id);
                        [status, tsec] = run_one(@() verifyResNet(net_mw, net_nnv, x0, ytrue, eps, ...
                                                                  tool, alg, opts.timeout), opts.timeout);
                        row = u.new_row(tool, model, instance_id, status, tsec, alg, opts.timeout);
                        u.append_to_mat(matFile, row);
                        fprintf("%s (%s s)\n", status, num2str(tsec));
                    end
                end
            end
        end
    end

    fprintf("\nmnist_resnet: done. Regenerate tables with:  make_mnist_resnet_table\n");
end

% -------------------------------------------------------------------------

function algs = algorithms_for(tool, model, override) %#ok<INUSD>
% Default algorithm set for the ResNet head-to-head. Matches the NNV 2.0
% paper methodology for residual networks: a single relax-star-area-50
% NNV configuration vs AIVL's deep-poly. Other relax/approx algorithms
% remain accepted via the 'algorithms' override but are not part of the
% reproducible Table C grid.
    switch tool
        case 'nnv',         algs = {'relax-star-area-50'};
        case 'mw_deeppoly', algs = {'deep-poly'};
        case 'mw_abc',      algs = {'alpha-beta-crown'};
        otherwise, error("Unknown tool: %s", tool);
    end
    if nargin >= 3 && ~isempty(override)
        algs = algs(ismember(algs, override));
    end
end

function [net_mw, net_nnv, X, y] = load_model_and_data(model, numPoints)
% Load the trained dlnetwork + a labeled test array from models/.
    mdir = models_dir();
    netFile = fullfile(mdir, sprintf("%s.mat", model));
    if ~exist(netFile, 'file')
        error("load_model_and_data: missing trained model %s", netFile);
    end
    S = load(netFile);
    if ~isfield(S, 'net')
        error("load_model_and_data: %s does not contain variable 'net'", netFile);
    end
    net_mw = S.net;
    % Strip softmax: AIVL's verifyNetworkRobustness rejects SoftmaxLayer.
    smIdx = find(arrayfun(@(L) isa(L,'nnet.cnn.layer.SoftmaxLayer'), net_mw.Layers), 1);
    if ~isempty(smIdx)
        net_mw = removeLayers(net_mw, net_mw.Layers(smIdx).Name);
        net_mw = initialize(net_mw);
    end
    net_nnv = matlab2nnv(net_mw);

    tsFile = fullfile(mdir, sprintf("%s_testset.mat", model));
    if ~exist(tsFile, 'file')
        error("load_model_and_data: missing testset %s", tsFile);
    end
    T = load(tsFile);
    Xall = T.Xtest;
    Yall = T.Ytest;
    nAll = size(Xall, 4);
    pick = 1:min(numPoints, nAll);
    X = Xall(:,:,:,pick);
    y = Yall(pick);
end

function [status, tsec] = verifyResNet(net_mw, net_nnv, x0, ytrue, eps, tool, alg, timeout) %#ok<INUSD>
% Verify L_inf robustness around x0 with radius eps. Returns canonical status.
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
            reachOpt = nnv_reach_opt(alg);
            R = net_nnv.reach(IS, reachOpt);
            tsec = toc(t);
            status = nnv_robust_status(R, ytrue);

        case 'mw_deeppoly'
            % Worker-side AIVL path setup (parfeval workers don't auto-run startup.m).
            sps = dir(fullfile(userpath(), 'SupportPackages', 'R*', 'toolbox', 'nnet', 'supportpackages', 'aivnv'));
            for kk = 1:numel(sps)
                addpath(fullfile(sps(kk).folder, sps(kk).name), '-begin');
            end
            XL = dlarray(XLow, "SSCB");
            XU = dlarray(XUp,  "SSCB");
            ytrueIdx = double(ytrue);
            t = tic;
            r = verifyNetworkRobustness(net_mw, XL, XU, ytrueIdx);
            tsec = toc(t);
            status = mw_result_to_status(r);

        case 'mw_abc'
            error('verifyResNet: mw_abc requires the R2026a alpha-beta-CROWN bridge.');

        otherwise
            error("verifyResNet: unknown tool %s", tool);
    end
end

function reachOpt = nnv_reach_opt(alg)
% NNV reachability options keyed by algorithm. Image-classification half
% uses 'relax-star-area' (area-based heuristic for ImageStars); the FC half
% (acas_rl_tll) uses 'relax-star-range'. The legacy alias 'relax-star-50'
% maps to canonical 'relax-star-area-50' for back-compat with earlier .mat.
    reachOpt = struct;
    reachOpt.reachOption = "single";
    reachOpt.numCores    = 1;
    reachOpt.relaxFactor = 0;
    switch alg
        case 'approx-star'
            reachOpt.reachMethod = 'approx-star';
        case {'relax-star-area-25'}
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 0.25;
        case {'relax-star-area-50','relax-star-50'}
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 0.5;
        case {'relax-star-area-75'}
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 0.75;
        case {'relax-star-area-100'}
            reachOpt.reachMethod = 'relax-star-area';
            reachOpt.relaxFactor = 1.0;
        case 'exact-star'
            reachOpt.reachMethod = 'exact-star';
            reachOpt.reachOption = "parallel";
            reachOpt.numCores    = min(8, maxNumCompThreads);
        otherwise
            error("nnv_reach_opt: unknown algorithm %s", alg);
    end
end

function status = nnv_robust_status(R, ytrue)
% argmax(R) == ytrue across every output set?
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

function status = mw_result_to_status(r)
    s = string(r);
    switch s
        case "verified", status = "verified";
        case "violated", status = "violated";
        otherwise,       status = "unknown";
    end
end

function [status, tsec] = run_one(fn, timeout)
% Process-pool parfeval for wall-clock timeout.
    persistent havePar
    if isempty(havePar)
        havePar = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
    end
    status = "error"; tsec = NaN; %#ok<NASGU>
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

function d = models_dir()
    d = fullfile(fileparts(mfilename('fullpath')), 'models');
end
