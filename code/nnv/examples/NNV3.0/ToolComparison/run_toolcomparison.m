function run_toolcomparison(varargin)
%RUN_TOOLCOMPARISON Top-level orchestrator for the NNV vs MathWorks AIVL comparison.
%
%   Dispatches to the two halves of the comparison and renders the paper
%   tables. Designed to be invoked from NNV3.0/run_all.sh as a single MATLAB
%   batch command, alongside FairNNV / GNNV / ProbVer / VideoStar / ModelStar.
%
%   Halves:
%       acas_rl_tll/run_acas_rl_tll.m    -- ACAS Xu p3/p4, RL controllers,
%                                           OVAL21 (VNNCOMP'21), Collins RUL
%                                           (VNNCOMP'22). NNV vs AIVL
%                                           estimateNetworkOutputBounds.
%       mnist_resnet/run_mnist_resnet.m  -- MNIST-ResNet-8 argmax robustness
%                                           (NNV vs AIVL verifyNetworkRobustness
%                                           with DeepPoly).
%
%   After both halves finish, renders the ATVA 2026 paper's Tables 5, 6, 7:
%       tables/out/table_A.{tex,txt} (paper Tables 5+6: FC and CNN VNNLIB)
%       tables/out/table_C.{tex,txt} (paper Table 7: MNIST-ResNet-8)
%       tables/out/sanity_report.txt (CAV'23 exact-star cross-check)
%
%   Examples:
%     run_toolcomparison()                 % full grid (~3-5 h)
%     run_toolcomparison('mode','smoke')   % ~12 min smoke
%     run_toolcomparison('halves',{'acas_rl_tll'})  % skip ResNet half
%
%   Requires the AI Verification Library (AIVL) Support Package (the MW-side calls
%   `verifyNetworkRobustness` and `estimateNetworkOutputBounds`). NNV-only
%   runs work without it; pass 'tools',{'nnv'} to skip MW comparisons.

    % Default mode: env var TOOLCOMPARISON_MODE if set ('smoke'|'full'),
    % otherwise 'full'. run_all.sh sets it to 'smoke' so the full suite stays
    % under ~30 min total.
    defaultMode = getenv('TOOLCOMPARISON_MODE');
    if isempty(defaultMode), defaultMode = 'full'; end

    p = inputParser;
    addParameter(p, 'mode',   defaultMode);                   % 'full' | 'smoke'
    addParameter(p, 'halves', {'acas_rl_tll','mnist_resnet'});
    addParameter(p, 'tools',  {});                            % empty = each half's defaults
    addParameter(p, 'algorithms', {});                        % empty = each half's defaults
    addParameter(p, 'rerun', {});
    parse(p, varargin{:});
    opts = p.Results;

    here = fileparts(mfilename('fullpath'));
    addpath(genpath(here));                 % runner + utils/, acas_rl_tll/, mnist_resnet/, tables/

    fprintf("\n========== ToolComparison run start ==========\n");
    fprintf("mode:    %s\n", opts.mode);
    fprintf("halves:  {%s}\n\n", strjoin(opts.halves, ", "));

    if any(strcmp(opts.halves, 'acas_rl_tll'))
        switch opts.mode
            case 'smoke'
                run_acas_rl_tll('benchmarks',{'acas_p3'}, 'tools',{'nnv'}, ...
                                'algorithms',{'approx-star','relax-star-range-50'}, ...
                                'numNets',5, 'timeout',60);
            case 'full'
                if isempty(opts.tools) && isempty(opts.algorithms)
                    run_acas_rl_tll('rerun', opts.rerun);
                else
                    args = {};
                    if ~isempty(opts.tools),      args = [args, {'tools', opts.tools}]; end
                    if ~isempty(opts.algorithms), args = [args, {'algorithms', opts.algorithms}]; end
                    if ~isempty(opts.rerun),      args = [args, {'rerun', opts.rerun}]; end
                    run_acas_rl_tll(args{:});
                end
            otherwise, error("Unknown mode: %s", opts.mode);
        end
    end

    if any(strcmp(opts.halves, 'mnist_resnet'))
        switch opts.mode
            case 'smoke'
                run_mnist_resnet('models',{'mnist_resnet8'}, 'tools',{'nnv'}, ...
                                 'algorithms',{'approx-star','relax-star-area-50'}, ...
                                 'numPoints',5, 'epsilons',[1/255], 'timeout',60);
            case 'full'
                if isempty(opts.tools) && isempty(opts.algorithms)
                    run_mnist_resnet('rerun', opts.rerun);
                else
                    args = {};
                    if ~isempty(opts.tools),      args = [args, {'tools', opts.tools}]; end
                    if ~isempty(opts.algorithms), args = [args, {'algorithms', opts.algorithms}]; end
                    if ~isempty(opts.rerun),      args = [args, {'rerun', opts.rerun}]; end
                    run_mnist_resnet(args{:});
                end
        end
    end

    % Render tables (cheap; always do it, falls back gracefully if a half had
    % no rows).
    try
        make_acas_rl_tll_table();
    catch ME
        fprintf(2, "make_acas_rl_tll_table failed: %s\n", ME.message);
    end
    try
        make_mnist_resnet_table();
    catch ME
        fprintf(2, "make_mnist_resnet_table failed: %s\n", ME.message);
    end
    try
        export_csv();
    catch ME
        fprintf(2, "export_csv failed: %s\n", ME.message);
    end

    fprintf("\n========== ToolComparison run done ==========\n");
end
