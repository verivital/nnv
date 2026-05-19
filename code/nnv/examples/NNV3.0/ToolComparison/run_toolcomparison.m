function run_toolcomparison(varargin)
%RUN_TOOLCOMPARISON  Top-level entry for the NNV-vs-AIVL artifact-evaluation
%comparison.
%
%   run_toolcomparison()                              -- default mode, all tools, all algorithms
%   run_toolcomparison('mode','smoke')                -- ~10-min smoke
%   run_toolcomparison('mode','default')              -- ~3-5 h full grid
%   run_toolcomparison('mode','full')                 -- alias for 'default' (paper convention)
%   run_toolcomparison('mode','smoke','tools',{'nnv'}) -- NNV-only (skip AIVL)
%   run_toolcomparison('mode','default','algorithms', {'approx-star','deep-poly'}) -- restrict
%
%   If 'mode' is not passed, honors the TOOLCOMPARISON_MODE env var
%   (smoke|default|full); defaults to 'default'. This matches the
%   contract expected by ../run_all.sh (the NNV3.0 orchestrator).
%
%   On entry: addpath's the local utils, runs the vnnlib half then the
%   argmax half, and renders the consolidated table to
%   tables/out/table_main.{tex,txt}.

    % Resolve env-var fallback for mode before parsing args.
    envMode = getenv('TOOLCOMPARISON_MODE');
    if isempty(envMode), envMode = 'default'; end

    p = inputParser;
    p.addParameter('mode',        envMode,    @(x) ischar(x) || isstring(x));
    p.addParameter('tools',       {'nnv','aivl'}, @iscell);
    p.addParameter('algorithms',  {},        @iscell);
    p.addParameter('halves',      {'vnnlib','argmax'}, @iscell);
    p.addParameter('benchmarks',  {},        @iscell);  % empty = all in scope; else filter to listed
    p.parse(varargin{:});
    opts = p.Results;
    opts.mode = char(opts.mode);
    % 'full' is the legacy alias for 'default' used by the NNV3.0 paper /
    % top-level orchestrator (legacy run_toolcomparison.m supported 'full').
    if strcmp(opts.mode, 'full'), opts.mode = 'default'; end
    if ~any(strcmp(opts.mode, {'smoke','default'}))
        error('run_toolcomparison: mode must be ''smoke''|''default''|''full'' (got ''%s'').', opts.mode);
    end

    here = fileparts(mfilename('fullpath'));

    fprintf('=== ToolComparison (mode=%s) ===\n', opts.mode);

    % addpath local utils (tool_utils, rebuild_for_aivl, parse_argmax_vnnlib).
    addpath(fullfile(here, 'utils'));
    addpath_shared();
    addpath(fullfile(here, 'drivers'));

    % NNV must already be installed and on the path.
    if isempty(which('matlab2nnv'))
        error(['run_toolcomparison: NNV not on MATLAB path. ', ...
               'Run code/nnv/install.m or addpath(genpath(''code/nnv/'')) first.']);
    end

    % 4. tool_utils handle (single canonical schema).
    u = tool_utils();

    t0 = tic;
    if ismember('vnnlib', opts.halves)
        fprintf('\n--- vnnlib half ---\n');
        run_vnnlib_half(opts, u);
    end
    if ismember('argmax', opts.halves)
        fprintf('\n--- argmax half ---\n');
        run_argmax_half(opts, u);
    end
    walltime = toc(t0);

    % 5. Render consolidated table.
    fprintf('\n--- table render ---\n');
    addpath(fullfile(here, 'tables'));
    try
        make_table_main();
    catch ME
        warning('run_toolcomparison:table_render_failed', 'make_table_main threw: %s', ME.message);
    end

    fprintf('\n=== run_toolcomparison done (mode=%s, wall=%.1f min) ===\n', opts.mode, walltime/60);
end
