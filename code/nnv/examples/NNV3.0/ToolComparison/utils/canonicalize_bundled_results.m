function canonicalize_bundled_results()
%CANONICALIZE_BUNDLED_RESULTS Rewrite bundled .mat results to the current schema.
%
%   Migrations applied (idempotent — safe to re-run):
%     1. Algorithm string `relax-star-50` -> `relax-star-range-50` (acas_rl_tll)
%        or `relax-star-area-50` (mnist_resnet).
%     2. Tool string `mw_estimate` -> tool=`aivl`, algorithm=`estimate-bounds`
%        (the algorithm value is preserved if it was already `estimate-bounds`).
%     3. Tool string `mw_deeppoly` -> tool=`aivl`, algorithm=`deep-poly`.
%     4. Tool string `mw_abc`      -> tool=`aivl`, algorithm=`alpha-beta-crown`.
%     5. Dedupes (tool, benchmark, instance_id, algorithm) keeping the first
%        occurrence (bundled paper-grade rows precede any newer smoke rows).

    here = fileparts(fileparts(mfilename('fullpath')));
    addpath(here);
    u = tool_utils();

    acasFiles = {'results_acas_p3.mat','results_acas_p4.mat', ...
                 'results_rl.mat','results_oval21.mat','results_collins_rul.mat', ...
                 'results_tllverify.mat'};
    for k = 1:numel(acasFiles)
        f = fullfile(here, 'acas_rl_tll', 'results', acasFiles{k});
        if isfile(f), rewrite(f, 'relax-star-range-50', u); end
    end

    resFile = fullfile(here, 'mnist_resnet', 'results', 'expC_mnist_resnet8.mat');
    if isfile(resFile), rewrite(resFile, 'relax-star-area-50', u); end
end

function rewrite(matFile, canonicalRelaxName, u)
    R = u.load_results(matFile);
    n0 = height(R);

    % (1) Algorithm-rename: legacy "relax-star-50" -> canonical name.
    legacy = R.algorithm == "relax-star-50";
    nAlg = sum(legacy);
    R.algorithm(legacy) = string(canonicalRelaxName);

    % (2-4) Tool-rename: collapse mw_* to aivl, set algorithm if missing.
    %   mw_estimate  -> aivl + algorithm 'estimate-bounds'  (only set if alg string is empty/legacy)
    %   mw_deeppoly  -> aivl + algorithm 'deep-poly'
    %   mw_abc       -> aivl + algorithm 'alpha-beta-crown'
    [R, nT_est] = remap_tool(R, "mw_estimate",  "aivl", "estimate-bounds");
    [R, nT_dp]  = remap_tool(R, "mw_deeppoly",  "aivl", "deep-poly");
    [R, nT_abc] = remap_tool(R, "mw_abc",       "aivl", "alpha-beta-crown");

    % (5) Dedupe by (tool, benchmark, instance_id, algorithm).
    [~, ia] = unique(R(:, {'tool','benchmark','instance_id','algorithm'}), ...
                     'rows', 'stable');
    nDup = height(R) - numel(ia);
    R = R(ia, :);

    results = R; %#ok<NASGU>
    save(matFile, 'results');

    fprintf(['[canonicalize] %s: rows %d -> %d  (alg-renames %d, tool-renames mw_estimate=%d ' ...
             'mw_deeppoly=%d mw_abc=%d, dropped %d dup)\n'], ...
            matFile, n0, height(R), nAlg, nT_est, nT_dp, nT_abc, nDup);
end

function [R, n] = remap_tool(R, oldTool, newTool, defaultAlg)
% Rewrite rows where R.tool == oldTool to (newTool, defaultAlg) — but only
% override the algorithm column if it's empty or matches the legacy alg.
    sel = R.tool == oldTool;
    n = sum(sel);
    if n == 0, return; end
    R.tool(sel) = newTool;
    % Preserve alg if it's already set to something meaningful (e.g. an
    % override-driven row); otherwise set to defaultAlg.
    needsAlg = sel & (R.algorithm == "" | ismissing(R.algorithm));
    R.algorithm(needsAlg) = defaultAlg;
end
