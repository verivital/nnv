function canonicalize_bundled_results()
%CANONICALIZE_BUNDLED_RESULTS One-shot rewrite of bundled .mat results so the
%   legacy `relax-star-50` algorithm string maps to the canonical
%   `relax-star-range-50` (acas_rl_tll) or `relax-star-area-50` (mnist_resnet).
%   Dedupes (tool, benchmark, instance_id, algorithm) keeping the first
%   occurrence (bundled paper-grade rows precede any newer smoke-run rows).

    here = fileparts(fileparts(mfilename('fullpath')));
    addpath(here);
    u = tool_utils();

    acasFiles = {'results_acas_p3.mat','results_acas_p4.mat', ...
                 'results_rl.mat','results_tllverify.mat'};
    for k = 1:numel(acasFiles)
        f = fullfile(here, 'acas_rl_tll', 'results', acasFiles{k});
        if isfile(f), rewrite(f, 'relax-star-range-50', u); end
    end

    resFile = fullfile(here, 'mnist_resnet', 'results', 'expC_mnist_resnet8.mat');
    if isfile(resFile), rewrite(resFile, 'relax-star-area-50', u); end
end

function rewrite(matFile, canonicalName, u)
    R = u.load_results(matFile);
    n0 = height(R);
    legacy = R.algorithm == "relax-star-50";
    nLegacy = sum(legacy);
    R.algorithm(legacy) = string(canonicalName);

    [~, ia] = unique(R(:, {'tool','benchmark','instance_id','algorithm'}), ...
                     'rows', 'stable');
    nDup = height(R) - numel(ia);
    R = R(ia, :);

    results = R; %#ok<NASGU>
    save(matFile, 'results');

    fprintf("[canonicalize] %s: rows %d -> %d  (renamed %d legacy, dropped %d dup)\n", ...
            matFile, n0, height(R), nLegacy, nDup);
end
