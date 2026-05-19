function uninstall_patches()
%UNINSTALL_PATCHES  Reverse-apply the two NNV engine patches shipped with
%this artifact. Use to restore the working tree to its pre-install state.

    here = fileparts(mfilename('fullpath'));
    [status, repoRoot] = system('git rev-parse --show-toplevel');
    if status ~= 0
        error('ToolComparison:uninstall_patches:not_a_repo', ...
              'uninstall_patches must be run from within a git checkout of nnv.');
    end
    repoRoot = strtrim(repoRoot);

    patchFiles = {
        '01_FeatureInputLayer.patch'
        '02_NN_start_pool.patch'
    };

    nReverted = 0;
    nSkipped = 0;
    % Reverse order in case any patches have inter-file dependencies.
    for i = numel(patchFiles):-1:1
        pf = fullfile(here, patchFiles{i});
        cmd = sprintf('cd %s && git apply --reverse --check %s 2>&1', esc(repoRoot), esc(pf));
        [rc, ~] = system(cmd);
        if rc ~= 0
            fprintf('[uninstall_patches] SKIP %s (not currently applied)\n', patchFiles{i});
            nSkipped = nSkipped + 1;
            continue;
        end
        cmd = sprintf('cd %s && git apply --reverse %s 2>&1', esc(repoRoot), esc(pf));
        [rc, out] = system(cmd);
        if rc ~= 0
            error('ToolComparison:uninstall_patches:revert_failed', ...
                  'git apply --reverse failed for %s. Output:\n%s', patchFiles{i}, out);
        end
        fprintf('[uninstall_patches] REVERT %s\n', patchFiles{i});
        nReverted = nReverted + 1;
    end

    fprintf('[uninstall_patches] %d reverted, %d already-clean. Done.\n', ...
            nReverted, nSkipped);
end

function s = esc(p)
    s = ['''' strrep(p, '''', '''\''''') ''''];
end
