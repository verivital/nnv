function install_patches()
%INSTALL_PATCHES  Apply the two NNV engine patches shipped with this
%artifact. Idempotent: rerunning is safe.
%
%   01_FeatureInputLayer.patch    — restore pre-Oct-2025 reach_star_single_input
%   02_NN_start_pool.patch        — getCurrentTask() worker bail in start_pool
%
%   For each patch:
%     1. `git apply --reverse --check` to see if already applied. If yes, skip.
%     2. `git apply --check`         to see if it applies cleanly. If no, fail.
%     3. `git apply`                  to apply.
%
%   The repo root is detected via `git rev-parse --show-toplevel`. The user
%   must run this from a clone of the nnv repo with no other staged changes
%   to these three files.

    here = fileparts(mfilename('fullpath'));
    [status, repoRoot] = system('git rev-parse --show-toplevel');
    if status ~= 0
        error('ToolComparison:install_patches:not_a_repo', ...
              'install_patches must be run from within a git checkout of nnv.');
    end
    repoRoot = strtrim(repoRoot);

    patchFiles = {
        '01_FeatureInputLayer.patch'
        '02_NN_start_pool.patch'
    };

    nApplied = 0;
    nSkipped = 0;
    for i = 1:numel(patchFiles)
        pf = fullfile(here, patchFiles{i});
        if ~isfile(pf)
            error('ToolComparison:install_patches:missing', ...
                  'Patch file missing: %s', pf);
        end

        % Already applied?
        cmd = sprintf('cd %s && git apply --reverse --check %s 2>&1', ...
                      esc(repoRoot), esc(pf));
        [rc, ~] = system(cmd);
        if rc == 0
            fprintf('[install_patches] SKIP  %s (already applied)\n', patchFiles{i});
            nSkipped = nSkipped + 1;
            continue;
        end

        % Clean apply check
        cmd = sprintf('cd %s && git apply --check %s 2>&1', ...
                      esc(repoRoot), esc(pf));
        [rc, out] = system(cmd);
        if rc ~= 0
            error('ToolComparison:install_patches:would_conflict', ...
                  ['Patch %s does not apply cleanly. Output:\n%s\n', ...
                   'Likely cause: local edits to the target file. ', ...
                   'Resolve manually or stash your changes first.'], ...
                  patchFiles{i}, out);
        end

        % Apply
        cmd = sprintf('cd %s && git apply %s 2>&1', esc(repoRoot), esc(pf));
        [rc, out] = system(cmd);
        if rc ~= 0
            error('ToolComparison:install_patches:apply_failed', ...
                  'git apply failed for %s. Output:\n%s', patchFiles{i}, out);
        end
        fprintf('[install_patches] APPLY %s\n', patchFiles{i});
        nApplied = nApplied + 1;
    end

    fprintf('[install_patches] %d applied, %d already-applied. Done.\n', ...
            nApplied, nSkipped);
end

function s = esc(p)
    % shell-escape a path with single quotes
    s = ['''' strrep(p, '''', '''\''''') ''''];
end
