function aivl_install()
%AIVL_INSTALL  Install the AI Verification Library (AIVL) Support Package.
%   Resolves the tarball at <ToolComparison>/atva26-aivl.tar.gz, sets the
%   TOOLCOMPARISON_AIVL_TARBALL env var, and invokes the local
%   toolbox_install.m which extracts to userpath + wires startup.m.
%
%   Idempotent.

    here = fileparts(mfilename('fullpath'));
    tarball = fullfile(here, 'atva26-aivl.tar.gz');

    if ~isfile(tarball)
        error('ToolComparison:aivl_install:missing_tarball', ...
              ['AIVL tarball not found at %s.\n', ...
               'The AIVL Support Package is non-redistributable, so this tarball ', ...
               'is not committed to the repository. Acquisition paths:\n', ...
               '  (1) Install AIVL yourself via MATLAB Home -> Add-Ons -> Get Add-Ons ', ...
               '("AI Verification Library"); aivl_install is not needed in that case.\n', ...
               '  (2) ATVA 2026 AE reviewers: place the tarball from the private link ', ...
               '(in the HotCRP submission cover note) at the path above, then re-run.\n', ...
               '  (3) Code Ocean reviewers: AIVL is pre-installed; aivl_install is not needed.\n', ...
               'See ToolComparison/README.md section "Install AIVL" for details.'], tarball);
    end

    setenv('TOOLCOMPARISON_AIVL_TARBALL', tarball);

    installer = fullfile(here, 'toolbox_install.m');
    if ~isfile(installer)
        error('ToolComparison:aivl_install:missing_installer', ...
              'Expected local installer at %s', installer);
    end

    fprintf('[aivl_install] Delegating to %s (tarball: %s)\n', installer, tarball);
    run(installer);
end
