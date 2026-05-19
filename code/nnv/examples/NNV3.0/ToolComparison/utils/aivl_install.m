function aivl_install()
%AIVL_INSTALL  Install the AI Verification Library (AIVL) Support Package.
%   Resolves the tarball at <ToolComparison>/atva26-aivl.tar.gz, sets the
%   TOOLCOMPARISON_AIVL_TARBALL env var, and invokes the local
%   toolbox_install.m which extracts to userpath + wires startup.m.
%
%   Idempotent.

    here = fileparts(mfilename('fullpath'));
    v2_root = fileparts(here);
    tarball = fullfile(v2_root, 'atva26-aivl.tar.gz');

    if ~isfile(tarball)
        error('ToolComparison:aivl_install:missing_tarball', ...
              ['AIVL tarball not found at %s. ', ...
               'Either it was never staged or you do not have an AIVL license. ', ...
               'See README.md "AIVL setup" for the tarball recipe.'], tarball);
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
