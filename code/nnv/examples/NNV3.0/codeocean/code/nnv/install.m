fprintf('\nINSTALLING NNV....');
fprintf('\nInstalling tbxmanager (requires MATLAB R2023a or later) ...');
tbxmanager_folder = 'tbxmanager';

root_folder = pwd();

list = dir;
if ~isfolder(tbxmanager_folder)
    mkdir(tbxmanager_folder);
end

% install mpt toobox and other dependencies
cd(tbxmanager_folder);
urlwrite('https://raw.githubusercontent.com/verivital/tbxmanager/master/tbxmanager.m', 'tbxmanager.m');
%urlwrite('http://www.tbxmanager.com/tbxmanager.m', 'tbxmanager.m');
tbxmanager
savepath
fprintf('\nInstalling tbxmanager toolbox is done!');

cd(root_folder);

fprintf('\nInstalling mpt toolbox and other dependencies...\n');
tbxmanager install mpt mptdoc;
tbxmanager install lcp hysdel cddmex clpmex glpkmex fourier sedumi;
% tbxmanager install yalmip; % todo: error due to license, need to force acceptance
fprintf('\nInstalling dependencies is done!');
adjust_glpk;

startup_nnv; % adding dependencies and nnv to the path

%% Post-Install Verification
fprintf('\n--- Verifying Installation ---\n');
install_ok = true;

% Test 1: Check Star class exists
if exist('Star', 'class') == 8
    fprintf('[OK] Star class available\n');
else
    fprintf('[FAIL] Star class not found\n');
    install_ok = false;
end

% Test 2: Create a simple Star set
try
    S = Star([0;0], [1;1]);
    fprintf('[OK] Star set creation works\n');
catch
    fprintf('[FAIL] Star set creation failed\n');
    install_ok = false;
end

% Test 3: Check NN class exists
if exist('NN', 'class') == 8
    fprintf('[OK] NN class available\n');
else
    fprintf('[FAIL] NN class not found\n');
    install_ok = false;
end

% Test 4: Check a layer works
try
    W = eye(2);
    b = [0; 0];
    layer = FullyConnectedLayer(W, b);
    out = layer.evaluate([1; 2]);
    fprintf('[OK] FullyConnectedLayer works\n');
catch
    fprintf('[FAIL] FullyConnectedLayer failed\n');
    install_ok = false;
end

fprintf('--- End Verification ---\n\n');

if install_ok
    fprintf('Installation SUCCESSFUL!\n');
else
    fprintf('Installation completed with WARNINGS - some features may not work.\n');
    fprintf('Run check_nnv_setup() for detailed diagnostics.\n');
end

fprintf('\nPlease go to examples or test folders to run case studies and test examples.');
fprintf('\nTHANK YOU FOR TRYING NNV!\n');
fprintf('For issues or suggestions, visit: github.com/verivital/nnv/issues\n\n');

% show toolboxes
ver()

function adjust_glpk()
    % Platform-specific path handling for glpkmex
    try
        arch = computer('arch');  % 'glnxa64', 'win64', 'maci64', or 'maca64'

        % Construct platform-specific path
        glpk_base = fullfile('tbxmanager', 'toolboxes', 'glpkmex', '1.0', arch);
        glpk_folder = sprintf('glpkmex_1_0_%s', arch);
        filename = fullfile(glpk_base, glpk_folder, 'glpk.m');

        % Check if file exists for this platform
        if ~isfile(filename)
            fprintf('Note: glpkmex not found for platform %s (may not be needed)\n', arch);
            return;
        end

        fid = fopen(filename);
        if fid == -1
            return;  % Could not open file
        end
        cac = textscan(fid, '%s', 'Delimiter', '\n', 'whitespace', '');
        fclose(fid);

        filename2 = fullfile(glpk_base, glpk_folder, 'glpk2.m');
        fid = fopen(filename2, 'w');
        if fid == -1
            return;  % Could not create temp file
        end

        change_here = 372;
        for jj = 1 : min(change_here-1, length(cac{1}))
            fprintf(fid, '%s\n', cac{1}{jj});
        end
        fprintf(fid, '%s\n', '%clear glpkcc;');
        for jj = change_here+1 : length(cac{1})
            fprintf(fid, '%s\n', cac{1}{jj});
        end
        fclose(fid);
        movefile(filename2, filename, 'f');
    catch
        % Silently ignore errors - glpk adjustment is optional
    end
end
