%% Create sppFile.zip for CodeOcean
% Run this script in MATLAB after installing the ONNX converter support package
%
% Steps:
% 1. In MATLAB, go to Add-Ons > Get Add-Ons
% 2. Search for "Deep Learning Toolbox Converter for ONNX Model Format"
% 3. Install it
% 4. Run this script to create sppFile.zip
% 5. Upload sppFile.zip to CodeOcean at /code/support_packages/ or /data/support_packages/

%% Find support package location
disp('Looking for ONNX support package...');

% Check common locations
possiblePaths = {
    fullfile(matlabroot, 'toolbox', 'nnet', 'supportpackages', 'onnx')
    fullfile(matlabroot, 'SupportPackages')
    fullfile(userpath, '..', 'Add-Ons', 'Toolboxes')
    fullfile(getenv('HOME'), 'Documents', 'MATLAB', 'SupportPackages')
    fullfile(getenv('HOME'), 'Library', 'Application Support', 'MathWorks', 'MATLAB Add-Ons')
};

% Try to find using MATLAB's built-in function
try
    sppRoot = matlabshared.supportpkg.getSupportPackageRoot;
    possiblePaths = [{sppRoot}; possiblePaths];
catch
    % Ignore if this fails
end

foundPath = '';
for i = 1:length(possiblePaths)
    p = possiblePaths{i};
    if exist(p, 'dir')
        disp(['Checking: ' p]);
        % Look for ONNX-related files
        onnxFiles = dir(fullfile(p, '**', '*onnx*'));
        if ~isempty(onnxFiles)
            foundPath = p;
            disp(['Found ONNX files in: ' p]);
            break;
        end
    end
end

if isempty(foundPath)
    error(['Could not find ONNX support package. Please install it first:\n' ...
           '1. Go to MATLAB Add-Ons > Get Add-Ons\n' ...
           '2. Search for "Deep Learning Toolbox Converter for ONNX Model Format"\n' ...
           '3. Install it\n' ...
           '4. Run this script again']);
end

%% Create zip file
outputFile = fullfile(pwd, 'sppFile.zip');

disp(['Creating support package archive from: ' foundPath]);
disp(['Output file: ' outputFile]);

% Create the zip
zip(outputFile, foundPath);

disp(' ');
disp('===========================================');
disp('sppFile.zip created successfully!');
disp('===========================================');
disp(' ');
disp('Next steps:');
disp('1. Upload sppFile.zip to CodeOcean');
disp('2. Place it in /code/support_packages/ or /data/support_packages/');
disp('3. Re-run the capsule');
disp(' ');
disp(['File location: ' outputFile]);
disp(['File size: ' num2str(dir(outputFile).bytes / 1024 / 1024, '%.1f') ' MB']);
