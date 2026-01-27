%% convert_onnx_to_mat.m
% Run this script LOCALLY (on your machine with working ONNX support)
% to convert ONNX models to .mat files for CodeOcean.
%
% This bypasses the need for the ONNX support package on CodeOcean.
%
% Usage:
%   1. Run this script in MATLAB (R2024a or R2024b with ONNX support)
%   2. Upload the generated .mat files to CodeOcean data folder
%   3. The test runners will detect and use .mat files instead of ONNX

%% FairNNV Models (Adult Census)
disp('Converting FairNNV models...');

fairnnvOnnxDir = 'data/ICAIF24/adult_onnx';
fairnnvMatDir = 'data/ICAIF24/adult_mat';

if ~exist(fairnnvMatDir, 'dir')
    mkdir(fairnnvMatDir);
end

% Convert AC-1 and AC-3 (the models used in the test)
modelList = {'AC-1', 'AC-3'};
for i = 1:length(modelList)
    modelName = modelList{i};
    onnxFile = fullfile(fairnnvOnnxDir, [modelName '.onnx']);
    matFile = fullfile(fairnnvMatDir, [modelName '.mat']);

    if exist(onnxFile, 'file')
        disp(['  Converting ' modelName '...']);
        try
            netONNX = importONNXNetwork(onnxFile, 'OutputLayerType', 'classification', 'InputDataFormats', {'BC'});
            net = matlab2nnv(netONNX);
            net.OutputSize = 2;
            save(matFile, 'net');
            disp(['    Saved: ' matFile]);
        catch ME
            disp(['    Error: ' ME.message]);
        end
    else
        disp(['    Not found: ' onnxFile]);
    end
end

%% VideoStar Model (ZoomIn)
disp('Converting VideoStar model...');

videostarOnnxDir = 'data/FORMALISE2025/models';
videostarMatDir = 'data/FORMALISE2025/models_mat';

if ~exist(videostarMatDir, 'dir')
    mkdir(videostarMatDir);
end

onnxFile = fullfile(videostarOnnxDir, 'zoomin_4f.onnx');
matFile = fullfile(videostarMatDir, 'zoomin_4f.mat');

if exist(onnxFile, 'file')
    disp('  Converting zoomin_4f...');
    try
        netonnx = importONNXNetwork(onnxFile, "InputDataFormats", "TBCSS", "OutputDataFormats", "BC");
        net = matlab2nnv(netonnx);
        save(matFile, 'net');
        disp(['    Saved: ' matFile]);
    catch ME
        disp(['    Error: ' ME.message]);
    end
else
    disp(['    Not found: ' onnxFile]);
end

%% ProbVer Model (TinyYOLO)
disp('Converting ProbVer model...');

probverOnnxDir = 'data/ProbVer/yolo_2023/onnx';
probverMatDir = 'data/ProbVer/yolo_2023/mat';

if ~exist(probverMatDir, 'dir')
    mkdir(probverMatDir);
end

onnxFile = fullfile(probverOnnxDir, 'TinyYOLO.onnx');
matFile = fullfile(probverMatDir, 'TinyYOLO.mat');

if exist(onnxFile, 'file')
    disp('  Converting TinyYOLO...');
    try
        % Note: TinyYOLO may have custom layers - try without specific format first
        netonnx = importONNXNetwork(onnxFile);
        net = matlab2nnv(netonnx);
        save(matFile, 'net');
        disp(['    Saved: ' matFile]);
    catch ME
        disp(['    Error converting to NNV: ' ME.message]);
        disp('    Saving raw MATLAB network instead...');
        try
            netonnx = importONNXNetwork(onnxFile);
            save(matFile, 'netonnx');
            disp(['    Saved raw network: ' matFile]);
        catch ME2
            disp(['    Error: ' ME2.message]);
        end
    end
else
    disp(['    Not found: ' onnxFile]);
end

disp(' ');
disp('Conversion complete!');
disp('Upload the _mat folders to CodeOcean and update the test runners to use .mat files.');
