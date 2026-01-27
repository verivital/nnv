% Convert VideoStar ONNX model to .mat file
addpath(genpath('/Users/ben/Library/Mobile Documents/com~apple~CloudDocs/Active Work/nnv/code/nnv'));

% Create output directories
if ~exist('data/FORMALISE2025/models_mat', 'dir')
    mkdir('data/FORMALISE2025/models_mat');
end

disp('Converting VideoStar model...');
onnxFile = 'data/FORMALISE2025/models/zoomin_4f.onnx';
matFile = 'data/FORMALISE2025/models_mat/zoomin_4f.mat';

disp('  Converting zoomin_4f...');
try
    netonnx = importONNXNetwork(onnxFile, "InputDataFormats", "TBCSS", "OutputDataFormats", "BC");
    net = matlab2nnv(netonnx);
    save(matFile, 'net');
    disp(['    Saved: ' matFile]);
catch ME
    disp(['    Error with matlab2nnv: ' ME.message]);
    disp('    Saving raw MATLAB network instead...');
    try
        netonnx = importONNXNetwork(onnxFile, "InputDataFormats", "TBCSS", "OutputDataFormats", "BC");
        save(matFile, 'netonnx');
        disp(['    Saved raw network: ' matFile]);
    catch ME2
        disp(['    Error: ' ME2.message]);
    end
end

disp('Done with VideoStar.');
