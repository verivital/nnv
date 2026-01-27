% Convert FairNNV ONNX models to .mat files
addpath(genpath('/Users/ben/Library/Mobile Documents/com~apple~CloudDocs/Active Work/nnv/code/nnv'));

% Create output directories
if ~exist('data/ICAIF24/adult_mat', 'dir')
    mkdir('data/ICAIF24/adult_mat');
end

% Convert FairNNV models
disp('Converting FairNNV models...');
modelList = {'AC-1', 'AC-3'};
for i = 1:length(modelList)
    modelName = modelList{i};
    onnxFile = ['data/ICAIF24/adult_onnx/' modelName '.onnx'];
    matFile = ['data/ICAIF24/adult_mat/' modelName '.mat'];

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
end

disp('Done with FairNNV.');
