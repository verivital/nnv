% Convert ProbVer ONNX model to .mat file for CodeOcean
% Uses importNetworkFromONNX and removes custom ONNX layers for portability
addpath(genpath('/Users/ben/Library/Mobile Documents/com~apple~CloudDocs/Active Work/nnv/code/nnv'));

% Create output directories
if ~exist('data/ProbVer/yolo_2023/mat', 'dir')
    mkdir('data/ProbVer/yolo_2023/mat');
end

disp('Converting ProbVer model...');
onnxFile = 'data/ProbVer/yolo_2023/onnx/TinyYOLO.onnx';
matFile = 'data/ProbVer/yolo_2023/mat/TinyYOLO.mat';

disp('  Converting TinyYOLO...');
try
    % Load the ONNX network using importNetworkFromONNX (returns dlnetwork)
    net = importNetworkFromONNX(onnxFile, "InputDataFormats", "BCSS");

    % Remove custom ONNX layers that cause portability issues
    % These layers (like VerifyBatchSizeLayer) aren't needed for inference
    % and cause errors when loaded without the ONNX support package
    disp('    Removing custom ONNX layers for portability...');

    layerGraph = layerGraph(net);
    layersToRemove = {};

    for i = 1:numel(net.Layers)
        layerClass = class(net.Layers(i));
        if contains(layerClass, 'nnet.onnx.layer')
            layersToRemove{end+1} = net.Layers(i).Name;
            disp(['      Removing: ' net.Layers(i).Name ' (' layerClass ')']);
        end
    end

    % Remove the custom layers and reconnect
    for i = 1:length(layersToRemove)
        layerName = layersToRemove{i};

        % Find connections to/from this layer
        conns = layerGraph.Connections;
        srcLayers = conns.Source(strcmp(conns.Destination, layerName));
        dstLayers = conns.Destination(strcmp(conns.Source, layerName));

        % Remove the layer
        layerGraph = removeLayers(layerGraph, layerName);

        % Reconnect if there was a passthrough
        if ~isempty(srcLayers) && ~isempty(dstLayers)
            for j = 1:length(dstLayers)
                try
                    layerGraph = connectLayers(layerGraph, srcLayers{1}, dstLayers{j});
                catch
                    % Connection might already exist or be incompatible
                end
            end
        end
    end

    % Convert back to dlnetwork
    net = dlnetwork(layerGraph);

    % Save the cleaned dlnetwork
    save(matFile, 'net');
    disp(['    Saved dlnetwork: ' matFile]);
    disp(['    Network type: ' class(net)]);
    disp(['    Input size: ' mat2str(net.Layers(1).InputSize)]);
    disp(['    Number of layers: ' num2str(numel(net.Layers))]);

catch ME
    disp(['    Error: ' ME.message]);
    disp('    Stack trace:');
    for i = 1:length(ME.stack)
        disp(['      ' ME.stack(i).name ' (line ' num2str(ME.stack(i).line) ')']);
    end
end

disp('Done with ProbVer.');
