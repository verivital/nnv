function nn = onnx2nnv(onnxFile, loadOptions)
% nn = onnx2nnv(onnxFile, loadOptions)
% @nn: NNV's neural network 
% @onnxFile: neural network to import in ONNX format
% @loadOptions (optional): see importONNXLayers.m for reference on optional arguments

% Import ONNX neural networks into NNV
% Output : NN -> NNV neural network class

%% Step 1. Load the network into MATLAB

switch nargin
    % Only onnx network as input
    case 1
        try
            net = importONNXNetwork(onnxFile, 'GenerateCustomLayers', false);
        catch
            warning('Using default options. Could not load the neural network with no custom layers');
            % If error when no custom layer generation options, try setting
            % input and output type, if still error, just load the layers
            try
                net = importONNXNetwork(onnxFile, 'GenerateCustomLayers', false, 'InputDataFormat', 'BSSC', 'OutputDataFormat', 'BC');
            catch
                try
                    net = importONNXLayers(onnxFile, 'InputDataFormat', 'BSSC', 'OutputDataFormat', 'BC', 'FoldConstants', "deep");
                catch
                    try 
                        net = importONNXLayers(onnxFile, 'OutputDataFormat', 'BC', 'FoldConstants', "deep");
                    catch
                        net = importONNXLayers(onnxFile);
                    end
                end
            end
        end
    % Onnx network + loading options as inputs (Parsing inputs)
    case 2
        if ~isstruct(loadOptions)
            error('Wrong input type for input 2. loadOptions must be a struct.')
        end
        if isfield(loadOptions,'InputDataFormat')
            InputDataFormat = loadOptions.InputDataFormat;
        else
            InputDataFormat = []; % automatically detected by impotONNXLayers
        end
        if isfield(loadOptions, 'OutputDataFormat')
            OutputDataFormat = loadOptions.OutputDataFormat;
        else
            OutputDataFormat = []; % automatically detected by impotONNXLayers
        end
        if isfield(loadOptions, 'TargetNetwork')
            targetNetwork = loadOptions.TargetNetwork;
        else
            targetNetwork = 'dagnetwork'; % default
        end
        if isfield(loadOptions, 'GenerateCustomLayers')
            GenerateCustomLayers = loadOptions.GenerateCustomLayers;
        else
            GenerateCustomLayers = true;
        end
        if isfield(loadOptions, 'FoldConstants')
            foldConstants = loadOptions.FoldConstants;
        else
            foldConstants = 'deep';
        end        
        % Inputs has been parsed, now try loading the network
        try
            net = importONNXLayers(onnxFile, 'GenerateCustomLayers', GenerateCustomLayers, 'FoldConstants', foldConstants, ...
                'TargetNetwork', targetNetwork, 'InputDataFormats', InputDataFormat, 'OutputDataFormats', OutputDataFormat);
        catch
            warning('Could not load the neural network with defined loading options. Trying default options for NNV.');
            try
                net = importONNXLayers(onnxFile, 'GenerateCustomLayers', false, 'InputDataFormat', 'BSSC', 'OutputDataFormat', 'BC');
            catch
                net = importONNXLayers(onnxFile, FoldConstants="deep");
            end
        end
end

% This function may not be perfect yet, we will define everything as a NN,
% return a succesful or unseccusful note, may even work with just the list
% of layers and connections. If fullyconnected (layer 1 -> layer 2, layer
% 2 -> layer 3, ...)

% The main function to be called will be importONNXLayers(...)
% There are different arguments to be passed in the onnx importers, please
% see importONNXNetwork and importONNXLayers for more info.

%% Step 2. Convert network into NNV format
nn = matlab2nnv(net); % Can use this separate function and add it to utils, or use a NN.parse function

%% Other notes
% Something more complicated but that may be more useful in the long term
% is directly using the function
%      nnet.internal.cnn.onnx.ModelProto(Filename);
% found it inside nnet.internal.cnn.onnx.importONNXLayers.
% This returns the "raw" onnx model into matlab, then we can create the
% model ourselves. 
% Cons: will be harder to debug, to understand. It will take longer to develop
% Pros: do not rely on MATLAB to add support to some of these operations, should be more robust

end

