% Examples: load neural networks into NNV

%% Example 1: Load a DAGNetwork (trained in MATLAB)
load('../../../tests/io/models/triangle_net.mat'); % semantic segmentation
net1 = matlab2nnv(net);

%% Example 2: Load a SeriesNetwork (trained in MATLAB)
load('../../../tests/io/models/TEST_NET.mat');
net2 = matlab2nnv(net);

%% Example 3: Load onnx model
% Check if ONNX support package is available
if exist('importNetworkFromONNX', 'file')
    acas_path = [nnvroot(), filesep, 'data', filesep, 'ACASXu', filesep];
    networks = dir(acas_path+"onnx/*.onnx");

    net_acas = importNetworkFromONNX([networks(1).folder filesep networks(1).name], "InputDataFormats","BCSS");
    net3 = matlab2nnv(net_acas);
else
    warning('load_models:MissingONNX', ...
        'Skipping ONNX example: Deep Learning Toolbox Converter for ONNX Model Format not installed');
    net3 = [];
end

%% Example 4: Load tensorflow/keras model
% Check if Keras/TensorFlow support package is available
if exist('importKerasNetwork', 'file')
    try
        net_h5 = importKerasNetwork('../../../tests/io/models/final_model.h5');
        net4 = matlab2nnv(net_h5);
    catch ME
        if strcmp(ME.identifier, 'nnet_cnn:supportpackages:InstallRequired')
            warning('load_models:MissingKeras', ...
                'Skipping Keras example: Deep Learning Toolbox Converter for TensorFlow Models not installed');
            net4 = [];
        else
            rethrow(ME);
        end
    end
else
    warning('load_models:MissingKeras', ...
        'Skipping Keras example: Deep Learning Toolbox Converter for TensorFlow Models not installed');
    net4 = [];
end