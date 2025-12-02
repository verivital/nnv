% Testing different ways to load neural networks into NNV

%% Test 1: Load a DAGNetwork 
load('models/triangle_net.mat');
net = matlab2nnv(net);

%% Test 2: Load a SeriesNetwork 
load('models/TEST_NET.mat');
net = matlab2nnv(net);


%% Test 3: Load onnx model
if ~is_github_actions() % importers (support packages) are not installed in MATLAB actions
    % Check if ONNX support package is available
    if exist('importNetworkFromONNX', 'file')
        net_onnx = importNetworkFromONNX('mobilenetv2-1.0.onnx');
        % net = matlab2nnv(net_onnx);
        % This should get error:
        % nnet.cnn.layer.GroupedConvolution2DLayer unsupported
    else
        warning('test_io:MissingONNX', 'Skipping ONNX test: Support package not installed');
    end

    %% Test 4: Load tensorflow/keras model
    % Check if Keras/TensorFlow support package is available
    if exist('importKerasNetwork', 'file')
        try
            net_h5 = importKerasNetwork('final_model.h5');
            net = matlab2nnv(net_h5);
        catch ME
            if strcmp(ME.identifier, 'nnet_cnn:supportpackages:InstallRequired')
                warning('test_io:MissingKeras', 'Skipping Keras test: Support package not installed');
            else
                rethrow(ME);
            end
        end
    else
        warning('test_io:MissingKeras', 'Skipping Keras test: Support package not installed');
    end
end