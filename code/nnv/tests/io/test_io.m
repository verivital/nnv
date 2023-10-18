% Testing different ways to load neural networks into NNV

%% Test 1: Load a DAGNetwork 
load('models/triangle_net.mat');
net = matlab2nnv(net);

%% Test 2: Load a SeriesNetwork 
load('models/TEST_NET.mat');
net = matlab2nnv(net);


%% Test 3: Load onnx model
if ~is_github_actions() % importers (support packages) are not installed in MATLAB actions
    net_onnx = importONNXNetwork('mobilenetv2-1.0.onnx');
    % net = matlab2nnv(net_onnx);
    % This should get error: 
    % nnet.cnn.layer.GroupedConvolution2DLayer unsupported 
    
    %% Test 4: Load tensorflow/keras model
    net_h5 = importKerasNetwork('final_model.h5');
    net = matlab2nnv(net_h5);
end