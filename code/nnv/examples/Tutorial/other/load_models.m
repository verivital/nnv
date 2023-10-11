% Examples: load neural networks into NNV

%% Example 1: Load a DAGNetwork (trained in MATLAB)
load('../../../tests/io/models/triangle_net.mat'); % semantic segmentation
net1 = matlab2nnv(net);

%% Example 2: Load a SeriesNetwork (trained in MATLAB)
load('../../../tests/io/models/TEST_NET.mat');
net2 = matlab2nnv(net);

%% Example 3: Load onnx model
acas_path = [nnvroot(), filesep, 'data', filesep, 'ACASXu', filesep];
networks = dir(acas_path+"onnx/*.onnx");

net_acas = importNetworkFromONNX([networks(1).folder filesep networks(1).name], "InputDataFormats","BCSS");
net3 = matlab2nnv(net_acas);

%% Example 4: Load tensorflow/keras model
net_h5 = importKerasNetwork('../../../tests/io/models/final_model.h5');
net4 = matlab2nnv(net_h5);