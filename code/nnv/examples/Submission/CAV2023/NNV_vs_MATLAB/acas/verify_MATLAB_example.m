%% Verify ACAS Xu Networks in MATLAB
vnnFolder = "/home/dieman95/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";

%% Load Model
% Load the network
acasFile = "acasxu/onnx/ACASXU_run2a_2_4_batch_2000.onnx";
acas24 = importONNXNetwork(vnnFolder+acasFile, InputDataFormats='BCSS');
% These layers are not supported, need to redo the network
%      nnet.onnx.layer.ElementwiseAffineLayer    (equivalent to bias, add to previous layer)
%      nnet.onnx.layer.FlattenInto2dLayer            (simply remove it)
%      nnet.onnx.layer.VerifyBatchSizeLayer        (can simply remove it)

% Format the layers to fit MATLAB's algorithm
Layers = acas24.Layers;
Layers = Layers(5:end-1); % remove input and output layers
input_layer = featureInputLayer(5); % number of inputs to acas xu
elem_idxs = 2:3:length(Layers); % elementwiseLayers position
for k=elem_idxs
    Layers(k-1).Bias = Layers(k).Offset; % add offset as bias on previous layer
end
Layers(elem_idxs) = []; % remove elementwiselayer
acas24 = dlnetwork([input_layer; Layers]); % create dlnetwork for verification

%% Verification
% Run a verification example
X = [0;0;0;0;0];
perturbation = 0.01;
XLower = dlarray(X - perturbation, "CB");
XUpper = dlarray(X + perturbation, "CB");

% result = verifyNetworkRobustness(net, XLower, XUpper, label);
t = tic;
[YLower, YUpper] = estimateNetworkOutputBounds(acas24, XLower, XUpper);
mT = toc(t);
