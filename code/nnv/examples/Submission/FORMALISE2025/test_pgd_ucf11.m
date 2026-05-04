%% Load the model and the data

modelName = sprintf("ucf11_c3d_%df_reducedparams_8outchannels.onnx", 16);
netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "BCSSS", "OutputDataFormats", "BC");
net = matlab2nnv(netonnx);
net.OutputSize = 11;
disp("Finished loading model: " + modelName);

index = 5;
data = readNPY(sprintf("data/UCF11/ucf11_grayscale_%df_verification_data.npy", 16)); % 100 16 112 112 3
labels = readNPY(sprintf("data/UCF11/ucf11_grayscale_%df_verification_labels.npy", 16));
s = data(index,:,:,:,:);
s = squeeze(s);
label = labels(index) + 1;

% get the true output (we know this is coreect)
output = net.evaluate(s);
[~, P] = max(output);

fprintf("Predicted output: %d\n", P);
fprintf("True output: %d\n", label);

%% Convert the model to dlnetwork and the data to dlarray

% prepare the model
lgraph = layerGraph(netonnx);
lastL = lgraph.Layers(end);
lgraph = removeLayers(lgraph, lastL.Name);
dlnet = dlnetwork(lgraph);

% prepare the data
s_pgd = reshape(s, 1, 1, 16, 112, 112);
s_pgd = dlarray(s_pgd, 'BCSSS');

% create the one-hot encoded vector
T = dlarray(onehotencode(label, 1, 'ClassNames', 1:11));

%% Generate the adversarial sample

% first create the one-hot encoded label vector so that the PGD knows what label to maximize loss

Xadv = PGD(dlnet, s_pgd, T, ...
    'Epsilon', 3/255, 'Alpha', 3/1020, 'Steps', 100, ...
    'DataMin', 0, 'DataMax', 1);

fprintf("Size of s: [ %d %d %d ]\n", size(s, 1), size(s, 2), size(s, 3));
fprintf("Size of Xadv: [ %d %d %d %d %d ]\n", size(Xadv, 1), size(Xadv, 2), size(Xadv, 3), size(Xadv, 4), size(Xadv, 5));

fprintf("First 5 of s: [ %f %f %f %f %f ]\n", s(8,55,55), s(4,12,76), s(2,100,100), s(12,95,15), s(16,56,15));
fprintf("First 5 of Xadv: [ %f %f %f %f %f ]\n", Xadv(8,55,55), Xadv(4,12,76), Xadv(2,100,100), Xadv(12,95,15), Xadv(16,56,15));



advOutput = net.evaluate(Xadv);
[~, advP] = max(advOutput);

fprintf("Adversarial output: %d\n", advP);

%% Helper functions
function Xadv = PGD(net, X, T, varargin)
    p = inputParser;
    addParameter(p, 'Epsilon', 1/255);
    addParameter(p, 'Alpha', 2/255);
    addParameter(p, 'Steps', 10);
    addParameter(p, 'DataMin', 0);
    addParameter(p, 'DataMax', 1);

    parse(p, varargin{:});
    eps = p.Results.Epsilon;
    alpha = p.Results.Alpha;
    K = p.Results.Steps;
    xmin = p.Results.DataMin;
    xmax = p.Results.DataMax;

    X0 = X;

    % random start
    R = (2*rand(size(X0), 'like', extractdata(X0)) - 1) * eps;
    Xadv = X0 + dlarray(R, dimsLike(X0));
    Xadv = clipToRange(Xadv, xmin, xmax);
    Xadv = projectToLinfBall(Xadv, X0, eps);

    for k=1:K
        % compute gradient of loss wrt inputs
        [loss, gX] = dlfeval(@gradWrtInput, net, Xadv, T);
        
        step = alpha * sign(gX);
        Xadv = Xadv + step;

        Xadv = projectToLinfBall(Xadv, X0, eps);
        Xadv = clipToRange(Xadv, xmin, xmax);
    end
end

function [loss, gX] = gradWrtInput(dlnet, X, T)
    Y = forward(dlnet, X);
    loss = crossentropy(Y, T, 'TargetCategories', 'independent');
    gX = dlgradient(loss, X);
end

function Xc = clipToRange(X, a, b)
    Xc = min(max(X, a), b);
end

function Xp = projectToLinfBall(X, X0, eps)
    Xp = min(max(X, X0 - eps), X0 + eps);
end

function d = dimsLike(X)
    if hasdims(X)
        d = dims(X);
    else
        d = 'SSSCB';
    end
end


% This section includes some notes about loading the model into nnv versus just as an onnx model.
% This is relevant because we need access to the gradients for PGD and I don't believe we have access through nnv.
% Although, we might have access through ONNX.
%
%
% modelName = sprintf("ucf11_c3d_%df_reducedparams_8outchannels.onnx", 16);
% netonnx = importONNXNetwork("models/" + modelName, "InputDataFormats", "BCSSS", "OutputDataFormats", "BC");
% net = matlab2nnv(netonnx);
% net.OutputSize = 11;
% 
% data = readNPY(sprintf("data/UCF11/uc11_grayscale_%df_verification_data.npy", 16));
% s = data(1,:,:,:);
%
% By this point, we've loaded both the NNV and ONNX models and the data, as well as grabbed our sample.
% 
% To do inference on the NNV model is straightforward.
% s_nnv = squeeze(s);
% output = net.evaluate(s_nnv);
%
% And we're done... now for the ONNX model. First, we need to permute the input
% data type from [1, 16, 112, 112] to [16, 112, 112, 1].
%
% s_onnx = permute(s, [2, 3, 4, 1])
% onnx_output = predict(netonnx, s_onnx); 
%
% With this, we have gotten the output logits for both the ONNX and NNV versions of the models.

