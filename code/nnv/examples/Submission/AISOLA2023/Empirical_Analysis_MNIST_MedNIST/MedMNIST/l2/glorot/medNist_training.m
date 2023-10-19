%
% 2D classification example: MedNIST
% 6 classes: [Abdominal CT, Breast MRI, Chest X-ray, chest CT, Hand, Head CT]
%   labels : [      1          2             3           4       5      6   ]           
%              

%% 1) Load data 

% Download if necessary
dataFolder = "../../../../../../../../../data/MedNIST";

% Go through every folder (label) and load all images
categs = dir(dataFolder);
% Initialize vars
XData = zeros(64,64,1,58954); % height = width = 64, 10k imgs per class (greyscale)
YData = zeros(58954, 1);        % except BreatMRI, that has only 8954
% Load images
count = 1;
for i = 3:length(categs)-1
    label = dir(dataFolder + "/"+ string(categs(i).name));
    for k = 3:length(label)
        XData(:, :, :, count) = imread([label(k).folder '/' label(k).name]);
        YData(count) = i-2;
        count = count + 1;
    end
end
YData = categorical(YData);

% Prepare data
rng(0); % fix random seed to select training data
% Let's use 3500 images from each class: 2500 (train), 500 (val), 500 (test)
N = 3500;
idxs = randperm(10000, N);
idxs2 = randperm(8954, N);
ntrain = 2500;
ntest = 500;
% Get all indexes into their corresponding classes
idxs_train = [idxs(1:ntrain), idxs2(1:ntrain)+10000, idxs(1:ntrain)+18954,...
    idxs(1:ntrain)+28954, idxs(1:ntrain)+38954, idxs(1:ntrain)+48954];
idxs_test = [idxs(ntrain+1:ntrain+ntest), idxs2(ntrain+1:ntrain+ntest)+10000, idxs(ntrain+1:ntrain+ntest)+18954,...
    idxs(ntrain+1:ntrain+ntest)+28954, idxs(ntrain+1:ntrain+ntest)+38954, idxs(ntrain+1:ntrain+ntest)+48954];
idxs_val = [idxs(ntrain+ntest+1:end), idxs2(ntrain+ntest+1:end)+10000, idxs(ntrain+ntest+1:end)+18954,...
    idxs(ntrain+ntest+1:end)+28954, idxs(ntrain+ntest+1:end)+38954, idxs(ntrain+ntest+1:end)+48954];
% Allocate data to corresponding variables
Xtrain = XData(:,:,:,idxs_train); Ytrain = YData(idxs_train);
Xtest  = XData(:,:,:, idxs_test); Ytest  = YData(idxs_test);
Xval   = XData(:,:, :, idxs_val); Yval   = YData(idxs_val);


%% 2) Create model

% Strat with a small model and see how it performs

% Define layers
layers = [imageInputLayer([64 64 1])
    convolution2dLayer(3, 3,'Stride',1,'WeightsInitializer','glorot', 'BiasInitializer','zeros');
    batchNormalizationLayer
    reluLayer
    averagePooling2dLayer(2,'Stride',2,'Padding',[0 0 0 1])
    flattenLayer
    fullyConnectedLayer(6, 'WeightsInitializer','glorot', 'BiasInitializer','zeros')
    softmaxLayer
    classificationLayer];


%% 3) Train models

seeds = [0,1,2,3,4];

% Specify training options
options = trainingOptions('sgdm', ...
    'MaxEpochs',5, ...
    'ExecutionEnvironment','auto', ...
    'Shuffle','every-epoch', ...
    'ValidationData',{Xval,Yval},...
    'OutputNetwork', 'best-validation-loss', ...
    'L2Regularization',0.0001); % default

disp("Begin training models");
for sd = seeds

    rng(sd); % parameters should initialize differently
    net = trainNetwork(Xtrain,Ytrain,layers,options);
    
    %% 4) Test model
    
    y = classify(net,Xtest);
    accuracy = sum(y == Ytest)/numel(Ytest);
    disp("Test accuracy = " + string(accuracy));
    
    % Save model
    save("models/model_l2_glorot_"+string(sd)+".mat", 'net', 'accuracy');

end
