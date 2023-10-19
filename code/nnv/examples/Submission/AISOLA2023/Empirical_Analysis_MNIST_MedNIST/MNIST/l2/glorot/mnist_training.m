%% 1) Load data

dataPath = "../../../../../../../../../data/MNIST/";
filenameImagesTrain = 'train-images-idx3-ubyte.gz';
filenameLabelsTrain = 'train-labels-idx1-ubyte.gz';

XData = processImagesMNIST(dataPath + filenameImagesTrain);
YData = processLabelsMNIST(dataPath + filenameLabelsTrain);

% Prepare data
rng(0); % fix random seed to select training data
% Let's use 3500 images from each class: 2500 (train), 500 (val), 500 (test)
N = 3500;
idxs = randperm(6000, N);
ntrain = 2500;
ntest = 500;
% Get all indexes into their corresponding classes
idxs_train = [];
idxs_test = [];
idxs_val = [];
for i=1:10 % 10 classes
    class_train = idxs(1:ntrain) + 6000*(i-1);
    idxs_train = [idxs_train, class_train];
    class_test = idxs(ntrain+1:ntrain+ntest) + 6000*(i-1);
    idxs_test = [idxs_test, class_test];
    class_val = idxs(ntrain+ntest+1:end) + 6000*(i-1);
    idxs_val= [idxs_val, class_val];
end

% Allocate data to corresponding variables
Xtrain = XData(:,:,:,idxs_train); Ytrain = YData(idxs_train);
Xtest  = XData(:,:,:, idxs_test); Ytest  = YData(idxs_test);
Xval   = XData(:,:, :, idxs_val); Yval   = YData(idxs_val);


%% 2) Create model

% Strat with a small model and see how it performs

% Define layers
layers = [imageInputLayer([28 28 1])
    convolution2dLayer(3, 3,'Stride',1,'WeightsInitializer','glorot', 'BiasInitializer','zeros');
    batchNormalizationLayer
    reluLayer
    averagePooling2dLayer(2,'Stride',2,'Padding',[0 0 0 1])
    flattenLayer
    fullyConnectedLayer(10, 'WeightsInitializer','glorot', 'BiasInitializer','zeros')
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
