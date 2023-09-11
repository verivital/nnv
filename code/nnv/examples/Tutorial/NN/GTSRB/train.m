%% Training of a GTSRB classifier (CNN)
% Code based on MathWorks example
% https://www.mathworks.com/help/deeplearning/ug/create-simple-deep-learning-network-for-classification.html

t = tic; % track total time for the training

% Load data
gtsrb_path = [nnvroot(), filesep, 'data', filesep, 'GTSRB', filesep];
imds = imageDatastore(gtsrb_path, 'IncludeSubfolders',true,'LabelSource','foldernames');

inputSize = [30 29];
imds.ReadFcn = @(loc)imresize(imread(loc),inputSize);

% Display some of the data
figure;
perm = randperm(10000,20);
for i = 1:20
    subplot(4,5,i);
    imshow(imds.Files{perm(i)});
end

% Show the number of images (per class) in the dataset
labelCount = countEachLabel(imds);
disp(labelCount);

% Select training and validation data
numTrainFiles = 160; % min class has 210 instances, so use ~3/4 of these
[imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,'randomize');

imdsTrain.ReadFcn = @(loc)imresize(imread(loc),inputSize);
imdsValidation.ReadFcn = @(loc)imresize(imread(loc),inputSize);

numClasses = height(labelCount); % number of classes in dataset
imgSize = [30 29 3]; % size of the images

% Create the neural network model
layers = [
    imageInputLayer(imgSize) % image size = [28 28 1] (Height x Width x Channels)
    
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(numClasses) % 10 = number of classes
    softmaxLayer
    classificationLayer];

% Training options
options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',10, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',30, ...
    'Verbose',true);

% Train network
net = trainNetwork(imdsTrain,layers,options);

% Validate network (accuracy)
YPred = classify(net,imdsValidation);
YValidation = imdsValidation.Labels;

accuracy = sum(YPred == YValidation)/numel(YValidation);
disp ("Validation accuracy = "+string(accuracy));

% Save model
disp("Saving model...");
save('gtsrb_model.mat', 'net', 'accuracy');

toc(t);
