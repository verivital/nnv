%% Training of an MNIST classifier (CNN)
% Code based on MathWorks example
% https://www.mathworks.com/help/deeplearning/ug/create-simple-deep-learning-network-for-classification.html

t = tic; % track total time for the training

% Load data (no download necessary)
digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
    'nndatasets','DigitDataset');
% Images
imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

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


% Select training procedure
numTrainFiles = 750; %total number of images to use
% Split the data
[imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,'randomize');

numClasses = height(labelCount); % number of classes in dataset
imgSize = [28 28 1]; % size of the images

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
    'MaxEpochs',4, ...
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
save('mnist_model.mat', 'net', 'accuracy');

toc(t);
