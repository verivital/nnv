%to run this as a test, use results_cnn_cnn=runtests('test_cnn_cnn')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables




%___________________________________________________________________________________________________
%tests below originally taken from test_CNN_constructor.m


%% test 1: CNN constructor

% construct FullyConnectedLayer objects
net = CNN;




%___________________________________________________________________________________________________
%tests below originally taken from test_CNN_parse.m


%% test 2: CNN parse

digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
    'nndatasets','DigitDataset');
imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

numTrainFiles = 750;
[imdsTrain,imdsValidation] = splitEachLabel(imds,numTrainFiles,'randomize');

layers = [
    imageInputLayer([28 28 1])
    
    convolution2dLayer(3,8,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer];

options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.01, ...
    'MaxEpochs',4, ...
    'Shuffle','every-epoch', ...
    'ValidationData',imdsValidation, ...
    'ValidationFrequency',30, ...
    'Verbose',false, ...
    'Plots','training-progress');

MatlabNet = trainNetwork(imdsTrain,layers,options);

nnvNet = CNN.parse(MatlabNet, 'MatlabNet');

save test_Nets.mat imds imdsTrain imdsValidation MatlabNet nnvNet;


%___________________________________________________________________________________________________
%tests below originally taken from construct_VGG16.m


%% test 3: CNN construct vgg16


% Manually construct CNN VGG16 objects
net = vgg16();













