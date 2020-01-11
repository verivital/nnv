digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
    'nndatasets','DigitDataset');
imds = imageDatastore(digitDatasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');
figure;
perm = randperm(10000,20);
for i = 1:20
    subplot(4,5,i);
    imshow(imds.Files{perm(i)});
end
labelCount = countEachLabel(imds);

img = readimage(imds,1);
size(img);

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
    
    maxPooling2dLayer(2,'Stride',2);
    
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

net = trainNetwork(imdsTrain,layers,options);


YPred = classify(net,imdsValidation);
YValidation = imdsValidation.Labels;

accuracy = sum(YPred == YValidation)/numel(YValidation)

save Small_ConvNet.mat net; 

N = 2000; % get 100 images and its labels from the imdsValidation
IM_data = zeros(28, 28, N);
IM_labels = zeros(N, 1);
for i=1:N
    if YPred(i) == YValidation(i)
        IM_data(:,:, i) = imread(imdsValidation.Files{i});
        IM_labels(i) = imdsValidation.Labels(i);
    end
end

save images.mat IM_data IM_labels
