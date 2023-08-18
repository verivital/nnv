M2NIST_PreProcess

% Define Segmentation Network
numClasses = 11;
numFilters = 64;
filterSize = 3;
imageSize = [64,84,1];
layers = [
    imageInputLayer(imageSize)
    convolution2dLayer(filterSize,32,'DilationFactor',1,'Padding','same')
    batchNormalizationLayer
    reluLayer
    averagePooling2dLayer(2, 'Stride', 1,'Padding','same');%maxPooling2dLayer
    
    convolution2dLayer(filterSize,64,'DilationFactor',2,'Padding','same')
    batchNormalizationLayer
    reluLayer
    averagePooling2dLayer(2, 'Stride', 1,'Padding','same');
    
    convolution2dLayer(1,numClasses);
    softmaxLayer()
    pixelClassificationLayer('Name','labels','Classes',tbl.Name,'ClassWeights',classWeights)
    ];


% define optimizer
opts = trainingOptions('sgdm', ...
    'InitialLearnRate',1e-3, ...
    'MaxEpochs',30, ...
    'ExecutionEnvironment','parallel',...
    'MiniBatchSize',64,...
	'ValidationData',pximdsVal,...
    'Plots','training-progress',...
    'ValidationPatience', 4);

% train the network 
net = trainNetwork(pximdsTrain,layers,opts);

