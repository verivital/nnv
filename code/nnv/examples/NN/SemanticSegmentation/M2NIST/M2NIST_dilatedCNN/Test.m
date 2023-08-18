M2NIST_PreProcess

load('M2NIST_dilatedCNN_avgpool.mat');
%load('M2NIST_dilatedCNN_maxpool.mat');

%Prediction on Test dataset
pxdsPred = semanticseg(pximdsTest,net,'MiniBatchSize', 64, 'WriteLocation',tempdir);
metrics = evaluateSemanticSegmentation(pxdsPred,pxdsTest);

%Prediction on Total dataset
pxdsPredTotal = semanticseg(imds,net,'MiniBatchSize', 64, 'WriteLocation',tempdir);
metrics = evaluateSemanticSegmentation(pxdsPredTotal,pxds);
