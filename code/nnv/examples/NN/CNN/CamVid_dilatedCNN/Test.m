CamVid_PreProcess

load('Camvid_dilatedCNN.mat');

%Prediction on Test dataset
pxdsPred = semanticseg(imdsTest,net,'MiniBatchSize',32,'WriteLocation',tempdir);
metrics = evaluateSemanticSegmentation(pxdsPred,pxdsTest);

%    GlobalAccuracy    MeanAccuracy    MeanIoU    WeightedIoU    MeanBFScore
 %   ______________    ____________    _______    ___________    ___________

  %     0.67056          0.61645       0.34294      0.57816        0.34279 

%Prediction on Total dataset
pxdsPred_Total = semanticseg(imds,net,'MiniBatchSize',32,'WriteLocation',tempdir);
metrics = evaluateSemanticSegmentation(pxdsPred_Total,pxds);

%    GlobalAccuracy    MeanAccuracy    MeanIoU    WeightedIoU    MeanBFScore
 %   ______________    ____________    _______    ___________    ___________

  %     0.68468          0.66265       0.35693      0.59102        0.34934 