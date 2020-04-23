% load images
imds = imageDatastore('images');

% class names
classNames =["zero","one","two","three","four","five","six","seven","eight","nine","ten"];
pixelLabelID = [0,1,2,3,4,5,6,7,8,9,10];
pxds = pixelLabelDatastore('masks',classNames,pixelLabelID);

% I = readimage(imds,1);
% C = readimage(pxds,1);
% B = labeloverlay(I,C);

% count the pixels 
tbl = countEachLabel(pxds);

% fix class weighting imbalance
numberPixels = sum(tbl.PixelCount);
frequency = tbl.PixelCount / numberPixels;
classWeights = 1./ frequency;

% % Visualize by pixel counts 
% bar(1:numel(classNames),frequency);
% xticks(1:numel(classNames));
% xticklabels(tbl.Name)
% xtickangle(45);
% ylabel('Frequency');

[imdsTrain, imdsVal, imdsTest, pxdsTrain, pxdsVal, pxdsTest] = partitionData(imds,pxds);
pximdsTrain = pixelLabelImageDatastore(imdsTrain,pxdsTrain);
pximdsVal = pixelLabelImageDatastore(imdsVal,pxdsVal);
pximdsTest = pixelLabelImageDatastore(imdsTest,pxdsTest);