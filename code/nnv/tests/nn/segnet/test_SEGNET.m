%to run this as a test, use results_SEGNET=runtests('test_SEGNET')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%% test 1: SEGNET evaluate

load SegNet.mat;
nnvSegNet = matlab2nnv(net);

% Load data
dataSetDir = fullfile(toolboxdir('vision'),'visiondata','triangleImages');
imageDir = fullfile(dataSetDir,'trainingImages');
labelDir = fullfile(dataSetDir,'trainingLabels');
imds = imageDatastore(imageDir);
% classNames = ["triangle","background"];
% labelIDs   = [255 0];
% pxds = pixelLabelDatastore(labelDir,classNames,labelIDs);

im = readimage(imds, 1); 
tic;
y = nnvSegNet.evaluate(im);
toc;
tic;
y1 = activations(net, im, net.Layers(31).Name);
toc;

tic;
y2 = nnvSegNet.evaluate(im, 30);
toc;
