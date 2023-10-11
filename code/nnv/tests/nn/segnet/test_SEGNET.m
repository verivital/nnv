%% test 1: SEGNET evaluate

load('../../io/models/triangle_net.mat');
nnvSegNet = matlab2nnv(net);

% Load data
dataSetDir = fullfile(toolboxdir('vision'),'visiondata','triangleImages');
imageDir = fullfile(dataSetDir,'trainingImages');
labelDir = fullfile(dataSetDir,'trainingLabels');
imds = imageDatastore(imageDir);

im = readimage(imds, 1);
im = single(im);

tic;
y = nnvSegNet.evaluate(im);
toc;

tic;
y1 = activations(net, im, net.Layers(31).Name);
toc;
