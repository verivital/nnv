%to run this as a test, use results_layers_batchNorm=runtests('test_layers_batchNorm')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables




%___________________________________________________________________________________________________
%tests below originally taken from test_BatchNormalizationLayer_evaluate.m


%% test 1: CNN constructor
load test_Nets.mat;

img = readimage(imdsValidation, 1);
imshow(img);

y1 = activations(MatlabNet, img, 13);

nnvNet.evaluate(img);
y2 = nnvNet.features{13};
