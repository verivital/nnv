%to run this as a test, use results_SEGNET=runtests('test_SEGNET')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%shared variables




%___________________________________________________________________________________________________
%tests below originally taken from test_SEGNET_evaluate.m


%% test 1: SEGNET evaluate

load SegNet.mat;
nnvSegNet = SEGNET.parse(net, 'SetNet');
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




%___________________________________________________________________________________________________
%tests below originally taken from test_SEGNET_parse.m


%% test 2: SEGNET parse

load SegNet.mat;
nnvSegNet = SEGNET.parse(net, 'SetNet');
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
