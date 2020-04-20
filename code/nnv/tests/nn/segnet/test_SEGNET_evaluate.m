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