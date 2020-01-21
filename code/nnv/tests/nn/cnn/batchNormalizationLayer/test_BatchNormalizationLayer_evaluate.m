load test_Nets.mat;

img = readimage(imdsValidation, 1);
imshow(img);

y1 = activations(MatlabNet, img, 13);
nnvNet.evaluate(img);
y2 = nnvNet.features{13};


