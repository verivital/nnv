% load images
imds = imageDatastore('images');

% class names
classNames =["zero","one","two","three","four","five","six","seven","eight","nine","ten"];
pixelLabelID = [0,1,2,3,4,5,6,7,8,9,10];
pxds = pixelLabelDatastore('masks',classNames,pixelLabelID);

%load network
load('M2NIST_dilatedCNN_avgpool.mat');

%ground_tructh image
I=readimage(imds,1129);
imshow(I)

%pixelLabeled image
C = readimage(pxds,1129);

%pixelLabeled image, overlaying it on top of ground_truth image.
B = labeloverlay(I,C);
figure
imshow(B)

%predicted pixelLabeled image
Pred = semanticseg(I,net);

%predicted pixelLabeled image, overlaying it on top of ground_truth image.
B_pred = labeloverlay(I,Pred);
figure
imshow(B_pred)