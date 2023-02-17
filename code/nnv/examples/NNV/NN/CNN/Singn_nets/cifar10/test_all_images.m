function test_all_images(image_filename,net)
% tests all the images from a .csv file given by 'image_filename' (1st
% column is the Labels associated to each image and column 2 to end is the
% image
% prints total number correctly classiefied images out of all the input
% images
% Author: Neelanjana Pal 04/11/2020

load('Norm.mat');
%load('cifar_cnn.mat');
image_all = csvread(image_filename);
r = size(image_all,1);
correct=0;
%Y = zeros(10, r);
for image=1:r
    Image=image_all(image,:);
    [Label, Label_pred, Y_cnn] = test_single_image(Image,net);
    if Label == Label_pred
        correct = correct+1;
%         fprintf('image : %d , Label: %d \n',image,Label);
    end
end
fprintf('total correctly classified image (out of %d) : %d\n',r,correct);
end