% Image=csvread('cifar10_test.csv');
% Image=Image(2:3073);
% load('cifar10_test.mat')
function [Label, Label_pred, Y_cnn]=test_single_image(input,net)
% tests single image given as 'input'(first column is the associated 'Label') 
% output: Y_cnn --> probabilies of each label being predicted
% Author: Neelanjana Pal 04/11/2020

Label = input(1);
Image = input(2:3073);
Image = Image/255;
Image = Image';
load('Norm.mat');
%load('cifar_cnn.mat');

if ~all(mean_data==0) && ~all(std_data==1)
    I=reshape(Image,3,1024);
    for i=1:3
        J(:,:,i)=reshape(I(i,:),32,32);
        J(:,:,i)=J(:,:,i)';
    end
    for i=1:3
      norm_im(:,:,i) = (J(:,:,i) - mean_data(i))/std_data(i);
    end
    normalized_image=norm_im;
else 
    normalized_image = Image-0.5;
end

Y_cnn=net.evaluate(normalized_image);
[~,max_idx] = max(Y_cnn);%(:,image)
Label_pred = max_idx-1;

% uncomment when running this file only
%     if Label == Label_pred
%         fprintf('Image correctly classified');
%     end
end
