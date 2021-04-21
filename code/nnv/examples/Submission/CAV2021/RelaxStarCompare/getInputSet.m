function [IS, Lbls] = getInputSet(eps)
% construct input set with specific eps
if eps==0.1
    load IS_01.mat IS_01 Labels;
    IS = IS_01;
    Lbls = Labels;
elseif eps==0.05
    load IS_005.mat IS_005 Labels;
    IS = IS_005;
    Lbls = Labels;
elseif eps==0.15
    load IS_015.mat IS_015 Labels;
    IS = IS_015;
    Lbls = Labels;
elseif eps==0.2
    load IS_02.mat IS_02 Labels;
    IS = IS_02;
    Lbls = Labels;
elseif eps==0.25
    load IS_025.mat IS_025 Labels;
    IS = IS_025;
    Lbls = Labels;
elseif eps==0.3
    load IS_03.mat IS_03 Labels;
    IS = IS_03;
    Lbls = Labels;
else
    T = csvread("mnist_test.csv");
    im_data = T(:,2:785);
    im_labels = T(:,1);
    N = size(T, 1);
    IS(N) = ImageStar;
    % number of imagestar input sets (100 imagestars)
    Lbls = im_labels + 1; % labels corresponding to the disturbanced images
    mean = 0.1307;
    std = 0.3018;
    for i=1:N
        fprintf("\nConstructing %d^th ImageStar input set", i);
        im = im_data(i,:);
        im = reshape(im, [28 28]);
        im = im';
        im = im/255;
        
        ub = im + eps;
        lb = im - eps;
        ub(ub>1) = 1;
        lb(lb<0) = 0;
        lb = (lb - mean)/std; 
        ub = (ub - mean)/std;    
        IS(i) = ImageStar(lb, ub);
    end

end

