%% Load images and construct imagestar input set

T = csvread("mnist_test.csv");
im_data = T(:,2:785);
im_labels = T(:,1);

N = size(T, 1);
% number of imagestar input sets (100 imagestars)
IS_01(N) = ImageStar;
IS_015(N) = ImageStar;
IS_02(N) = ImageStar; 
IS_025(N) = ImageStar;
IS_03(N) = ImageStar;
Labels = im_labels + 1; % labels corresponding to the disturbanced images
mean = 0.1307;
std = 0.3018;
for i=1:N
    fprintf("\nConstructing %d^th ImageStar input set", i);
    im = im_data(i,:);
    im = reshape(im, [28 28]);
    im = im';
    im = im/255;
    eps = 0.1; 
    ub = im + eps;
    lb = im - eps;
    ub(ub>1) = 1;
    lb(lb<0) = 0;
    lb = (lb - mean)/std; 
    ub = (ub - mean)/std;    
    IS_01(i) = ImageStar(lb, ub);
    
    eps = 0.15; 
    ub = im + eps;
    lb = im - eps;
    ub(ub>1) = 1;
    lb(lb<0) = 0;
    lb = (lb - mean)/std; 
    ub = (ub - mean)/std;    
    IS_015(i) = ImageStar(lb, ub);
    
    eps = 0.2; 
    ub = im + eps;
    lb = im - eps;
    ub(ub>1) = 1;
    lb(lb<0) = 0;
    lb = (lb - mean)/std; 
    ub = (ub - mean)/std;    
    IS_02(i) = ImageStar(lb, ub);
    
    eps = 0.25; 
    ub = im + eps;
    lb = im - eps;
    ub(ub>1) = 1;
    lb(lb<0) = 0;
    lb = (lb - mean)/std; 
    ub = (ub - mean)/std;    
    IS_025(i) = ImageStar(lb, ub);
    
    eps = 0.3; 
    ub = im + eps;
    lb = im - eps;
    ub(ub>1) = 1;
    lb(lb<0) = 0;
    lb = (lb - mean)/std; 
    ub = (ub - mean)/std;    
    IS_03(i) = ImageStar(lb, ub);
end

save IS_01.mat IS_01 Labels;
save IS_015.mat IS_015 Labels;
save IS_02.mat IS_02 Labels;
save IS_025.mat IS_025 Labels;
save IS_03.mat IS_03 Labels;
