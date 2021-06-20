% function to check Imagestar based robustness for a given epsilon value
% input: epsilon  -->  
%        mean_data--> mean considered in the network model (from ERAN)
%        std_data --> std considered in the network model (from ERAN)

% Dung Tran 5/4/2020


%% load data 
M = csvread('classified_images.csv');
Labels = csvread('classified_Labels.csv') + 1;
load('cifar_cnn.mat');
load Norm.mat;

if ~all(mean_data==0) && ~all(std_data==1)
    mean_data=mean_data;
    std_data=std_data;
else
    mean_data=[0.5,0.5,.5];
    std_data=[1,1,1];
end


%% construct input set;
r = size(M,1);
% delta = [0.005; 0.008; 0.01; 0.012; 0.015];
delta = [0.02]; % used for testing
N = length(delta);
inputSet = cell(1,N); % a cell of imagestar input set
labelSet = cell(1,N); % a cell of correct label corresponding to imagestar input set
avr_num_pixels_attacked = zeros(1, N); % average number of attacked pixels
fprintf("\nConstructing Input Sets...");
for i=1:N
    IS = [];
    label = [];
    ct = 0;
    fprintf("\nConstructing %d^th Input Set...", i);
    for j=1:r
        IS1 = attack(M(j,:), delta(i), mean_data, std_data);
        ct = ct + IS1.numPred;  
        label = [label Labels(j)];
        IS = [IS IS1];
    end
    inputSet{i} = IS;
    labelSet{i} = label;
    avr_num_pixels_attacked(i) = ct/r; 
end
fprintf("\nFinish Constructing Input Sets");
save ImageStarIntputSet.mat inputSet labelSet avr_num_pixels_attacked; 


%% evaluate the robustness
rb = zeros(N, 1); % robustness value = number of correctly classified image / total images
vt = zeros(N, 1); % verification time
n_correct = zeros(N, 1);
avr_num_att_pixels = round(avr_num_pixels_attacked');
numCores = 4; % number of cores used in computation
fprintf("\nEvaluating Robustness...");
for i=1:N
    t = tic;
    rb(i) = net.evaluateRobustness(inputSet{i}, labelSet{i}, 'approx-star', numCores);
    n_correct(i) = round(r*rb(i));
    vt(i) = toc(t);
end

%% print the results
fprintf("\n\n======================VERIFICATION RESULTS USING NNV==========================");
T = table(delta, avr_num_att_pixels, rb, n_correct, vt);
display(T);
writetable(T, "results.txt"); % use readtable to read the table back to matlab, e.g. T = readtable("results.txt");


%% brightenning attack
function IS = attack(in_image, delta, means, stds)

    IM = in_image/255;
    n = size(IM, 2);
    lb = zeros(1, n);
    ub = zeros(1, n);
    for i=1:n
        if IM(i) >= 1 - delta
            lb(i) = IM(i);
            ub(i) = 1;
        else
            lb(i) = IM(i);
            ub(i) = IM(i);
        end
    end
    
    lb1 = reshape(lb, 3, 1024);
    ub1 = reshape(ub, 3, 1024);
    
    for i=1:3
        lb2(:,:,i) = reshape(lb1(i,:), 32, 32);
        ub2(:,:,i) = reshape(ub1(i, :), 32, 32);
    end
    
    for i=1:3
        lb3(:,:,i) = (lb2(:,:,i)-means(i))/stds(i);
        ub3(:,:,i) = (ub2(:,:,i)-means(i))/stds(i);
    end
    
    IS = ImageStar(lb3, ub3);
    
end

function [im, im1] = normalize(in_image, means, stds)
    IM = in_image/255;
    IM_r = reshape(IM,3,1024);
    for j=1:3
        J_im(:,:,j) = reshape(IM_r(j,:),32,32);
    end
    im1 = J_im;
    for k=1:3
        im(:,:,k) = (J_im(:,:,k)-means(k))/stds(k);
    end
end
