clc; clear;
%% Load and parse networks into NNV
Nets = [];
load('m2nist_62iou_dilatedcnn_avgpool.mat');
net1 = SEGNET.parse(net, 'm2nist_62iou_dilatedcnn_avgpool');
Nets = [Nets net1];
load('m2nist_75iou_transposedcnn_avgpool.mat');
net2 = SEGNET.parse(net, 'm2nist_75iou_transposedcnn_avgpool');
Nets = [Nets net2];
load('m2nist_dilated_72iou_24layer.mat');
net3 = SEGNET.parse(net, 'm2nist_dilated_72iou_24layer.mat');
Nets = [Nets net3];
load('m2nist_6484_test_images.mat');

Nmax = 50; % maximum allowable number of attacked pixels
de = [0.00002; 0.00006];
Nt = 150;

%% create input set
N1 = length(de);  

IS(N1) = ImageStar;
GrTruth = cell(1,N1);
for l=1:N1
    ct = 0;
    flag = 0;
    im = im_data(:,:,l);
    at_im = im;
    for i=1:64
        for j=1:84
            if im(i,j) > Nt
                at_im(i,j) = 0;
                ct = ct + 1;
                if ct == Nmax
                    flag = 1;
                    break;
                end
            end
        end
        if flag == 1
            break;
        end
    end

    dif_im = im - at_im;
    noise = -dif_im;
    % Perform robustness analysis
    V(:,:,:,1) = double(im);
    V(:,:,:,2) = double(noise);
    C = [1; -1];
    d = [1; de(l)-1];
    S = ImageStar(V, C, d, 1-de(l), 1);
    IS(l) = S; 
    GrTruth{l} = {im};
end

%% compute reachability time of ReLU layers
N2 = length(Nets);
ReLU_ReachTime = zeros(N2, N1);
Other_ReachTime = zeros(N2, N1);
Total_ReachTime = zeros(N2, N1);

numCores = 1;

for i=1:N2
    for j=1:N1
        Nets(i).reach(IS(j), "relax-star-random", 1, 0);
        [ReLU_ReachTime(i, j), Other_ReachTime(i, j), Total_ReachTime(i, j)] = getReachTime(Nets(i));
    end
end


%% plot reachability time
fig = figure;
subplot(1,3,1);
y = [ReLU_ReachTime(1,1) Other_ReachTime(1,1); ReLU_ReachTime(1,2) Other_ReachTime(1,2)];
bar(y);
ylabel('Reachability Time');
xlabel('$\Delta_\epsilon$', 'Interpreter', 'latex');
set(gca, 'XTickLabel', {'0.00002'; '0.00006'});
title('N_4');


subplot(1,3,2);
y = [ReLU_ReachTime(2,1) Other_ReachTime(2,1); ReLU_ReachTime(2,2) Other_ReachTime(2,2)];
bar(y);
ylabel('Reachability Time');
xlabel('$\Delta_\epsilon$', 'Interpreter', 'latex');
set(gca, 'XTickLabel', {'0.00002'; '0.00006'});
title('N_5');


subplot(1,3,3);
y = [ReLU_ReachTime(3,1) Other_ReachTime(3,1); ReLU_ReachTime(3,2) Other_ReachTime(3,2)];
bar(y);
ylabel('Reachability Time');
xlabel('$\Delta_\epsilon$', 'Interpreter', 'latex');
set(gca, 'XTickLabel', {'0.00002'; '0.00006'});
title('N_6');


%% function for computing total ReLU reachTime and others reachTime
function [relu_reachTime, others_reachTime, total_reachTime] = getReachTime(net)
    total_reachTime = sum(net.reachTime);
    n = length(net.Layers);
    relu_reachTime = 0;
    others_reachTime = 0; 
    for i=1:n
        if isa(net.Layers{i}, 'ReluLayer')
            relu_reachTime = relu_reachTime + net.reachTime(i);
        else
            others_reachTime = others_reachTime + net.reachTime(i);
        end
    end
end