clc; clear;
%% Load and parse networks into NNV
load('m2nist_62iou_dilatedcnn_avgpool.mat');
Nets = SEGNET.parse(net, 'm2nist_62iou_dilatedcnn_avgpool');
load('m2nist_75iou_transposedcnn_avgpool.mat');
N2 = SEGNET.parse(net, 'm2nist_75iou_transposedcnn_avgpool');
Nets = [Nets N2];
load('m2nist_dilated_72iou_24layer.mat');
N3 = SEGNET.parse(net, 'm2nist_dilated_72iou_24layer.mat');
Nets = [Nets N3];

load('m2nist_6484_test_images.mat');

Nmax = [5 10 15 20 25]; % maximum allowable number of attacked pixels 
% Nmax = 30;
M = length(Nmax);
L = length(Nets);

%% create input set
de = 0.00001; % size of input
N = 20; % number of tested images 
IS = cell(1, M); % imagestar input set caused by the attack
GrTruth = cell(1, M); % ground truth images
for k=1:M
    IMS(N) = ImageStar;
    GrTr = cell(1,N);
    for l=1:N
        ct = 0;
        flag = 0;
        im = im_data(:,:,l);
        at_im = im;
        for i=1:64
            for j=1:84
                if im(i,j) > 150
                    at_im(i,j) = 0;
                    ct = ct + 1;
                    if ct == Nmax(k)
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
        d = [1; de-1];
        S = ImageStar(V, C, d, 1-de, 1);
        IMS(l) = S; 
        GrTr{l} = im;
    end
    IS{k} = IMS;
    GrTruth{k} = GrTr;
end

%% Verify networks
avg_RIoU = zeros(L, M); % average IoU corresponding to the number of attacked pixels
avg_RV = zeros(L, M); % average RV corresponding to number of attacked pixels
avg_RS = zeros(L, M); % average RS corresponding to number of attacked pixels
avg_numRbPixels = zeros(L, M); % average number of robust pixels
avg_numMisPixels = zeros(L, M); % average number of misclassified pixels
avg_numAttPixels = zeros(L, M); % average number of attacked pixels
avg_numUnkPixels = zeros(L, M); % average number of unknown pixels
VT = zeros(L, M); % verification time

c = parcluster('local');
numCores = c.NumWorkers;
poolobj = gcp('nocreate');
delete(poolobj); % reset parpool

% verify L networks in the Nets array
t1 = tic;
for i=1:L
    for k=1:M
        t = tic;
        Nets(i).verify(IS{k}, GrTruth{k}, 'approx-star', numCores);
        avg_RIoU(i, k) = sum(Nets(i).RIoU)/(N);
        avg_RV(i, k) = sum(Nets(i).RV)/(N);
        avg_RS(i, k) = sum(Nets(i).RS)/(N);
        avg_numRbPixels(i,k) = sum(Nets(i).numRbPixels)/(N);
        avg_numMisPixels(i,k) = sum(Nets(i).numMisPixels)/(N);
        avg_numAttPixels(i,k) = sum(Nets(i).numAttPixels)/(N);
        avg_numUnkPixels(i,k) = sum(Nets(i).numUnkPixels)/(N);
        VT(i, k) = toc(t);
    end
end
total_VT = toc(t1);
 
%% Plot Robustness Statistics

% plot (average) robustness value, robustness sensitivity, and verification time
fig = figure; 
color = ['b', 'r', 'g', 'c', 'k', 'y'];
markers= {'--*', '--o', '--x', '--s', '--p', '--d'};
titles = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)'};
labels = {'$N_4$', '$N_5$', '$N_6$', '$N_7$', '$N_8$', '$N_{9}$'};

subplot(2,3,1); % avg_RV
for i=1:L
    p = plot(avg_numAttPixels(i,:), avg_RV(i,:), markers{i}, 'Color', color(i));
    xlabel('$\overline{N}_{attackedpixels}$', 'interpreter', 'latex');
    ylabel('$\overline{RV}$', 'interpreter', 'latex');
    xticks(Nmax);
    xlim([Nmax(1) Nmax(M)]);
    title(titles{1});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,2); % avg_RS
for i=1:L
    p = plot(avg_numAttPixels(i,:), avg_RS(i,:), markers{i}, 'Color', color(i));
    xlabel('$\overline{N}_{attackedpixels}$', 'interpreter', 'latex');
    ylabel('$\overline{RS}$', 'interpreter', 'latex');
    xticks(Nmax);
    xlim([Nmax(1) Nmax(M)]);
    title(titles{2});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,3); % avg_RIoU
for i=1:L
    p = plot(avg_numAttPixels(i,:), avg_RIoU(i,:), markers{i}, 'Color', color(i));
    xlabel('$\overline{N}_{attackedpixels}$', 'interpreter', 'latex');
    ylabel('$\overline{R}_{IoU}$', 'interpreter', 'latex');
    xticks(Nmax);
    xlim([Nmax(1) Nmax(M)]);
    title(titles{3});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,4); % avg_numRbPixels
for i=1:L
    p = plot(avg_numAttPixels(i,:), avg_numRbPixels(i,:), markers{i}, 'Color', color(i));
    xlabel('$\overline{N}_{attackedpixels}$', 'interpreter', 'latex');
    ylabel('$\overline{N}_{robustpixels}$', 'interpreter', 'latex');
    xticks(Nmax);
    xlim([Nmax(1) Nmax(M)]);
    title(titles{4});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,5); % avg_numMisPixels (unrobust pixels under the attack)
for i=1:L
    p = plot(avg_numAttPixels(i,:), avg_numMisPixels(i,:), markers{i}, 'Color', color(i));
    xlabel('$\overline{N}_{attackedpixels}$', 'interpreter', 'latex');
    ylabel('$\overline{N}_{unrobustpixels}$', 'interpreter', 'latex');
    xticks(Nmax);
    xlim([Nmax(1) Nmax(M)]);
    title(titles{5});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,6); % avg_numUnkPixels (average number of unknowned pixels)
for i=1:L
    p = plot(avg_numAttPixels(i,:), avg_numUnkPixels(i,:), markers{i}, 'Color', color(i));
    xlabel('$\overline{N}_{attackedpixels}$', 'interpreter', 'latex');
    ylabel('$\overline{N}_{unknownpixels}$', 'interpreter', 'latex');
    xticks(Nmax);
    xlim([Nmax(1) Nmax(M)]);
    title(titles{6});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

saveas(fig, 'compare_m2nist_nets_vs_num_attackedpixels.pdf');
%% plot reach time

N1 = Nets(1);
N2 = Nets(2);
N3 = Nets(3);

N1_relu_reachTime = N1.reachTime(4) + N1.reachTime(8) + N1.reachTime(12);
N1_others_reachTime = sum(N1.reachTime) - N1_relu_reachTime; 

relu_id = [3 5 8 10 13 15 17 19];
N2_relu_reachTime = 0;
for i=1:length(relu_id)
    N2_relu_reachTime = N2_relu_reachTime + N2.reachTime(relu_id(i));
end
N2_others_reachTime = sum(N2.reachTime) - N2_relu_reachTime;

relu_id = [4 8 12 16 20];
N3_relu_reachTime = 0;
for i=1:length(relu_id)
    N3_relu_reachTime = N3_relu_reachTime + N3.reachTime(relu_id(i));
end
N3_others_reachTime = sum(N3.reachTime) - N3_relu_reachTime;

fig = figure;
subplot(1,2,1);
for i=1:L
    p = plot(avg_numAttPixels(i,:), VT(i,:), markers{i}, 'Color', color(i));
    xlabel('$\overline{N}_{attackedpixels}$', 'interpreter', 'latex');
    ylabel('$VT$', 'interpreter', 'latex');
    xticks(Nmax);
    xlim([Nmax(1) Nmax(M)]);
    title(titles{1});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;
subplot(1,2,2);
y = [N1_relu_reachTime N1_others_reachTime; N2_relu_reachTime N2_others_reachTime; N3_relu_reachTime N3_others_reachTime];
bar(y)
set(gca, 'XTickLabel', ['N_4'; 'N_5'; 'N_6']);
title(titles{2});

%% plot verified output set

Nets(1).plotVerifiedOutputSet(4);
Nets(2).plotVerifiedOutputSet(4);