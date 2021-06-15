clc;clear;
%% Load and parse networks into NNV
load('net_mnist_3_relu.mat');
Nets = SEGNET.parse(net, 'net_mnist_3_relu_avgpool');
load('net_mnist_3_relu_maxpool.mat');
N2 = SEGNET.parse(net, 'net_mnist_3_relu_maxpool');
Nets = [Nets N2];
load('mnist_dilated_net_21_later_83iou');
N3 = SEGNET.parse(net, 'mnist_dilated_net_21_later_83iou');
Nets = [Nets N3];
load('test_images.mat');
poolobj = gcp('nocreate');
delete(poolobj); % reset parpool

Nmax = 20; % maximum allowable number of attacked pixels
de = [0.001 0.0015 0.002 0.0025 0.003]; % size of input set
M = length(de);
L = length(Nets);

%% create input set
N = 10; % number of tested images 
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
        for i=1:28
            for j=1:28
                if im(i,j) > 150
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
        d = [1; de(k)-1];
        S = ImageStar(V, C, d, 1-de(k), 1);
        IMS(l) = S; 
        GrTr{l} = im;
    end
    IS{k} = IMS;
    GrTruth{k} = GrTr;
end

%% Verify networks
avg_RIoU = zeros(L, M); % average RIoU corresponding to the input size
avg_RV = zeros(L, M); % average RV corresponding to the input size
avg_RS = zeros(L, M); % average RS corresponding to the input size
avg_numRbPixels = zeros(L, M); % average number of robust pixels
avg_numMisPixels = zeros(L, M); % average number of misclassified pixels
avg_numAttPixels = zeros(L, M); % average number of attacked pixels
avg_numUnkPixels = zeros(L, M); % average number of unknown pixels
VT = zeros(L, M); % verification time

c = parcluster('local');
numCores = c.NumWorkers;

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
labels = {'$N_1$', '$N_2$', '$N_3$', '$N_4$', '$N_5$', '$N_6$'};

subplot(2,3,1); % avg_RV
for i=1:L
    p = plot(de, avg_RV(i,:), markers{i}, 'Color', color(i));
    ylabel('$\overline{RV}$', 'interpreter', 'latex');
    xlabel('$\Delta_{\epsilon}$', 'interpreter', 'latex');
    xticks(de);
    xlim([de(1) de(M)]);
    title(titles{1});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,2); % avg_RS
for i=1:L
    p = plot(de, avg_RS(i,:), markers{i}, 'Color', color(i));
    ylabel('$\overline{RS}$', 'interpreter', 'latex');
    xlabel('$\Delta_{\epsilon}$', 'interpreter', 'latex');
    xticks(de);
    xlim([de(1) de(M)]);
    title(titles{2});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,3); % avg_RIoU
for i=1:L
    p = plot(de, avg_RIoU(i,:), markers{i}, 'Color', color(i));
    ylabel('$\overline{R}_{IoU}$', 'interpreter', 'latex');
    xlabel('$\Delta_{\epsilon}$', 'interpreter', 'latex');
    xticks(de);
    xlim([de(1) de(M)]);
    title(titles{3});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,4); % avg_numRbPixels
for i=1:L
    p = plot(de, avg_numRbPixels(i,:), markers{i}, 'Color', color(i));
    ylabel('$\overline{N}_{robustpixels}$', 'interpreter', 'latex');
    xlabel('$\Delta_{\epsilon}$', 'interpreter', 'latex');
    xticks(de);
    xlim([de(1) de(M)]);
    title(titles{4});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,5); % avg_numMisPixels (unrobust pixels under the attack)
for i=1:L
    p = plot(de, avg_numMisPixels(i,:), markers{i}, 'Color', color(i));
    ylabel('$\overline{N}_{unrobustpixels}$', 'interpreter', 'latex');
    xlabel('$\Delta_{\epsilon}$', 'interpreter', 'latex');
    xticks(de);
    xlim([de(1) de(M)]);
    title(titles{5});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

subplot(2,3,6); % avg_numUnkPixels (average number of unknowned pixels)
for i=1:L
    p = plot(de, avg_numUnkPixels(i,:), markers{i}, 'Color', color(i));
    ylabel('$\overline{N}_{unknownpixels}$', 'interpreter', 'latex');
    xlabel('$\Delta_{\epsilon}$', 'interpreter', 'latex');
    xticks(de);
    xlim([de(1) de(M)]);
    title(titles{6});
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex');
hold off;

saveas(fig, 'compare_mnist_nets_vs_inputsize.pdf');

fig2 = figure;
for i=1:L
    p = plot(de, VT(i,:), markers{i}, 'Color', color(i));
    ylabel('$VT$', 'interpreter', 'latex');
    xlabel('$\Delta_{\epsilon}$', 'interpreter', 'latex');
    xticks(de);
    xlim([de(1) de(M)]);
    hold on;
end
legend(labels{1:L}, 'interpreter', 'latex', 'FontSize', 13);
hold off;
ax = gca;
ax.FontSize = 13;
saveas(fig2, 'VT_mnist_nets_vs_inputsize.pdf');

%% plot verified output set

Nets(1).plotVerifiedOutputSet(2);
Nets(2).plotVerifiedOutputSet(2);