% Compare output ranges for different reachability methods
load images.mat;
load('Small_ConvNet.mat');
nnvNet = matlab2nnv(net);

% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9

% computing reachable set

%% Reachability analysis

% Adversarial attack
IM = IM_data(:,:,1);
IM = reshape(IM, [numel(IM) 1]);
label = IM_labels(1);
% Brightening attack
d = 250; % threshold to modify pixels
lb = IM;
ub = IM;
for i=1:numel(IM)
    if  IM(i) >= d
        lb(i) = 0;
        ub(i) = 0.05*IM(i);
    end
end
% Create input set
lb = reshape(lb, [28 28 1]);
ub = reshape(ub, [28 28 1]);
inputSet = ImageZono(lb, ub); 
inputSet1 = inputSet.toImageStar;

% Define reachability parameters and compute reachable set

% 1) approx-star
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
t = tic;
outputSet_Star = nnvNet.reach(inputSet1, reachOptions);
reachTime_star = toc(t);

% 2) approx-zono
reachOptions.reachMethod = 'approx-zono';
t = tic;
outputSet_Zono = nnvNet.reach(inputSet, reachOptions);
reachTime_zono = toc(t);

% 3) abs-dom
reachOptions.reachMethod = 'abs-dom';
t = tic;
outputSet_AbsDom = nnvNet.reach(inputSet1, reachOptions);
reachTime_AbsDom = toc(t);

% Get output ranges
[lb1, ub1] = outputSet_Star.getRanges;
[lb2, ub2] = outputSet_Zono.getRanges;
[lb3, ub3] = outputSet_AbsDom.getRanges;

% Reshape ranges
lb1 = reshape(lb1, [10 1]);
ub1 = reshape(ub1, [10 1]);
lb2 = reshape(lb2, [10 1]);
ub2 = reshape(ub2, [10 1]);
lb3 = reshape(lb3, [10 1]);
ub3 = reshape(ub3, [10 1]);

% Create error bars for plotting
im_center1 = (lb1 + ub1)/2;
err1 = (ub1 - lb1)/2;
x1 = 0:1:9;
y1 = im_center1;

im_center2 = (lb2 + ub2)/2;
err2 = (ub2 - lb2)/2;
x2 = 0:1:9;
y2 = im_center2;

im_center3 = (lb3 + ub3)/2;
err3 = (ub3 - lb3)/2;
x3 = 0:1:9;
y3 = im_center3;

%% Visualization
% plot ranges of different approaches
figure;

subplot(1,3,1);
e = errorbar(x1,y1,err1);
e.LineStyle = 'none';
e.LineWidth = 1;
e.Color = 'red';
xlabel('Output', 'FontSize', 11);
ylabel('Ranges', 'FontSize', 11);
xlim([0 9]);
title('ImageStar', 'FontSize', 11);
xticks([0 5 9]);
xticklabels({'0', '5', '9'});
set(gca, 'FontSize', 10);

subplot(1,3,2);
e = errorbar(x2,y2,err2);
e.LineStyle = 'none';
e.LineWidth = 1;
e.Color = 'red';
xlabel('Output', 'FontSize', 11);
ylabel('Ranges', 'FontSize', 11);
xlim([0 9]);
title('Zonotope', 'FontSize', 11);
xticks([0 5 9]);
xticklabels({'0', '5', '9'});
set(gca, 'FontSize', 10);

subplot(1,3,3);
e = errorbar(x3,y3,err3);
e.LineStyle = 'none';
e.LineWidth = 1;
e.Color = 'red';
xlabel('Output', 'FontSize', 11);
ylabel('Ranges', 'FontSize', 11);
xlim([0 9]);
title('Polytope', 'FontSize', 11);
xticks([0 5 9]);
xticklabels({'0', '5', '9'});
set(gca, 'FontSize', 10);

