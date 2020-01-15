%% This script will be used to test the blackboox systems with the LEC1 in closed-loop
% The environment data will be fed by ground truth data recorded using the
% UUV Simulator
clc;clear;close all;

addpath(genpath(pwd))

% Load blackbox model
load('UUV_model.mat');
% % Load data to test the closed-loop system
load('dataSysVal.mat')

%% Compare the outputs of whole system (Cont+Norm+sys)
[Y, FIT, X0] = compare(sys,ida);
% Compute errors
dist = sqrt((Y.y(:,1)-ida.y(:,1)).^2 + (Y.y(:,2)-ida.y(:,2)).^2);
Mdis = max(dist); % max distance between simulated and recorded trajectories
ave = mean(dist); % average distance between simulated and recorded trajectories
er = mse(Y.y(:,1:2),ida.y(:,1:2)); % compute mean square error (mse)

% Create figure and plot results
f = figure('units','normalized','outerposition',[0 0 0.7 0.7]);
plot(Y.y(:,1),Y.y(:,2),'LineWidth',3);
hold on;
plot(ida.y(:,1),ida.y(:,2),'LineWidth',3);

ylabel('Y Position (m)');
xlabel('X Position (m)');
legend('SysID','UUV');
annotation('textbox', [0.15, 0.7, 0.1, 0.1], 'String', 'Average distance = ' + string(ave) + ' m', 'FitBoxToText','on', 'FontSize', 20);
annotation('textbox', [0.15, 0.8, 0.1, 0.1], 'String', 'Max distance = '+ string(Mdis) + ' m', 'FitBoxToText','on', 'FontSize', 20);
annotation('textbox', [0.15, 0.6, 0.1, 0.1], 'String','MSE = ' + string(er), 'FitBoxToText','on', 'FontSize', 20);
title('UUV Model Validation');
grid;
set(gca,'FontSize',20);
saveas(f,'figures/sys_val','png');
