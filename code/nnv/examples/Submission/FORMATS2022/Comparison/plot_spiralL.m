%% Plot one scenario for the linear spiral 2D


%% Gotube (no results)


%% NNV
% sp1_nnv = load("results/spiral_0.01.mat");
% sp2_nnv = load("results/spiral_0.05.mat");
sp3_nnv = load("results/spiral_0.1.mat");

%% JuliaReach
julia_path = "juliareach/results/";
% sp1_julia = load(julia_path + "spiral4.mat");
% sp2_julia = load(julia_path + "spiral5.mat");
sp3_julia = load(julia_path + "spiral6.mat");


%% Generate Plots

figure;
hold on;
Star.plotBoxes_2D_noFill(sp3_nnv.Rall(1:end-1),1,2,'k');
% plot_juliareach(sp3_julia,1,2,'b');
plot_juliareach(sp3_julia,1,2,'b');
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
legend('off')
grid;
% pg = plot(2, 0,'r');
% pn = plot(2, 0,'k');
% pj = plot(2, 0,'b');
% pf = plot(2, 0,'Color',[0 0.4 0]);
ax = gca; % Get current axis
ax.XAxis.FontSize = 14; % Set font size of axis
ax.YAxis.FontSize = 14;
% legend([pf pn pj],{'Flow*', 'NNV (ours)', 'Juliareach'},"Location","best",'FontSize',13);
run flowstar/results/spiralL3.m
legend('off')
hold on;
pn = plot(2, 0,'k');
pj = plot(2, 0,'b');
pf = plot(2, 0,'Color',[0 0.4 0]);
legend([pf pn pj],{'Flow*', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
saveas(gcf,'spiralL_compare_0.1.pdf');
pause(0.01);
xlim([0.23 0.49]);
ylim([0.48 0.7]);
saveas(gcf,'spiralL_compare_0.1_zoom.pdf');

% f1 = figure;
% hold on;
% plot_gotube(sp3_gt(:,end),1,2,2,'r'); 
% Star.plotBoxes_2D_noFill(sp3_nnv.Rall(end-1),1,2,'k');
% plot_juliareach(sp3_julia,1,2,'b',"last");
% % Add legend and labels
% xlabel('x_1');
% ylabel('x_2');
% legend('off')
% grid;
% pg = plot(2, 0,'r');
% pn = plot(2, 0,'k');
% pj = plot(2, 0,'b');
% % pf = plot(2, 0,'Color',[0 0.4 0]);
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 14; % Set font size of axis
% ax.YAxis.FontSize = 14;
% xlim([0.2 0.9])
% ylim([0.2 0.8])
% legend([pg pn pj],{'GoTube', 'NNV (ours)', 'Juliareach'},"Location","best",'FontSize',13);
% saveas(f1,'spiralL_compare_0.1_last.pdf');

% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------

% NNV
sp2_nnv = load("results/spiral_0.05.mat");

% JuliaReach
julia_path = "juliareach/results/";
sp2_julia = load(julia_path + "spiral5.mat");

figure;
hold on;
% plot_gotube(sp2_gt,1,2,2,'r'); 
Star.plotBoxes_2D_noFill(sp2_nnv.Rall(1:end-1),1,2,'k');
plot_juliareach(sp2_julia,1,2,'b');
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
% legend('off')
grid;
% pg = plot(2, 0,'r');
% pn = plot(2, 0,'k');
% pj = plot(2, 0,'b');
% pf = plot(2, 0,'Color',[0 0.4 0]);
ax = gca; % Get current axis
ax.XAxis.FontSize = 14; % Set font size of axis
ax.YAxis.FontSize = 14;
% legend([pf pn pj],{'Flow*', 'NNV (ours)', 'Juliareach'},"Location","best",'FontSize',13);
run flowstar/results/spiralL2.m
legend('off')
hold on;
pn = plot(2, 0,'k');
pj = plot(2, 0,'b');
pf = plot(2, 0,'Color',[0 0.4 0]);
legend([pf pn pj],{'Flow*', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
saveas(gcf,'spiralL_compare_0.05.pdf');
pause(0.01);
xlim([0.23 0.49]);
ylim([0.48 0.7]);
saveas(gcf,'spiralL_compare_0.05_zoom.pdf');

% f1 = figure;
% hold on;
% plot_gotube(sp2_gt(:,end),1,2,2,'r'); 
% Star.plotBoxes_2D_noFill(sp2_nnv.Rall(end-1),1,2,'k');
% plot_juliareach(sp2_julia,1,2,'b',"last");
% % Add legend and labels
% xlabel('x_1');
% ylabel('x_2');
% legend('off')
% grid;
% pg = plot(2, 0,'r');
% pn = plot(2, 0,'k');
% pj = plot(2, 0,'b');
% % pf = plot(2, 0,'Color',[0 0.4 0]);
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 14; % Set font size of axis
% ax.YAxis.FontSize = 14;
% xlim([0.45 0.7])
% ylim([0.41 0.56])
% legend([pg pn pj],{'GoTube', 'NNV (ours)', 'Juliareach'},"Location","best",'FontSize',13);
% saveas(f1,'spiralL_compare_0.05_last.pdf');

% -----------------------------------------------------------------------------
% -----------------------------------------------------------------------------

% NNV
sp1_nnv = load("results/spiral_0.01.mat");

% JuliaReach
julia_path = "juliareach/results/";
sp1_julia = load(julia_path + "spiral4.mat");

figure;
hold on;
% plot_gotube(sp1_gt,1,2,2,'r'); 
Star.plotBoxes_2D_noFill(sp1_nnv.Rall(1:end-1),1,2,'k');
plot_juliareach(sp1_julia,1,2,'b');
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
legend('off')
grid;
% pg = plot(2, 0,'r');
% pn = plot(2, 0,'k');
% pj = plot(2, 0,'b');
ax = gca; % Get current axis
ax.XAxis.FontSize = 14; % Set font size of axis
ax.YAxis.FontSize = 14;
% pf = plot(2, 0,'Color',[0 0.4 0]);
% legend([pf pn pj],{'Flow*', 'NNV (ours)', 'Juliareach'},"Location","best",'FontSize',13);
run flowstar/results/spiralL1.m
legend('off')
hold on;
pn = plot(2, 0,'k');
pj = plot(2, 0,'b');
pf = plot(2, 0,'Color',[0 0.4 0]);
legend([pf pn pj],{'Flow*', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
saveas(gcf,'spiralL_compare_0.01.pdf');
pause(0.01);
xlim([0.23 0.49]);
ylim([0.48 0.7]);
saveas(gcf,'spiralL_compare_0.01_zoom.pdf');

% f1 = figure;
% hold on;
% plot_gotube(sp1_gt(:,end),1,2,2,'r'); 
% Star.plotBoxes_2D_noFill(sp1_nnv.Rall(end),1,2,'k');
% plot_juliareach(sp1_julia,1,2,'b',"last");
% % Add legend and labels
% xlabel('x_1');
% ylabel('x_2');
% legend('off')
% grid;
% pg = plot(2, 0,'r');
% pn = plot(2, 0,'k');
% pj = plot(2, 0,'b');
% % pf = plot(2, 0,'Color',[0 0.4 0]);
% xlim([0.52 0.62])
% ylim([0.46 0.51])
% ax = gca; % Get current axis
% ax.XAxis.FontSize = 14; % Set font size of axis
% ax.YAxis.FontSize = 14;
% legend([pg pn pj],{'GoTube', 'NNV (ours)', 'Juliareach'},"Location","best",'FontSize',13);
% saveas(f1,'spiral_compare_0.01_last.pdf');

% %% Run Flow* at the end
% figure;
% run flowstar/results/spiralL1.m
% figure;
% run flowstar/results/spiralL2.m
% figure;
% run flowstar/results/spiralL3.m