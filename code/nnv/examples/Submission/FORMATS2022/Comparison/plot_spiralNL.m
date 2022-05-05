%% Plot one scenario for the nonlinear spiral 2D


%% Gotube
gotube_path = "/home/manzand/Documents/Python/GoTube/saved_outputs/";
% Nonlinear spiral (case 1)
sp1_gt = "spiralNL_10.0_0.01_1000_0.01_0.01_1.5_GoTube.txt";
sp1_gt = gotube_path+sp1_gt; 
fileID = fopen(sp1_gt,'r');
sp1_gt = fscanf(fileID, "%f", [8 1001]);
fclose(fileID);
% Nonlinear spiral (case 2)
sp2_gt = "spiralNL_10.0_0.01_1000_0.05_0.01_1.5_GoTube.txt";
sp2_gt = gotube_path+sp2_gt; 
fileID = fopen(sp2_gt,'r');
sp2_gt = fscanf(fileID, "%f", [8 1001]);
fclose(fileID);
% Nonlinear spiral (case 3)
sp3_gt = "spiralNL_10.0_0.01_1000_0.1_0.01_1.5_GoTube.txt";
sp3_gt = gotube_path+sp3_gt; 
fileID = fopen(sp3_gt,'r');
sp3_gt = fscanf(fileID, "%f", [8 1001]);
fclose(fileID);

%% NNV
sp1_nnv = load("results/spiral_nl_0.01.mat");
sp2_nnv = load("results/spiral_nl_0.05.mat");
sp3_nnv = load("results/spiral_nl_0.1.mat");

%% JuliaReach
julia_path = "juliareach/results/";

sp1_julia = load(julia_path + "spiral1.mat");
sp2_julia = load(julia_path + "spiral2.mat");
sp3_julia = load(julia_path + "spiral3.mat");


%% Generate Plots

f1 = figure;
hold on;
plot_gotube(sp3_gt,1,2,2,'r'); 
Star.plotBoxes_2D_noFill(sp3_nnv.Rall(1:end-1),1,2,'k');
plot_juliareach(sp3_julia,1,2,'b');
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
legend('off')
grid;
pg = plot(2, 0,'r');
pn = plot(2, 0,'k');
pj = plot(2, 0,'b');
% pf = plot(2, 0,'Color',[0 0.4 0]);
ax = gca; % Get current axis
ax.XAxis.FontSize = 14; % Set font size of axis
ax.YAxis.FontSize = 14;
legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
saveas(f1,'spiral_compare_0.1_nl.pdf');

f1 = figure;
hold on;
plot_gotube(sp3_gt(:,end),1,2,2,'r'); 
Star.plotBoxes_2D_noFill(sp3_nnv.Rall(end-1),1,2,'k');
plot_juliareach(sp3_julia,1,2,'b',"last");
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
legend('off')
grid;
pg = plot(2, 0,'r');
pn = plot(2, 0,'k');
pj = plot(2, 0,'b');
% pf = plot(2, 0,'Color',[0 0.4 0]);
ax = gca; % Get current axis
ax.XAxis.FontSize = 14; % Set font size of axis
ax.YAxis.FontSize = 14;
xlim([0.2 0.9])
ylim([0.2 0.8])
legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
saveas(f1,'spiral_compare_0.1_nl_last.pdf');

f1 = figure;
hold on;
plot_gotube(sp2_gt,1,2,2,'r'); 
Star.plotBoxes_2D_noFill(sp2_nnv.Rall(1:end-1),1,2,'k');
plot_juliareach(sp2_julia,1,2,'b');
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
legend('off')
grid;
pg = plot(2, 0,'r');
pn = plot(2, 0,'k');
pj = plot(2, 0,'b');
% pf = plot(2, 0,'Color',[0 0.4 0]);
ax = gca; % Get current axis
ax.XAxis.FontSize = 14; % Set font size of axis
ax.YAxis.FontSize = 14;
legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
saveas(f1,'spiral_compare_0.05_nl.pdf');

f1 = figure;
hold on;
plot_gotube(sp2_gt(:,end),1,2,2,'r'); 
Star.plotBoxes_2D_noFill(sp2_nnv.Rall(end-1),1,2,'k');
plot_juliareach(sp2_julia,1,2,'b',"last");
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
legend('off')
grid;
pg = plot(2, 0,'r');
pn = plot(2, 0,'k');
pj = plot(2, 0,'b');
% pf = plot(2, 0,'Color',[0 0.4 0]);
ax = gca; % Get current axis
ax.XAxis.FontSize = 14; % Set font size of axis
ax.YAxis.FontSize = 14;
xlim([0.45 0.7])
ylim([0.41 0.56])
legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
saveas(f1,'spiral_compare_0.05_nl_last.pdf');

f1 = figure;
hold on;
plot_gotube(sp1_gt,1,2,2,'r'); 
Star.plotBoxes_2D_noFill(sp1_nnv.Rall(1:end-1),1,2,'k');
plot_juliareach(sp1_julia,1,2,'b');
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
legend('off')
grid;
pg = plot(2, 0,'r');
pn = plot(2, 0,'k');
pj = plot(2, 0,'b');
ax = gca; % Get current axis
ax.XAxis.FontSize = 14; % Set font size of axis
ax.YAxis.FontSize = 14;
% pf = plot(2, 0,'Color',[0 0.4 0]);
legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
saveas(f1,'spiral_compare_nl_0.01.pdf');

f1 = figure;
hold on;
plot_gotube(sp1_gt(:,end),1,2,2,'r'); 
Star.plotBoxes_2D_noFill(sp1_nnv.Rall(end),1,2,'k');
plot_juliareach(sp1_julia,1,2,'b',"last");
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
legend('off')
grid;
pg = plot(2, 0,'r');
pn = plot(2, 0,'k');
pj = plot(2, 0,'b');
% pf = plot(2, 0,'Color',[0 0.4 0]);
xlim([0.52 0.62])
ylim([0.46 0.51])
ax = gca; % Get current axis
ax.XAxis.FontSize = 14; % Set font size of axis
ax.YAxis.FontSize = 14;
legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',13);
saveas(f1,'spiral_compare_0.01_nl_last.pdf');