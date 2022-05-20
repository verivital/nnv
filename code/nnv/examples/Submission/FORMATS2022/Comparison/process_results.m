%% Process and plot all reachability results
% Load all reachability results and generate plots for paper
%% 1 - Load JuliaReach
julia_path = "juliareach/results/";
% Cartpole
% cartpole_julia_1 = load(julia_path + "cartpole1.mat");
cartpole_julia_2 = load(julia_path + "cartpole2.mat");
% % Damped Forced Pendulum
% dfp_julia_1 = load(julia_path+'dfp1.mat');
dfp_julia_2 = load(julia_path+'dfp2.mat');

%% 2 - Load Gotube
gotube_path = "/home/manzand/Documents/Python/GoTube/saved_outputs/";
% Damped Forced Pendulum
dfp1_gotube = "dfpCTRNN_10.0_0.01_1000_0.01_0.01_1.5_GoTube.txt";
dfp1_gotube = gotube_path+dfp1_gotube; 
fileID = fopen(dfp1_gotube,'r');
dfp_GT = fscanf(fileID, "%f", [32, 1001]);
fclose(fileID);
% Cartpole
cartpole_gotube = "cartpoleCTRNN_2.0_0.02_10000_0.0001_0.01_1.5_GoTube.txt";
cartpole_gotube = gotube_path+cartpole_gotube;
fileID= fopen(cartpole_gotube,'r');
% cp_GT = fscanf(fileID, "%f");
cp_GT = fscanf(fileID, "%f", [158, 101]);
fclose(fileID);


%% 3 - Load NNV results
nnv_path = "results/";
% Damped Forced Pendulum
dfp_nnv = "dfp_reach.mat";
dfp_nnv = load(nnv_path+dfp_nnv);
% Cartpole
cartpole_nnv = "cartpole_reach.mat";
cartpole_nnv = load(nnv_path+cartpole_nnv);


%% 5 - Plot Damped Forced Pendulum reach sets
f_dfp = figure;
hold on;
plot_gotube(dfp_GT,1,2,5,'r'); 
Star.plotBoxes_2D_noFill(dfp_nnv.Rall,1,2,'k');
plot_juliareach(dfp_julia_2,1,2,'b');
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
% Add initial state to generate legends x0 = [0.21535, -0.58587]
legend('off')
grid;
pg = plot(0.21535, -0.58587,'r');
pn = plot(0.21535, -0.58587,'k');
pj = plot(0.21535, -0.58587,'b');
ax = gca; % Get current axis
ax.XAxis.FontSize = 15; % Set font size of axis
ax.YAxis.FontSize = 15;
legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',14);
saveas(f_dfp,'dfp_compare.pdf');
pause(0.1); % Make sure figure is saved before modifying axis (zooming in)
xlim([-1.5 -1.3])
ylim([-2.14 -1.98])
% legend('off');
% legend([pg pn pj],{'GoTube', 'NNV (ours)', 'Juliareach'},"Location","best");
saveas(f_dfp,'dfp_compare_zoom.pdf');


%% 6 - Plot Cartpole reachsets
f_cp = figure;
hold on;
plot_gotube(cp_GT,1,2,5,'r'); 
Star.plotBoxes_2D_noFill(cartpole_nnv.Rall,1,2,'k');
plot_juliareach(cartpole_julia_2,1,2,'b');
% Add legend and labels
xlabel('x_1');
ylabel('x_2');
legend('off');
grid;
pg = plot(0,0,'r');
pn = plot(0,0,'k');
pj = plot(0,0,'b');
ax = gca; % Get current axis
ax.XAxis.FontSize = 15; % Set font size of axis
ax.YAxis.FontSize = 15;
legend([pg pn pj],{'GoTube', 'NNVODE (ours)', 'Juliareach'},"Location","best",'FontSize',14);
saveas(f_cp,'cartpole_compare.pdf');
pause(0.1); % Make sure figure is saved before modifying axis (zooming in)
xlim([0.6 0.68])
ylim([0.27 0.37])
% legend('off');
% legend([pg pn pj],{'GoTube', 'NNV (ours)', 'Juliareach'},"Location","best");
saveas(f_cp,'cartpole_compare_zoom.pdf');

%% 7 - Plot Damped Oscillator reachsets
zonoF0 = load("results/DampedOsc_ilnodenl_true_0_zonoF.mat");
zonoA0 = load("results/DampedOsc_ilnodenl_true_0_zonoA.mat");
polyF0 = load("results/DampedOsc_ilnodenl_true_0_polyF.mat");
polyA0 = load("results/DampedOsc_ilnodenl_true_0_polyA.mat");
zonoF1 = load("results/DampedOsc_ilnodenl_true_1_zonoF.mat");
zonoA1 = load("results/DampedOsc_ilnodenl_true_1_zonoA.mat");
polyF1 = load("results/DampedOsc_ilnodenl_true_1_polyF.mat");
polyA1 = load("results/DampedOsc_ilnodenl_true_1_polyA.mat");

fosc0 = figure;
hold on;
% Star.plotBoxes_2D_noFill(zonoF0.Rall(end),1,2,'k');
% Star.plotBoxes_2D_noFill(zonoA0.Rall(end),1,2,'r');
% Star.plotBoxes_2D_noFill(polyF0.Rall(end),1,2,'b');
% Star.plotBoxes_2D_noFill(polyA0.Rall(end),1,2,'g');
Star.plots(zonoF0.Rall(end),'k');
Star.plots(polyF0.Rall(end),'b');
Star.plots(zonoA0.Rall(end),'r');
Star.plots(polyA0.Rall(end),'g');
xlabel('x_1');
ylabel('x_2');
legend('off');
grid;
x0 = [-1.4996; -0.4609];
p1 = plot(x0(1), x0(2),'k');
p2 = plot(x0(1), x0(2),'r');
p3 = plot(x0(1), x0(2),'b');
p4 = plot(x0(1), x0(2),'g');
ax = gca; % Get current axis
ax.XAxis.FontSize = 15; % Set font size of axis
ax.YAxis.FontSize = 15;
legend([p1 p2 p3 p4],{'ZonoF','ZonoA', 'PolyF', 'PolyA'},"Location","best",'FontSize',14);
saveas(fosc0,'dampedOsc_compare0.pdf');
pause(0.1);
xlim([-1.1891 -1.18898])
ylim([0.93732 0.93759])
saveas(fosc0,'dampedOsc_compare0_zoomed.pdf');

fosc1 = figure;
hold on;
% Star.plotBoxes_2D_noFill(zonoF0.Rall(end),1,2,'k');
% Star.plotBoxes_2D_noFill(zonoA0.Rall(end),1,2,'r');
% Star.plotBoxes_2D_noFill(polyF0.Rall(end),1,2,'b');
% Star.plotBoxes_2D_noFill(polyA0.Rall(end),1,2,'g');
% Star.plots(zonoF1.Rall(end),'k');
Star.plots(zonoA1.Rall(end),'r');
Star.plots(polyF1.Rall(end),'b');
Star.plots(zonoF1.Rall(end),'k');
Star.plots(polyA1.Rall(end),'g');
xlabel('x_1');
ylabel('x_2');
legend('off');
grid;
x0 = [-1.4996; -0.4609];
p1 = plot(x0(1), x0(2),'k');
p2 = plot(x0(1), x0(2),'r');
p3 = plot(x0(1), x0(2),'b');
p4 = plot(x0(1), x0(2),'g');
ax = gca; % Get current axis
ax.XAxis.FontSize = 15; % Set font size of axis
ax.YAxis.FontSize = 15;
legend([p1 p2 p3 p4],{'ZonoF','ZonoA', 'PolyF', 'PolyA'},"Location","best",'FontSize',14);
saveas(fosc1,'dampedOsc_compare1.pdf');
pause(0.1);
xlim([-1.17825 -1.17815])
ylim([0.98400 0.98414])
saveas(fosc1,'dampedOsc_compare1_zoomed.pdf');