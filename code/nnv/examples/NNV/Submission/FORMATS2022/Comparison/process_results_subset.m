%% Process and plot all reachability results
% Load all reachability results and generate plots for paper
%% 1 - Load JuliaReach
if is_codeocean
    julia_path = "/results/logs/juliaresults/";
else
    julia_path = "juliareach/results/";
end
% FPA
fpa_julia_2 = load(julia_path+'fpa2.mat');

%% 2 - Load Gotube
if is_codeocean
    gotube_path = '/results/logs/tuberesults/'
else
    gotube_path = 'Gotube/results/'
end
% FPA
fpa1_gotube = "fpaCTRNN_10.0_0.01_1000_0.01_0.01_1.5_GoTube.txt";
fpa1_gotube = gotube_path+fpa1_gotube; 
fileID = fopen(fpa1_gotube,'r');
fpa_GT = fscanf(fileID, "%f", [32, 1001]);
fclose(fileID);


%% 3 - Load NNV results
nnv_path = "nnvresults/";
% Damped Forced Pendulum
fpa_nnv = "fpa_reach.mat";
fpa_nnv = load(nnv_path+fpa_nnv);


%% 5 - Plot FPA reach sets
f_fpa = figure;
hold on;
plot_gotube(fpa_GT,1,2,5,'r'); 
Star.plotBoxes_2D_noFill(fpa_nnv.Rall,1,2,'k');
plot_juliareach(fpa_julia_2,1,2,'b');
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
if is_codeocean
    saveas(f_fpa,'/results/logs/fpa_compare.pdf');
else
    saveas(f_fpa,'fpa_compare.pdf');
end
pause(0.1); % Make sure figure is saved before modifying axis (zooming in)
xlim([-1.5 -1.3])
ylim([-2.14 -1.98])
if is_codeocean
    saveas(f_fpa,'/results/logs/fpa_compare_zoom.pdf');
else
    saveas(f_fpa,'fpa_compare_zoom.pdf');
end