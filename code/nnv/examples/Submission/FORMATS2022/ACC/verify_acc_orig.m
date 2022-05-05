%% Reachability analysis of ACC 
% Taken from arch-comp 2021
% Load components and set reachability parameters
net = Load_nn('controller_5_20.mat');
reachStep = 0.01;
controlPeriod = 0.1;
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
plant = NonLinearODE(6,1,@dynamicsACC, reachStep, controlPeriod, output_mat);

%% Reachability analysis
tF = 5; % seconds
time = 0:controlPeriod:5;
steps = length(time);
input_ref = [30;1.4];
% Initial set
lb = [90; 32; 0; 10; 30; 0];
ub = [110; 32.2; 0; 11; 30.2; 0];
init_set = Star(lb,ub);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
% Execute reachabilty analysis
nncs = NonlinearNNCS(net,plant);
reachPRM.ref_input = input_ref;
reachPRM.numSteps = 50;
reachPRM.init_set = init_set;
reachPRM.numCores = 1;
reachPRM.reachMethod = 'approx-star';
[R,rT] = nncs.reach(reachPRM);
disp("Time to compute ACC reach sets: " +string(rT));
save("results_original.mat","R","rT");

%% Visualize results
t_gap = 1.4;
D_default = 10;
outAll = [];
safe_dis = [];
% for i=1:length(plant.intermediate_reachSet)
%     outAll = [outAll plant.intermediate_reachSet(i).affineMap(output_mat,[])];
%     safe_dis = [safe_dis plant.intermediate_reachSet(i).affineMap([0 0 0 0 t_gap 0], D_default)];
% end
% times = reachStep:reachStep:tF;
for i=1:length(R)
    outAll = [outAll R(i).affineMap(output_mat,[])];
    safe_dis = [safe_dis R(i).affineMap([0 0 0 0 t_gap 0], D_default)];
end
times = 0:controlPeriod:tF;
f = figure;
hold on;
pb = plot(0,85,'m');
pr = plot(0,85,'r');
Star.plotRanges_2D(outAll,2,times,'r');
hold on;
Star.plotRanges_2D(safe_dis,1,times,'m');
ax = gca; % Get current axis
ax.XAxis.FontSize = 15; % Set font size of axis
ax.YAxis.FontSize = 15;
xlabel('Time (s)');
ylabel('Distance (m)')
legend([pr,pb],{'rel dist','safe dist'},"Location","best",'FontSize',14);
% saveas(f,'reach_orig.png');
exportgraphics(f,'reach_orig.pdf','ContentType','vector');
% Visualize all variables
% trajR = plant.intermediate_reachSet;
trajR = R;
f = figure;
subplot(2,3,1)
Star.plotRanges_2D(trajR,1,times,'r');
xlabel('Time (s)');
ylabel('xlead');
subplot(2,3,2)
Star.plotRanges_2D(trajR,2,times,'b');
xlabel('Time (s)')
ylabel('vlead')
subplot(2,3,3)
Star.plotRanges_2D(trajR,3,times,'b');
xlabel('Time (s)')
ylabel('alead')
subplot(2,3,4)
Star.plotRanges_2D(trajR,4,times,'b');
xlabel('Time (s)')
ylabel('xego')
subplot(2,3,5)
Star.plotRanges_2D(trajR,5,times,'b');
xlabel('Time (s)')
ylabel('vego')
subplot(2,3,6)
Star.plotRanges_2D(trajR,6,times,'b');
xlabel('Time (s)')
ylabel('aego')
saveas(f,'reach_orig_all.png');