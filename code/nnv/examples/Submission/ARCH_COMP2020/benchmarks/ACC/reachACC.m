%% Reachability analysis of ACC
% Load components and set reachability parameters
% net = Load_nn('controller_5_20.onnx');
net = Load_nn('controller_5_20_nnv.mat');
reachStep = 0.01;
controlPeriod = 0.1;
output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
plant = NonLinearODE(6,1,@dynamicsACC, reachStep, controlPeriod, output_mat);
% plant.set_taylorTerms(5);
% plant.set_zonotopeOrder(200);
% plant.set_polytopeOrder(20);
% error = 0.1;
% plant.options.maxError = [error; error; error; error; error; error];

%% Reachability analysis
tF = 5; % seconds
time = 0:controlPeriod:5;
steps = length(time);
input_ref = [30;1.4];
% Initial set
lb = [90; 32; 0; 10; 30; 0];
ub = [110; 32.2; 0; 0; 30.2; 11];
init_set = Star(lb,ub);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Create NNCS 
nncs = NonlinearNNCS(net,plant);
reachPRM.ref_input = input_ref;
reachPRM.numSteps = 50;
reachPRM.init_set = init_set;
reachPRM.numCores = 1;
reachPRM.reachMethod = 'approx-star';
[R,rT] = nncs.reach(reachPRM);

%% Set output path
path_out_acc = ['..' filesep path_results() filesep 'ACC' filesep];
mkdir(path_out_acc);
%% Visualize results
t_gap = 1.4;
D_default = 10;
outAll = [];
safe_dis = [];
for i=1:length(R)
    outAll = [outAll R(i).affineMap(output_mat,[])];
    safe_dis = [safe_dis R(i).affineMap([0 0 0 0 t_gap 0], D_default)];
end
save([path_out_acc 'sets'],'R','rT','outAll','safe_dis','-v7.3');
f = figure('visible','off');
Star.plotRanges_2D(outAll,2,time,'r');
hold on;
Star.plotRanges_2D(safe_dis,1,time,'b');
title('Relative distance (red) vs. Safe distance (blue)');
xlabel('Time (s)');
ylabel('Distance (m)')
% saveas(f,'../../results/ACC_distance_sets.jpg');
saveas(f, [path_out_acc 'plot_distance.jpg']);


