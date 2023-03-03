%% Reachability analysis of ACC
% Load components and set reachability parameters
net = Load_nn('controller_5_20.mat');
reachStep = 0.02;
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

%% Visualize results
t_gap = 1.4;
D_default = 10;
outAll = [];
safe_dis = [];
for i=1:length(plant.intermediate_reachSet)
    outAll = [outAll plant.intermediate_reachSet(i).affineMap(output_mat,[])];
    safe_dis = [safe_dis plant.intermediate_reachSet(i).affineMap([0 0 0 0 t_gap 0], D_default)];
end
times = reachStep:reachStep:tF;
mkdir([path_results(), filesep, 'ACC'])
save([path_results(), filesep, '/ACC/reach.mat'],'R','rT','-v7.3');
f = figure;
Star.plotRanges_2D(outAll,2,times,'b');
hold on;
Star.plotRanges_2D(safe_dis,1,times,'r');
xlabel('Time (s)');
ylabel('Distance (m)')
saveas(f,[path_results(), filesep, '/ACC/reach.pdf']);



