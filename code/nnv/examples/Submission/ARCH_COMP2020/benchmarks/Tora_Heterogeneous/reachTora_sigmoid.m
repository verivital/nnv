%% Reachability analysis of TORA (benchmark 9)
% Load components and set reachability parameters
% net = Load_nn('nn_tora_sigmoid.txt');
net = Load_nn('nn_tora_sigmoid.mat');
net.Layers(1).f = 'logsig';
net.Layers(2).f = 'logsig';
net.Layers(3).f = 'logsig';
net.Layers(4).f = 'logsig';

reachStep = 0.01;
controlPeriod = 0.5;
plant = NonLinearODE(4,1,@dynamicsTora, reachStep, controlPeriod, eye(4));
% noise = Star(-0.0001, 0.0001);
plant.set_taylorTerms(10);
plant.set_zonotopeOrder(100);
% plant.set_polytopeOrder(5);% error = 0.001;
error = 0.01;
plant.options.maxError = [error; error; error; error];
time = 0:controlPeriod:20;
steps = length(time);
% Initial set
lb = [-0.77; -0.45; 0.51; -0.3];
ub = [-0.75; -0.43; 0.54; -0.28];
% ub = [0.61; -0.69; -0.39; 0.51];
% init_set = Box(lb,ub);
offset = 0;
scale_factor = 11;

%% Reachability analysis
B = Box(lb, ub);
Bs = B.partition([1 2 3 4], [8 8 8 8]); % partition the box into smaller boxes at index 1, 2, 3, 4. At each index, we split one time

n = length(Bs);
IS = [];
for i=1:n
    IS = [IS Bs(i).toStar]; 
end

init_set = IS(1);
%init_set = Star(lb,lb);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
% Execute reachabilty analysis
% for i =1:steps
t = tic;
for i =1:10
    % Compute plant reachable set
    init_set = plant.stepReachStar(init_set(1), input_set(1));
    reachAll = [reachAll init_set(1)];
    % Compute controller output set
    input_set = net.reach(init_set(1),'approx-star');
    input_set = input_set(1).affineMap(scale_factor,-offset);
end
timing = toc(t);
%% Set output path
path_out_t = ['..' filesep path_results() filesep 'tora_heterogeneous' filesep];
mkdir(path_out_t);
save([path_out_t 'sets_sigmoid'],'reachAll','timing','-v7.3');

%% Visualize results
goal = Box([-0.1;-0.9],[0.2;-0.6]);
f = figure('visible','off');
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachAll,1,2,'m');
grid;
hold on;
goal.plot;
title('Reachable sets for dimensions 1 and 2')
xlabel('x1');
ylabel('x2');
saveas(f,[path_out_t 'sigmoid_plot.jpg']);