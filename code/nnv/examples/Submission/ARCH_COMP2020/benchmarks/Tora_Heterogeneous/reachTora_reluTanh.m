%% Reachability analysis of TORA (benchmark 9)
% Load components and set reachability parameters
net = Load_nn('nn_tora_relu_tanh.txt');
net.Layers(4).f = 'tansig';
reachStep = 0.003;
controlPeriod = 0.5;
% controlPeriod = 1;
plant = NonLinearODE(4,1,@dynamicsTora, reachStep, controlPeriod, eye(4));
% noise = Star(-0.0001, 0.0001);
plant.set_taylorTerms(10);
plant.set_zonotopeOrder(200);
% plant.set_polytopeOrder(50);% error = 0.001;
error = 0.0005;
plant.options.maxError = [error; error; error; error];
time = 0:controlPeriod:20;
steps = length(time);
% Initial set
lb = [-0.77; -0.45; 0.51; -0.3];
% ub = [-0.75; -0.43; 0.54; -0.28];
ub = [-0.76; -0.44; 0.52; -0.29];
init_set = Box(lb,ub);
goal = Box([-0.1;-0.9],[0.2;-0.6]);
offset = 0;
scale_factor = 11;

%% Reachability analysis
init_set = Star(lb,ub);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
% Execute reachabilty analysis
% for i =1:steps
t = tic;
for i = 1:9
    % Compute plant reachable set
    init_set = plant.stepReachStar(init_set(1), input_set(1));
    reachAll = [reachAll init_set(1)];
    % Compute controller output set
    input_set = net.reach(init_set(1),'approx-star',1);
    input_set = input_set(1).affineMap(scale_factor,-offset);
end
init_set = plant.stepReachStar(init_set(1), input_set(1));
reachAll = [reachAll init_set(1)];
timing = toc(t);
%% Set output path
path_out_t = ['..' filesep path_results() filesep 'tora_heterogeneous' filesep];
mkdir(path_out_t);
save([path_out_t 'sets_reluTanh'],'reachAll','timing','-v7.3');

%% Visualize results
f = figure('visible','off');
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachAll,1,2,'m');
grid;
Box.plotBoxes_2D(goal,1,2,'r');
grid;
title('Reachable sets for dimensions 1 and 2')
xlabel('x1');
ylabel('x2');
saveas(f,[path_out_t 'reluTanh_plot.jpg']);