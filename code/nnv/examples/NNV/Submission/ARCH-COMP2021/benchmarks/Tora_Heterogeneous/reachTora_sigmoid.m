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
% plant.set_taylorTerms(10);
% plant.set_zonotopeOrder(100);
% error = 0.01;
% plant.options.maxError = [error; error; error; error];
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
init_set = Star(lb,ub);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
reachLin = init_set;
% Execute reachabilty analysis
t = tic;
for i =1:10
    % Compute plant reachable set
    init_set = plant.stepReachStar(init_set, input_set,'lin');
    reachLin = [reachLin init_set];
    % Compute controller output set
    input_set = net.reach(init_set,'approx-star');
    input_set = input_set.affineMap(scale_factor,-offset);
end
tlin = toc(t);

% Reach 2
plant2 = NonLinearODE(4,1,@dynamicsTora, reachStep, controlPeriod, eye(4));
lb = [-0.77; -0.45; 0.51; -0.3];
ub = [-0.75; -0.43; 0.54; -0.28];
init_set = Star(lb,ub);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
reachPoly = init_set;
% Execute reachabilty analysis
t = tic;
for i =1:10
    % Compute plant reachable set
    init_set = plant2.stepReachStar(init_set, input_set,'poly');
    reachPoly = [reachPoly init_set];
    % Compute controller output set
    input_set = net.reach(init_set,'approx-star');
    input_set = input_set.affineMap(scale_factor,-offset);
end
tpoly = toc(t);
path_out = [path_results(), filesep, 'Tora_Heterogeneous', filesep];
mkdir(path_out)
save([path_out, 'reach_sigmoid.mat'],'tlin','tpoly','reachPoly','reachLin','-v7.3');

%% Visualize results
goal = Box([-0/1;-0.9],[0.2;-0.6]);
f = figure;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachLin,1,2,'m');
grid;
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach_sigmoid.pdf']);

f = figure;
Star.plotBoxes_2D_noFill(plant2.intermediate_reachSet,1,2,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachPoly,1,2,'m');
grid;
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach_sigmoid_poly.pdf']);