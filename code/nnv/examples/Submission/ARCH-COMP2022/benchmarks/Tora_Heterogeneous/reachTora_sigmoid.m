%% Reachability analysis of TORA (benchmark 9)
% Load components and set reachability parameters
% net = Load_nn('nn_tora_sigmoid.txt');
net = Load_nn('nn_tora_sigmoid.mat');
net.Layers(1).f = 'logsig';
net.Layers(2).f = 'logsig';
net.Layers(3).f = 'logsig';
% net.Layers(4).f = 'logsig'; % purelin

reachStep = 0.01;
controlPeriod = 0.5;
plant = NonLinearODE(4,1,@dynamicsTora, reachStep, controlPeriod, eye(4));
% noise = Star(-0.0001, 0.0001);
% plant.set_taylorTerms(10);
% plant.set_zonotopeOrder(100);
% error = 0.01;
% plant.options.maxError = [error; error; error; error];
time = 0:controlPeriod:10;
steps = length(time);
% Initial set
% lb = [-0.77; -0.45; 0.51; -0.3];
% ub = [-0.75; -0.43; 0.54; -0.28];
lb = [-0.751; -0.45; 0.535; -0.3];
ub = [-0.75; -0.449; 0.54; -0.295];
% init_set = Box(lb,ub);
offset = 0;
scale_factor = 11;

%% Reachability analysis
% init_set = Box(lb,ub);
init_set = Star(lb,ub);
% init = init_set.partition([1 2 3 4],[20 20 6 4]);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
% reachAll = cell(length(init),1);
% Execute reachabilty analysis
% for i =1:steps
t = tic;
% for j = 1:length(init)
% for j = 1:30
%     init_set = init(j).toStar;
    reachSub = init_set;
    for i = 1:10
        % Compute controller output set
        input_set = net.reach(init_set,'approx-star',1);
        input_set = input_set.affineMap(scale_factor,-offset);
        % Compute plant reachable set
        init_set = plant.stepReachStar(init_set, input_set,'lin');
        reachSub = [reachSub init_set];
    end
%     reachAll{j} = reachSub;
% end
% init_set = plant.stepReachStar(init_set, input_set,'lin');
% reachAll = [reachAll init_set];
t = toc(t);

path_out = [path_results(), filesep, 'Tora_Heterogeneous', filesep];
mkdir(path_out)
save([path_out, 'reach_sigmoid.mat'],'t','reachAll','-v7.3');

%% Visualize results
f = figure;
rectangle('Position',[-0.1,-0.9,0.3,0.3],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
% for r=1:length(init)
%     Star.plotBoxes_2D_noFill(reachAll{r},1,2,'b');
% end
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach_sigmoid.pdf']);

f = figure;
rectangle('Position',[-0.1,-0.9,0.3,0.3],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
hold on;
for r=1:length(init)
    Star.plotBoxes_2D_noFill(reachAll{r}(end),1,2,'b');
end
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach_sigmoid_last.pdf']);
