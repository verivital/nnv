%% Reachability analysis of TORA (benchmark 9)
% Load components and set reachability parameters
% net = Load_nn('nn_tora_relu_tanh.txt');
net = Load_nn('nn_tora_relu_tanh.mat');
net.Layers(4).f = 'tansig';

reachStep = 0.05;
controlPeriod = 0.5;
% controlPeriod = 1;
plant = NonLinearODE(4,1,@dynamicsTora, reachStep, controlPeriod, eye(4));
% noise = Star(-0.0001, 0.0001);
% plant.set_taylorTerms(8);
% plant.set_zonotopeOrder(20);
% % plant.set_polytopeOrder(50);% error = 0.001;
% error = 0.0005;
% plant.options.maxError = [error; error; error; error];
time = 0:controlPeriod:20;
steps = length(time);
% Initial set
lb = [-0.77; -0.45; 0.51; -0.3];
ub = [-0.75; -0.43; 0.54; -0.28];
% lb = [-0.77; -0.45; 0.535; -0.3];
% ub = [-0.765; -0.4475; 0.54; -0.295];
goal = Box([-0.1;-0.9],[0.2;-0.6]);
offset = 0;
scale_factor = 11;

%% Reachability analysis
init_set = Box(lb,ub);
% init_set = Star(lb,ub);
init = init_set.partition([1 2 3 4],[4 8 6 4]);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
reachAll = cell(length(init),1);
% Execute reachabilty analysis
% for i =1:steps
t = tic;
for j = 1:length(init)
    init_set = init(j).toStar;
    reachSub = init_set;
    for i = 1:10
        % Compute controller output set
        input_set = net.reach(init_set,'approx-star',1);
        input_set = input_set.affineMap(scale_factor,-offset);
        % Compute plant reachable set
        init_set = plant.stepReachStar(init_set, input_set,'lin');
        reachSub = [reachSub init_set];
    end
    reachAll{j} = reachSub;
end
% init_set = plant.stepReachStar(init_set, input_set,'lin');
% reachAll = [reachAll init_set];
t = toc(t);

path_out = [path_results(), filesep, 'Tora_Heterogeneous', filesep];
mkdir(path_out)
save([path_out, 'reach_reluTanh.mat'],'t','reachAll','-v7.3');
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
saveas(f,[path_out, 'reach_reluTanh.pdf']);

f = figure;
rectangle('Position',[-0.1,-0.9,0.3,0.3],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
hold on;
for r=1:length(init)
    Star.plotBoxes_2D_noFill(reachAll{r}(end),1,2,'b');
end
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach_reluTanh_last.pdf']);