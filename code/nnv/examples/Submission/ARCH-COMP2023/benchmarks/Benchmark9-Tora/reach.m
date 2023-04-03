%% Reachability analysis of TORA (benchmark 9)
% Load components and set reachability parameters
% net = Load_nn('controllerTora.onnx');
net = load_NN_from_mat('controllerTora.mat');
% net = onnx2nnv('controllerTora.onnx');
reachStep = 0.01;
controlPeriod = 1;
plant = NonLinearODE(4,1,@dynamics9, reachStep, controlPeriod, eye(4));
time = 0:controlPeriod:20;
steps = length(time);
% Initial set
lb = [0.69; -0.7; -0.4; 0.5];
% lb = [0.6; -0.7; -0.4; 0.5];
% ub = [0.7; -0.6; -0.3; 0.6];
% ub = [0.61; -0.69; -0.39; 0.51];
ub = [0.7; -0.69; -0.39; 0.51];
% init_set = Box(lb,ub);
offset = 10; % Applied this in the dynamics9.m function
scale_factor = 1;

% Custom plant reachability options
plant.options.taylorTerms = 4;
plant.options.zonotopeOrder = 200;
alg = 'poly';
plant.options.alg = alg;
plant.options.tensorOrder = 3;
plant.options.errorOrder = 10;
plant.options.intermediateOrder = 50;

%% Reachability analysis
init_set = Star(lb,ub);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
% Execute reachabilty analysis
reachOptions.reachMethod = 'approx-star';
tFinal = 4;
t = tic;
for i = 1:tFinal % CORA verifies computing U beforehand
    % Compute controller output set
    input_set = net.reach(init_set,reachOptions);
    % Compute plant reachable set
%     if length(input_set) > 1
%         input_set = Star.get_hypercube_hull(input_set);
%         input_set = input_set.toStar();
%     end
    init_set = plant.stepReachStar(init_set, input_set,alg);
    reachAll = [reachAll init_set];
end
tlin = toc(t);

% Compute using NNCS
% nncs = NonlinearNNCS(net,plant);
% reachPRM.numSteps = 5;
% reachPRM.init_set = init_set;
% reachPRM.numCores = 1;
% reachPRM.reachMethod = 'approx-star';
% [R,rT] = nncs.reach(reachPRM);

% Save results
% if is_codeocean
%     save('/results/logs/tora.mat', 'reachAll','t','-v7.3');
% else
%     save('tora.mat', 'reachAll','t','-v7.3');
% end

%% Visualize results
plant.get_interval_sets;

f = figure;
grid;
hold on;
rectangle('Position',[-2,-2,4,4],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
% Star.plotBoxes_2D_noFill(reachAll,1,2,'m');
grid;
xlabel('x1');
ylabel('x2');
% saveas(f,[path_out, 'reach1vs3.pdf']);

f1 = figure;
grid;
hold on;
rectangle('Position',[-2,-2,4,4],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
% Star.plotBoxes_2D_noFill(reachAll,3,4,'m');
grid;
xlabel('x2');
ylabel('x4');
% saveas(f1,[path_out, 'reach3vs4.pdf']);

%% Helper function
% function init_set = plantReach(plant,init_set,input_set,algoC)
%     nS = length(init_set);
%     nL = length(input_set);
%     ss = [];
%     for k=1:nS
%         for l=1:nL
%             ss =[ss plant.stepReachStar(init_set(k), input_set(l),algoC)];
%         end
%     end
%     init_set = ss;
% end