%% Reachability analysis of TORA (benchmark 9)
% Load components and set reachability parameters
% net = Load_nn('controllerTora.onnx');
net = Load_nn('controllerTora_nnv.mat');
% reachStep = 0.005;
reachStep = 0.01;
controlPeriod = 1;
plant = NonLinearODE(4,1,@dynamics9, reachStep, controlPeriod, eye(4));
% noise = Star(-0.0001, 0.0001);
% plant.set_taylorTerms(10);
% plant.set_zonotopeOrder(100);
% plant.set_polytopeOrder(5);% error = 0.001;
% error = 0.01;
% plant.options.maxError = [error; error; error; error];
time = 0:controlPeriod:20;
steps = length(time);
% Initial set
lb = [0.69; -0.7; -0.4; 0.5];
% ub = [0.7; -0.6; -0.3; 0.6];
% ub = [0.61; -0.69; -0.39; 0.51];
ub = [0.7; -0.69; -0.39; 0.51];
% init_set = Box(lb,ub);
offset = 10;
scale_factor = 1;

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
for i =1:6
    % Compute plant reachable set
    init_set = plantReach(plant,init_set,input_set);
%     init_set = plant.stepReachStar(init_set, input_set);
    reachAll = [reachAll init_set];
    % Compute controller output set
    inp_set = net.reach(init_set,'approx-star');
    input_set = [];
    for k = 1:length(inp_set)
        input_set = [input_set inp_set(k).affineMap(1,-offset)];
    end
end
disp(' ');
timing = toc(t);
%% Set output path
path_out_b9 = ['..' filesep path_results() filesep 'benchmark9' filesep];
mkdir(path_out_b9);
save([path_out_b9 'sets'],'reachAll','timing','-v7.3')

%% Visualize results
f = figure('visible','off');
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachAll,1,2,'m');
grid;
title('Reachable sets Benchmark 9 - TORA')
xlabel('x1');
ylabel('x3');
saveas(f,[path_out_b9 'plot1v3.jpg']);

f1 = figure('visible','off');
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,2,4,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachAll,2,4,'m');
grid;
title('Reachable sets Benchmark 9 - TORA')
xlabel('x2');
ylabel('x4');
saveas(f1,[path_out_b9 'plot3v4.jpg']);

%% Helper function
function init_set = plantReach(plant,init_set,input_set)
    nS = length(init_set);
    nL = length(input_set);
    ss = [];
    for k=1:nS
        for l=1:nL
            ss =[ss plant.stepReachStar(init_set(k), input_set(l))];
        end
    end
    init_set = ss;
end