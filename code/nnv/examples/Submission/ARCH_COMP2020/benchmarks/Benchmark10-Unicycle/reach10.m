%% Reachability analysis of the Unicycle (benchmark 10)
% net = Load_nn('controllerB.onnx');
net = Load_nn('controllerB_nnv.mat');
controlPeriod = 0.2;
% controlPeriod = 0.5;
reachstep = 0.001;
plant = NonLinearODE(4,2,@dynamics10, reachstep, controlPeriod, eye(4));
% noise = Star(-0.0001, 0.0001);
% plant.set_taylorTerms(2);
% plant.set_zonotopeOrder(50);
% plant.set_polytopeOrder(20);
% error = 0.01;
% plant.options.maxError = [error; error; error; error];
tF = 10;
time = 0:controlPeriod:tF;
steps = length(time);
offset = 30;
offsetM = offset*ones(2,1);
scale_factor = 1;

%% Reachability analysis
% Initial set
lb = [9.5; -4.5; 2.1; 1.5];
ub = [9.51; -4.49; 2.11; 1.51];
% ub = [9.55; -4.45; 2.11; 1.51];
init_set = Star(lb,ub);
% Input set
lb = [0;0];
ub = [0;0];
input_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
% Execute reachabilty analysis
steps = 10;
nI = 1;
t = tic;
for i=1:steps
    % Compute plant reachable set
%     init_set = plant.stepReachStar(init_set, input_set);
    init_set = plantReach(plant,init_set(nI),input_set(1));
    nI = length(init_set);
    reachAll = [reachAll init_set(nI)];
    % Compute controller output set
    inp_set = net.reach(init_set(nI),'approx-star');
    input_set = [];
    input_set = inp_set(1).affineMap(eye(2),-offsetM);
%     for k = 1:length(inp_set)
%         input_set = [input_set inp_set(k).affineMap(eye(2),-offsetM)];
%     end
end
timing = toc(t);

%% Set output path
path_out_b10 = ['..' filesep path_results() filesep 'benchmark10' filesep];
mkdir(path_out_b10);
save([path_out_b10 'sets'],'reachAll','timing','-v7.3')

%% Visualize results
f = figure('visible','off');
% Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
hold on;
Star.plotBoxes_2D_noFill(reachAll,1,2,'b');
grid;
title('Benchmark 10 - Unicycle');
xlabel('x1');
ylabel('x2');
saveas(f,[path_out_b10 'plot1v2.jpg']);

f1 = figure('visible','off');
% Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
hold on;
Star.plotBoxes_2D_noFill(reachAll,3,4,'b');
grid;
title('Benchmark 10 - Unicycle');
xlabel('x1');
ylabel('x2');
saveas(f1,[path_out_b10 'plot3v4.jpg']);

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