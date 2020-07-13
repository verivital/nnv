%% Reachability analysis of Single Pendulum Benchmark

%% load the controller
net = Load_nn('controller_single_pendulum.mat');

% Specify the reach step, has to be smaller than the control period
reachStep = 0.001;
%% specify the control period as specified by the benchmark description
controlPeriod = 0.05;

% define the plant as specified by nnv
plant = NonLinearODE(3,1,@dynamics_sp, reachStep, controlPeriod, eye(3));
plant.set_zonotopeOrder(100);
% plant.set_polytopeOrder(20);
error = 0.01;
plant.options.maxError = [error; error; error];

%% Reachability analysis
% Initial set
% lb = [1.0; 0.0];
lb = [1.19; 0.19; 1];
ub = [1.2; 0.2; 1];
init_set = Star(lb,ub);
% Input set
lb = [0];
ub = [0];
input_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
reachCell = {};
reachCell{1} = init_set;
% Execute reachabilty analysis
% for i =1:steps
num_steps = 19;
t = tic;
for i=1:num_steps
    % Compute plant reachable set
    init_set = plantReach(plant,init_set, input_set);
    reachAll = [reachAll init_set];
    reachCell{i+1} = init_set;
    % Compute controller output set
    inp_nn = [];
    for i = 1:length(init_set)
        inp_nn = [inp_nn init_set(i).affineMap([1 0 0;0 1 0],[])];
    end
    input_set = net.reach(inp_nn,'approx-star');
end
timing = toc(t);
%% Set output path
path_out_sp = ['..' filesep path_results() filesep 'single_pendulum' filesep];
mkdir(path_out_sp);
save([path_out_sp 'sets'],'reachAll','timing','-v7.3');
%% Visualize results

f = figure('visible','off');
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,1,'b');
grid;
title('Single Pendulum reachable sets');
xlabel('Time steps');
ylabel('Theta');
saveas(f,[path_out_sp 'plot.jpg']);

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