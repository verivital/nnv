%% Reachability analysis of Double Pendulum Benchmark

%% load the controller
net = Load_nn('controller_double_pendulum_more_robust.mat');

% Specify the reach step, has to be smaller than the control period
reachStep = 0.001;
%% specify the control period as specified by the benchmark description
controlPeriod = 0.02;

% define the plant as specified by nnv
plant = NonLinearODE(4,2,@dynamics_dp, reachStep, controlPeriod, eye(4));
% plant.set_zonotopeOrder(50);
% plant.set_polytopeOrder(20);
% error = 0.01;
% plant.options.maxError = [error; error;error;error];

%% Reachability analysis
% Initial set
lb = [1.0; 1.0;1.0;1.0];
% ub = [1.3; 1.3;1.3;1.3];
ub = [1.01; 1.01;1.01;1.01];
init_set = Star(lb,ub);
% Input set
lb = [0;0];
ub = [0;0];
input_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
% Execute reachabilty analysis
% for i =1:steps
num_steps = 5;
t = tic;
for i=1:num_steps
    % Compute plant reachable set
    init_set = plantReach(plant,init_set, input_set);
    reachAll = [reachAll init_set];
    % Compute controller output set
    input_set = net.reach(init_set,'approx-star');
end
timing = toc(t);
%% Set output path
path_out_dp = ['..' filesep path_results() filesep 'double_pendulum' filesep];
mkdir(path_out_dp);
save([path_out_dp 'sets_more'],'reachAll','timing','-v7.3');

%% Visualize results
f = figure('visible','off');
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
title('Double Pendulum reachable sets');
xlabel('x1');
ylabel('x2');
saveas(f,[path_out_dp 'plot_more_1v2.jpg']);

f1 = figure('visible','off');
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
grid;
title('Double Pendulum reachable sets');
xlabel('x3');
ylabel('x4');
saveas(f1,[path_out_dp 'plot_more_3v4.jpg']);

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