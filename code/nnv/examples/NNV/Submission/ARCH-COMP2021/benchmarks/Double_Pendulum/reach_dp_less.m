%% Reachability analysis of Double Pendulum Benchmark

%% load the controller
net = Load_nn('controller_double_pendulum_less_robust.mat');

% Specify the reach step, has to be smaller than the control period
reachStep = 0.001;
%% specify the control period as specified by the benchmark description
controlPeriod = 0.05;

% define the plant as specified by nnv
plant = NonLinearODE(4,2,@dynamics_dp, reachStep, controlPeriod, eye(4));
% plant.set_zonotopeOrder(50);


%% Reachability analysis
% Initial set
% lb = [1.0; 1.0;1.0;1.0];
lb = [1.29; 1.29;1.29;1.29];
ub = [1.3; 1.3;1.3;1.3];
% ub = [1.01; 1.01;1.01;1.01];
init_set = Star(lb,ub);
% Input set
lb = [0;0];
ub = [0;0];
input_set = Star(lb,ub);
% Store all reachable sets
reachLin = init_set;
% Execute reachabilty analysis
% for i =1:steps
num_steps = 4;
t = tic;
for i=1:num_steps
    % Compute plant reachable set
    init_set = plantReach(plant,init_set, input_set,'lin');
    reachLin = [reachLin init_set];
    % Compute controller output set
    input_set = net.reach(init_set,'approx-star');
end
tlin = toc(t);

% Reach 2
plant2 = NonLinearODE(4,2,@dynamics_dp, reachStep, controlPeriod, eye(4));
% Initial set
lb = [1.29; 1.29;1.29;1.29];
ub = [1.3; 1.3;1.3;1.3];
init_set = Star(lb,ub);
% Input set
lb = [0;0];
ub = [0;0];
input_set = Star(lb,ub);
% Store all reachable sets
reachPoly = init_set;
% Execute reachabilty analysis
num_steps = 4;
t = tic;
for i=1:num_steps
    % Compute plant reachable set
    init_set = plantReach(plant2,init_set, input_set,'poly');
    reachPoly = [reachPoly init_set];
    % Compute controller output set
    input_set = net.reach(init_set,'approx-star');
end
tpoly = toc(t);
path_out = [path_results(), filesep, 'DoublePendulum', filesep];
mkdir(path_out)
save([path_out, 'reach_less.mat'],'reachLin','reachPoly','tlin','tpoly','-v7.3');

%% Visualize results
f = figure;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach_less_1v2.pdf']);

f1 = figure;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
grid;
xlabel('x3');
ylabel('x4');
saveas(f1,[path_out, 'reach_less_3v4.pdf']);

f = figure;
Star.plotBoxes_2D_noFill(plant2.intermediate_reachSet,1,2,'b');
grid;
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach_less_1v2_poly.pdf']);

f1 = figure;
Star.plotBoxes_2D_noFill(plant2.intermediate_reachSet,3,4,'b');
grid;
xlabel('x3');
ylabel('x4');
saveas(f1,[path_out, 'reach_less_3v4_poly.pdf']);

%% Helper function
function init_set = plantReach(plant,init_set,input_set,algoC)
    nS = length(init_set);
    nL = length(input_set);
    ss = [];
    for k=1:nS
        for l=1:nL
            ss =[ss plant.stepReachStar(init_set(k), input_set(l),algoC)];
        end
    end
    init_set = ss;
end