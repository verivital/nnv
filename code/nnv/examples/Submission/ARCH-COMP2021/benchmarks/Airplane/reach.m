%% Reachability analysis of Airplane Benchmark

%% load the controller
net = Load_nn('controller_airplane.mat');

% Specify the reach step, has to be smaller than the control period
reachStep = 0.001;
%% specify the control period as specified by the benchmark description
controlPeriod = 0.1;

% define the plant as specified by nnv
plant = NonLinearODE(12,6,@dynamics, reachStep, controlPeriod, eye(12));
plant.set_taylorTerms(20)
plant.set_zonotopeOrder(10);

%% Reachability analysis (lin)
% Initial set
lb = [0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0];
ub = [0.0;0.0;0.0;0.01;0.01;0.01;0.01;0.01;0.01;0.0;0.0;0.0];
init_set = Star(lb,ub);
% Input set
lb = [0;0;0;0;0;0];
ub = [0;0;0;0;0;0];
input_set = Star(lb,ub);
% Store all reachable sets
reachLin = init_set;
% Execute reachabilty analysis
num_steps = 13;
t = tic;
for i=1:num_steps
    % Compute plant reachable set
%     init_set = plant.stepReachStar(init_set(1),input_set(1));
    init_set = plantReach(plant,init_set,input_set,'lin');
    reachLin = [reachLin init_set];
    % Compute controller output set
    input_set = net.reach(init_set,'approx-star');
end
tlin = toc(t);

%% Reachability analysis (poly)
plant2 = NonLinearODE(12,6,@dynamics, reachStep, controlPeriod, eye(12));
plant2.set_taylorTerms(20)
plant2.set_zonotopeOrder(10);
% Initial set
lb = [0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0];
ub = [0.0;0.0;0.0;0.01;0.01;0.01;0.01;0.01;0.01;0.0;0.0;0.0];
init_set = Star(lb,ub);
% Input set
lb = [0;0;0;0;0;0];
ub = [0;0;0;0;0;0];
input_set = Star(lb,ub);
% Store all reachable sets
reachPoly = init_set;
% Execute reachabilty analysis
num_steps = 13;
t = tic;
for i=1:num_steps
    % Compute plant reachable set
%     init_set = plant.stepReachStar(init_set(1),input_set(1));
    init_set = plantReach(plant2,init_set,input_set,'poly');
    reachPoly = [reachPoly init_set];
    % Compute controller output set
    input_set = net.reach(init_set,'approx-star');
end
tpoly = toc(t);
path_out = [path_results(), filesep, 'Airplane', filesep];
mkdir(path_out)
save([path_out, 'reach.mat'],'reachPoly','reachLin','tlin','tpoly','-v7.3');
%% Visualize results

f2 = figure;
% Star.plotBoxes_2D_noFill(reachAll,2,5,'b');
grid;hold on;
Star.plotBoxes_2D_noFill(reachPoly,2,5,'m');
plot([-0.5 -0.5],[-2 2],'r');
plot([0.5 0.5],[-2 2],'r');
xlabel('y');
ylabel('v');
saveas(f2,[path_out, 'reach2vs5_cP.pdf']);

f4 = figure;
% Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,7,10,'b');
grid;hold on;
Star.plotBoxes_2D_noFill(reachPoly,7,10,'m');
plot([-1 -1],[-0.2 0.2],'r');
plot([1 1],[-0.2 0.2],'r');
xlabel('x_7');
ylabel('x_{10}');
saveas(f4,[path_out, 'reach7vs10.pdf']);

f5 = figure;
% Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,8,11,'b');
grid;hold on;
Star.plotBoxes_2D_noFill(reachPoly,8,11,'m');
plot([-1 -1],[-0.2 0.2],'r');
plot([1 1],[-0.2 0.2],'r');
xlabel('x_8');
ylabel('x_{11}');
saveas(f5,[path_out, 'reach8vs11.pdf']);

f6 = figure;
% Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,9,12,'b');
grid;hold on;
Star.plotBoxes_2D_noFill(reachPoly,9,12,'m');
plot([-1 -1],[-0.2 0.2],'r');
plot([1 1],[-0.2 0.2],'r');
xlabel('x_9');
ylabel('x_{12}');
saveas(f6,[path_out, 'reach9vs12.pdf']);

%%% Visualize Poly results
f2 = figure;
grid;hold on;
Star.plotBoxes_2D_noFill(reachPoly,2,5,'b');
plot([-0.5 -0.5],[-2 2],'r');
plot([0.5 0.5],[-2 2],'r');
xlabel('y');
ylabel('v');
saveas(f2,[path_out, 'reach2vs5_cP_poly.pdf']);

f4 = figure;
grid;hold on;
Star.plotBoxes_2D_noFill(reachPoly,7,10,'b');
plot([-1 -1],[-0.2 0.2],'r');
plot([1 1],[-0.2 0.2],'r');
xlabel('x_7');
ylabel('x_{10}');
saveas(f4,[path_out, 'reach7vs10_poly.pdf']);

f5 = figure;
grid;hold on;
Star.plotBoxes_2D_noFill(reachPoly,8,11,'b');
plot([-1 -1],[-0.2 0.2],'r');
plot([1 1],[-0.2 0.2],'r');
xlabel('x_8');
ylabel('x_{11}');
saveas(f5,[path_out, 'reach8vs11_poly.pdf']);

f6 = figure;
% Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,9,12,'b');
grid;hold on;
Star.plotBoxes_2D_noFill(reachPoly,9,12,'b');
plot([-1 -1],[-0.2 0.2],'r');
plot([1 1],[-0.2 0.2],'r');
title('Airplane x_9 vs. x_{12}');
xlabel('x_9');
ylabel('x_{12}');
saveas(f6,[path_out, 'reach9vs12_poly.pdf']);
toc(t);
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