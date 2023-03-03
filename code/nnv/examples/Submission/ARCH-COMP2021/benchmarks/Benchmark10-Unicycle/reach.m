%% Reachability analysis of the Unicycle (benchmark 10)
% net = Load_nn('controllerB.onnx');
net = Load_nn('controllerB_nnv.mat');
controlPeriod = 0.2;
% controlPeriod = 0.5;
reachstep = 0.01;
plant = NonLinearODE(4,2,@dynamics10, reachstep, controlPeriod, eye(4));
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
reachLin = init_set;
% Execute reachabilty analysis
steps = 10;
nI = 1;
t = tic;
% Disclaimer: Only choose a subset of all possible input and state sets to
% show that results are inconclusive (faster computation)
for i=1:steps
    % Compute plant reachable set
    init_set = plantReach(plant,init_set(nI),input_set(1),'lin');
    nI = length(init_set);
    reachLin = [reachLin init_set(nI)];
    % Compute controller output set
    inp_set = net.reach(init_set(nI),'approx-star');
    input_set = [];
    input_set = inp_set(1).affineMap(eye(2),-offsetM);
end
tlin = toc(t);

% Reach 2
plant2 = NonLinearODE(4,2,@dynamics10, reachstep, controlPeriod, eye(4));
% Initial set
lb = [9.5; -4.5; 2.1; 1.5];
ub = [9.51; -4.49; 2.11; 1.51];
init_set = Star(lb,ub);
% Input set
lb = [0;0];
ub = [0;0];
input_set = Star(lb,ub);
% Store all reachable sets
reachPoly = init_set;
% Execute reachabilty analysis
steps = 10;
nI = 1;
t = tic;
% Disclaimer: Only choose a subset of all possible input and state sets to
% show that results are inconclusive (faster computation)
for i=1:steps
    % Compute plant reachable set
    init_set = plantReach(plant2,init_set(nI),input_set(1),'poly');
    nI = length(init_set);
    reachPoly = [reachPoly init_set(nI)];
    % Compute controller output set
    inp_set = net.reach(init_set(nI),'approx-star');
    input_set = [];
    input_set = inp_set(1).affineMap(eye(2),-offsetM);
end
tpoly = toc(t);
path_out = [path_results(), filesep, 'Unicycle', filesep];
mkdir(path_out)
save([path_out, 'reach.mat'],'tlin','tpoly','reachPoly','reachLin','-v7.3')
%% Visualize results
f = figure;
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach1v2.pdf']);

f1 = figure;
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
grid;
xlabel('x1');
ylabel('x2');
saveas(f1,[path_out, 'reach3v4.pdf']);

f = figure;
hold on;
Star.plotBoxes_2D_noFill(plant2.intermediate_reachSet,1,2,'b');
grid;
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach1v2_poly.pdf']);

f1 = figure;
hold on;
Star.plotBoxes_2D_noFill(plant2.intermediate_reachSet,3,4,'b');
grid;
xlabel('x1');
ylabel('x2');
saveas(f1,[path_out, 'reach3v4_poly.pdf']);

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