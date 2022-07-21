%% Reachability analysis of TORA (benchmark 9)
% Load components and set reachability parameters
% net = Load_nn('controllerTora.onnx');
net = Load_nn('controllerTora_nnv.mat');
reachStep = 0.005;
controlPeriod = 1;
plant = NonLinearODE(4,1,@dynamics9, reachStep, controlPeriod, eye(4));
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
reachLin = init_set;
% Execute reachabilty analysis
% for i =1:steps
t = tic;
for i =1:6
    % Compute plant reachable set
    init_set = plantReach(plant,init_set,input_set,'lin');
    reachLin = [reachLin init_set];
    % Compute controller output set
    inp_set = net.reach(init_set,'approx-star');
    input_set = [];
    for k = 1:length(inp_set)
        input_set = [input_set inp_set(k).affineMap(1,-offset)];
    end
end
tlin = toc(t);

% Reach #2
plant2 = NonLinearODE(4,1,@dynamics9, reachStep, controlPeriod, eye(4));
plant2.options.taylorTerms = 10;
plant2.options.zonotopeOrder = 10;
% plant2.options.intermediateOrder = 10;
plant2.options.errorOrder = 5;

lb = [0.69; -0.7; -0.4; 0.5];
ub = [0.7; -0.69; -0.39; 0.51];
init_set = Star(lb,ub);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
reachPoly = init_set;
% Execute reachabilty analysis
% for i =1:steps
t = tic;
for i =1:6
    % Compute plant reachable set
    init_set = plantReach(plant2,init_set,input_set,'poly');
    reachPoly = [reachPoly init_set];
    % Compute controller output set
    inp_set = net.reach(init_set,'approx-star');
    input_set = [];
    for k = 1:length(inp_set)
        input_set = [input_set inp_set(k).affineMap(1,-offset)];
    end
end
tpoly = toc(t);

path_out = [path_results(), filesep, 'benchmark9-tora', filesep];
mkdir(path_out)
save([path_out, 'reach.mat'],'reachPoly','reachLin','tlin','tpoly','-v7.3')

%% Visualize results
f = figure;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachLin,1,2,'m');
grid;
xlabel('x1');
ylabel('x3');
saveas(f,[path_out, 'reach1vs3.pdf']);

f1 = figure;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,2,4,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachLin,2,4,'m');
grid;
xlabel('x2');
ylabel('x4');
saveas(f1,[path_out, 'reach3vs4.pdf']);

% Poly results
f = figure;
Star.plotBoxes_2D_noFill(plant2.intermediate_reachSet,1,2,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachLin,1,2,'m');
grid;
xlabel('x1');
ylabel('x3');
saveas(f,[path_out, 'reach1vs3_poly.pdf']);

f1 = figure;
Star.plotBoxes_2D_noFill(plant2.intermediate_reachSet,2,4,'b');
grid;
hold on;
Star.plotBoxes_2D_noFill(reachLin,2,4,'m');
grid;
xlabel('x2');
ylabel('x4');
saveas(f1,[path_out, 'reach3vs4_poly.pdf']);

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