%% Reachability analysis of TORA (benchmark 9)
% Load components and set reachability parameters
% net = Load_nn('controllerTora.onnx');
net = Load_nn('controllerTora_nnv.mat');
reachStep = 0.05;
controlPeriod = 1;
plant = NonLinearODE(4,1,@dynamics9, reachStep, controlPeriod, eye(4));
time = 0:controlPeriod:20;
steps = length(time);
% Initial set
lb = [0.6; -0.7; -0.4; 0.5];
% ub = [0.7; -0.6; -0.3; 0.6];
ub = [0.6; -0.7; -0.4; 0.5];
% ub = [0.7; -0.69; -0.39; 0.51];
% init_set = Box(lb,ub);
offset = 10;
scale_factor = 1;
% plant.set_tensorOrder(3);
plant.options.intermediateOrder = 20;
plant.options.zonotopeOrder = 20;
plant.options.taylorTerms  = 5;
% plant.set_maxError(0.001*ones(4,1))

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
for i=1:15
    % Compute controller output set
    inp_set = net.reach(init_set,'approx-star');
    input_set = [];
    for k = 1:length(inp_set)
        input_set = [input_set inp_set(k).affineMap(1,-offset)];
    end
    % Compute plant reachable set
    init_set = plantReach(plant,init_set,input_set,'lin');
    reachLin = [reachLin init_set];
end
tlin = toc(t);

path_out = [path_results(), filesep, 'benchmark9-tora', filesep];
mkdir(path_out)
save([path_out, 'reach.mat'],'reachLin','tlin','-v7.3')

%% Visualize results
f = figure;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
hold on;
% Star.plotBoxes_2D_noFill(reachLin,1,2,'m');
grid;
xlabel('x1');
ylabel('x2');
saveas(f,[path_out, 'reach1vs3.pdf']);

f1 = figure;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
grid;
hold on;
% Star.plotBoxes_2D_noFill(reachLin,3,4,'m');
grid;
xlabel('x3');
ylabel('x4');
saveas(f1,[path_out, 'reach3vs4.pdf']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Helper functions
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

%% Notes

% Even when we start with no uncertainty, it just builds up very quickly,
% unable to verify the system.

% Reduce the number of star sets by checking if one is contained in the others
% We'll do by using the ranges function
% function stars = reduce_stars(input_set)
% 
% 
% end