% function t = reach()

%% Reachability analysis of the Docking  benchmark

%% Load components

net = load_NN_from_mat('model.mat');
controlPeriod = 1;
reachstep = 0.05;
plant = NonLinearODE(4,2,@dynamics, reachstep, controlPeriod, eye(4));
plant.set_taylorTerms(4);
plant.set_zonotopeOrder(50);
plant.set_tensorOrder(2);
% plant.set_intermediateOrder(50);
% plant.set_polytopeOrder(50);% error = 0.001;
% error = 0.0005;
% plant.options.maxError = [error; error; error; error];

%% Reachability analysis

% Initial set
lb = [87;  87 ; -0.01; -0.01];
% ub = [106; 106;  0.28;  0.28];
ub = [89; 89;  0.01;  0.01];
init_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
reachOptions.reachMethod = 'approx-star';
% Execute reachabilty analysis
steps = 4;
i = 1;
prop = true;
t = tic;
while prop && i < steps
    % Compute controller output set
    input_set = net.reach(init_set,reachOptions);
    % Compute plant reachable set
    init_set = plantReach(plant,init_set,input_set,'lin');
    reachAll = [reachAll init_set];
    prop = verify_spec(init_set);
    i = i+1;
end
t = toc(t);


%% Visualize results
f = figure;
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
xlabel('x1');
ylabel('x2');
% saveas(f,[path_out, 'reach1v2.pdf']);

f1 = figure;
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
grid;
xlabel('x1');
ylabel('x2');
% saveas(f1,[path_out, 'reach3v4.pdf']);

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