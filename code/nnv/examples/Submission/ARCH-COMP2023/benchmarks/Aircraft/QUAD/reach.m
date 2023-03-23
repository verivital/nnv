%% Reachability analysis of the aircarft quad benchmark
% Logsig operations error out when trying to find the minimum value (lp error), not
% sure why. Can we do the getRanges to fix this?

net = Load_nn('model.mat');
controlPeriod = 0.1;
reachstep = 0.001;
plant = NonLinearODE(12,3,@dynamics, reachstep, controlPeriod, eye(12));
plant.set_taylorTerms(8);
plant.set_zonotopeOrder(50);
plant.set_tensorOrder(3);
% plant.set_intermediateOrder(50);
% plant.set_polytopeOrder(50);% error = 0.001;
% error = 0.0005;
% plant.options.maxError = [error; error; error; error];

%% Reachability analysis
% Initial set
% lb = [-0.4; -0.4 ; -0.4; -0.4; -0.4; -0.4; 0; 0; 0; 0; 0; 0];
% ub = [0.4; 0.4 ; 0.4; 0.4; 0.4; 0.4; 0; 0; 0; 0; 0; 0];
lb = [-0.01; -0.01; -0.01; -0.01; -0.01; -0.01; 0; 0; 0; 0; 0; 0];
ub = [0.01; 0.01 ; 0.01; 0.01; 0.01; 0.01; 0; 0; 0; 0; 0; 0];
init_set = Star(lb,ub);
% Input set
% lb = [0;0];
% ub = [0;0];
% input_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
% Execute reachabilty analysis
steps = 10;
t = tic;
prop = true;
for i = 1:steps
    % Compute controller output set
    input_set = net.reach(init_set,'approx-star');
    % Compute plant reachable set
    init_set = plantReach(plant,init_set,input_set,'lin');
    reachAll = [reachAll init_set];
end
t = toc(t);

path_out = [path_results(), filesep, 'quad', filesep];
mkdir(path_out)
save([path_out, 'reach.mat'],'t','reachAll','-v7.3')
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