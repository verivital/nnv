%% Reachability analysis of the aircarft attitude benchmark
% Logsig operations error out when trying to find the minimum value (lp error), not
% sure why. Can we do the getRanges to fix this?

net = Load_nn('model.mat');
controlPeriod = 0.1;
reachstep = 0.01;
plant = NonLinearODE(6,3,@dynamics, reachstep, controlPeriod, eye(6));
plant.set_taylorTerms(4);
plant.set_zonotopeOrder(20);
plant.set_tensorOrder(2);
% plant.set_intermediateOrder(50);
% plant.set_polytopeOrder(50);% error = 0.001;
% error = 0.0005;
% plant.options.maxError = [error; error; error; error];

%% Reachability analysis
% Initial set
% lb = [-0.45; -0.55; 0.65; -0.75; 0.85; -0.65];
% ub = [-0.44; -0.54; 0.66; -0.74; 0.86; -0.64];
lb = [-0.45; -0.55; 0.66; -0.75; 0.8595; -0.65];
ub = [-0.4495; -0.5495; 0.6595; -0.7495; 0.86; -0.6495];
% ub = [-0.449; -0.549; 0.66; -0.749; 0.86; -0.649];
% init_set = Box(lb,ub);
init_set = Star(lb,ub);
% init = init_set.partition([1 2 3 4 5 6],[50 50 50 50 100 50]);

% Store all reachable sets
% reachAll = cell(length(init),1);
% Execute reachabilty analysis
steps = 20;
t = tic;
% for j = 1:length(init)
% for j = 1:30
%     init_set = init(j).toStar;
reachAll = init_set;
for i = 1:steps
    % Compute controller output set
    input_set = net.reach(init_set,'approx-star',1);
    % Compute plant reachable set
%     init_set = plant.stepReachStar(init_set, input_set,'lin');
    init_set = plantReach(plant,init_set,input_set,'lin');
    reachAll = [reachAll init_set];
end
%     reachAll{j} = reachSub;
% end
t = toc(t);

path_out = [path_results(), filesep, 'attitude', filesep];
mkdir(path_out)
save([path_out, 'reach.mat'],'t','reachAll','-v7.3')
%% Visualize results

% Goal: avoid this region
% x1 ∈ [-0.2,0],     x2 ∈ [-0.5, -0.4],  x3 ∈ [0, 0.2]
% x4 ∈ [-0.7, -0.6]  x5 ∈ [0.7, 0.8],  x6 ∈ [-0.4, -0.2]

f = figure;
% rectangle('Position',[-0.2,-0.5,0.2,0.1],'FaceColor',[0.5 0 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
grid;
xlabel('x1');
ylabel('x2');
% saveas(f,[path_out, 'reach1v2.pdf']);

f1 = figure;
% rectangle('Position',[0,-0.7,0.2,0.1],'FaceColor',[0.5 0 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
grid;
xlabel('x3');
ylabel('x4');
% saveas(f1,[path_out, 'reach3v4.pdf']);

f1 = figure;
% rectangle('Position',[0.7,-0.4,0.1,0.2],'FaceColor',[0.5 0 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,5,6,'b');
grid;
xlabel('x5');
ylabel('x6');
% saveas(f1,[path_out, 'reach5v6.pdf']);

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