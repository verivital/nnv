%% Reachability analysis of Single Pendulum Benchmark
    
%% Load Components 

% Load the controller
net = load_NN_from_mat('controller_single_pendulum.mat');
% Load plant
reachStep = 0.01;
controlPeriod = 0.05;
plant = NonLinearODE(2,1,@dynamics_sp, reachStep, controlPeriod, eye(2));
% plant.set_tensorOrder(2);


%% Reachability analysis

% Initial set
lb =  [0; -0.1];
ub = [1; 0.1];
% lb = [1.1999; 0.1999; 0];
% ub = [1.2; 0.2; 0];
init_set = Box(lb,ub);
init_partitions = init_set.partition([1 2], [20, 8]);
% Store all reachable sets
outputPartitions = cell(length(init_partitions),1);
num_steps = 20;
reachOptions.reachMethod = 'approx-star';
for s = 1:length(init_partitions)
    % Get initial set for corresponding partition
    init_set = init_partitions(s);
    init_set = Star(init_set.lb, init_set.ub);
    reachAll = init_set;
    t = tic;
    % Begin NNCS reachability for partition
    for i = 1:num_steps
        % Compute controller output set
        input_set = net.reach(init_set,reachOptions);
        input_set = Star.get_hypercube_hull(input_set);
        input_set = input_set.toStar;
        % Compute plant reachable set
        init_set = plant.stepReachStar(init_set, input_set,'lin');
        init_set = Star.get_hypercube_hull(init_set);
        init_set = init_set.toStar;
        reachAll = [reachAll init_set];
    end
    toc(t);
    % store results for partition
    outputPartitions{s} = reachAll;
end

%% Visualize results
timeVector = 0:controlPeriod:controlPeriod*num_steps;

f = figure;
hold on;
rectangle('Position',[0.5,1,1,1],'FaceColor',[1 0 0 0.5],'EdgeColor','r', 'LineWidth',0.1)
for i=1:length(outputPartitions)
    Star.plotRanges_2D(outputPartitions{i},1,timeVector,'b');
    hold on;
end
grid;
xlabel('Time (s)');
ylabel('\theta');
% xlim([0 0.6])
% ylim([0.95 1.25])
% Save figure
if is_codeocean
    exportgraphics(f,'/results/logs/singlePendulum.pdf', 'ContentType', 'vector');
else
    exportgraphics(f,'singlePendulum.pdf','ContentType', 'vector');
end
