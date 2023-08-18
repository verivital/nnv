% function t = reach()

%% Reachability analysis of the aircarft attitude benchmark

%% Load Components

    net = load_NN_from_mat('model.mat');
    controlPeriod = 0.1;
    reachstep = 0.01;
    plant = NonLinearODE(6,3,@dynamics, reachstep, controlPeriod, eye(6));
    plant.set_taylorTerms(4);
    plant.set_zonotopeOrder(50);
    plant.set_tensorOrder(2);

%% Reachability analysis
    % Initial set
    lb = [-0.45; -0.55; 0.65; -0.75; 0.85; -0.65];
%     ub = [-0.44; -0.54; 0.66; -0.74; 0.86; -0.64];
    % lb = [-0.45; -0.55; 0.66; -0.75; 0.8595; -0.65];
    ub = [-0.4495; -0.5495; 0.6595; -0.7495; 0.8595; -0.6495];
    % ub = [-0.449; -0.549; 0.66; -0.749; 0.86; -0.649];
%     init_set = Box(lb,ub);
    init_set = Star(lb,ub);
%     init = init_set.partition([1 2 3 4 5 6],[4 4 4 4 4 4]);
    
    % Store all reachable sets
    % reachAll = cell(length(init),1);
    % Execute reachabilty analysis
    steps = 16;
    
    % for j = 1:length(init)
    % for j = 1:30
    %     init_set = init(j).toStar;
    reachAll = init_set;
    reachOptions.reachMethod = 'approx-star';
    t = tic;
    for i = 1:steps
        % Compute controller output set
        input_set = net.reach(init_set, reachOptions);
        % Compute plant reachable set
    %     init_set = plant.stepReachStar(init_set, input_set,'lin');
        init_set = plantReach(plant,init_set,input_set,'lin');
        reachAll = [reachAll init_set];
    end
    %     reachAll{j} = reachSub;
    % end
    t = toc(t);

    % Save results
    if is_codeocean
        save('/results/logs/attitude.mat', 'reachAll','t','-v7.3');
    else
        save('attitude.mat', 'reachAll','t','-v7.3');
    end


    %% Visualize results

    plant.get_interval_sets;
    
    % Goal: avoid this region
    % x1 ∈ [-0.2,0],     x2 ∈ [-0.5, -0.4],  x3 ∈ [0, 0.2]
    % x4 ∈ [-0.7, -0.6]  x5 ∈ [0.7, 0.8],  x6 ∈ [-0.4, -0.2]
    
    f = figure;
    rectangle('Position',[-0.2,-0.5,0.2,0.1],'FaceColor',[0.5 0 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    hold on;
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
    grid;
    xlabel('x1');
    ylabel('x2');
    
    f1 = figure;
    rectangle('Position',[0,-0.7,0.2,0.1],'FaceColor',[0.5 0 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    hold on;
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
    grid;
    xlabel('x3');
    ylabel('x4');
    
    f2 = figure;
    rectangle('Position',[0.7,-0.4,0.1,0.2],'FaceColor',[0.5 0 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    hold on;
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,5,6,'b');
    grid;
    xlabel('x5');
    ylabel('x6');

    % Save figure
    if is_codeocean
        exportgraphics(f, '/results/logs/attitude_1v2.pdf', 'ContentType', 'vector');
        exportgraphics(f1,'/results/logs/attitude_3v4.pdf', 'ContentType', 'vector');
        exportgraphics(f2,'/results/logs/attitude_5v6.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'attitude_1v2.pdf','ContentType', 'vector');
        exportgraphics(f,'attitude_3v4.pdf','ContentType', 'vector');
        exportgraphics(f,'attitude_5v6.pdf','ContentType', 'vector');
    end

% end

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