function t = reach()

%% Reachability analysis of TORA (benchmark 9)

    % Load components and set reachability parameters
    net = load_NN_from_mat('controllerTora.mat');
    reachStep = 0.1;
    controlPeriod = 1;
    plant = NonLinearODE(4,1,@dynamics9, reachStep, controlPeriod, eye(4));
    tFinal = 20;
    % Initial set
    lb = [0.6; -0.7; -0.4; 0.5];
    ub = [0.7; -0.6; -0.3; 0.6];
%     offset = 10; % Applied this in the dynamics9.m function
%     scale_factor = 1;
    
    % Custom plant reachability options
    plant.options.taylorTerms = 4;
    plant.options.zonotopeOrder = 200;
    alg = 'poly';
    plant.options.alg = alg;
    plant.options.tensorOrder = 3;
    plant.options.errorOrder = 10;
    plant.options.intermediateOrder = 50;
    
    %% Reachability analysis
    init_set = Star(lb,ub);
    % Store all reachable sets
    reachAll = init_set;
    % Execute reachabilty analysis
    reachOptions.reachMethod = 'approx-star';
    t = tic;
    for i = 1:tFinal
        % Compute controller output set
        input_set = net.reach(init_set,reachOptions);
        % Compute plant reachable set
        init_set = plant.stepReachStar(init_set, input_set,alg);
        reachAll = [reachAll init_set];
    end
    t = toc(t);
    
    % Save results
    if is_codeocean
        save('/results/logs/tora.mat', 'reachAll','t','-v7.3');
    else
        save('tora.mat', 'reachAll','t','-v7.3');
    end

%% Visualize results
    plant.get_interval_sets;
    
    f = figure;
    grid;
    hold on;
    rectangle('Position',[-2,-2,4,4],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
    grid;
    xlabel('x1');
    ylabel('x2');
    
    f1 = figure;
    grid;
    hold on;
    rectangle('Position',[-2,-2,4,4],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
    grid;
    xlabel('x3');
    ylabel('x4');

    % Save figures
    if is_codeocean
        exportgraphics(f,'/results/logs/tora_1v2.pdf', 'ContentType', 'vector');
        exportgraphics(f1,'/results/logs/tora_3v4.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'tora_1v2.pdf','ContentType', 'vector');
        exportgraphics(f,'tora_3v4.pdf','ContentType', 'vector');
    end

end