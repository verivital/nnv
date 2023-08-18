function t = reachTora_reluTanh()

    %% Reachability analysis of TORA (benchmark 9)
    
    % Load components and set reachability parameters
    net = load_NN_from_mat('nn_tora_relu_tanh.mat');
    reachStep = 0.05;
    controlPeriod = 0.5;
    plant = NonLinearODE(4,1,@dynamicsTora, reachStep, controlPeriod, eye(4));
    time = 0:controlPeriod:20;
    steps = length(time);
    
    % Initial set
    lb = [-0.77; -0.45; 0.51; -0.3];
    ub = [-0.75; -0.43; 0.54; -0.28];
    goal = Box([-0.1;-0.9],[0.2;-0.6]);
    offset = 0;
    scale_factor = 11;
    
    %% Reachability analysis
    
    % Initial state (partitioned)
    init_set = Box(lb,ub);
    init = init_set.partition([1 2 3 4],[4 8 6 4]);
    
    % Input set
    lb = 0;
    ub = 0;
    input_set = Star(lb,ub);
    
    % Store all reachable sets
    reachAll = cell(length(init),1);
    
    % Execute reachabilty analysis
    reachOpt.reachMethod = 'approx-star'; % controller reach options
    t = tic;
    for j = 1:length(init)
        init_set = init(j).toStar;
        reachSub = init_set;
        for i = 1:10
            % Compute controller output set
            input_set = net.reach(init_set,reachOpt);
            input_set = input_set.affineMap(scale_factor,-offset);
            % Compute plant reachable set
            init_set = plant.stepReachStar(init_set, input_set,'lin');
            reachSub = [reachSub init_set];
        end
        reachAll{j} = reachSub;
    end
    t = toc(t); % get reach time
    
    % Save results
    if is_codeocean
        save('/results/logs/tora_relu_tanh.mat', 'reachAll','t','-v7.3');
    else
        save('tora_relu_tanh.mat', 'reachAll','t','-v7.3');
    end


    %% Visualize results

    f = figure;
    rectangle('Position',[-0.1,-0.9,0.3,0.3],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    hold on;
    grid;
    for r=1:length(init)
        Star.plotBoxes_2D_noFill(reachAll{r},1,2,'b');
    end
    xlabel('x1');
    ylabel('x2');

    
    % Last control period reach sets
    f2 = figure;
    rectangle('Position',[-0.1,-0.9,0.3,0.3],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    hold on;
    for r=1:length(reachAll)
        Star.plotBoxes_2D_noFill(reachAll{r}(end),1,2,'b');
    end
    xlabel('x1');
    ylabel('x2');
    
    % Save figure
    if is_codeocean
        exportgraphics(f,'/results/logs/tora_relu_tanh.pdf', 'ContentType', 'vector');
        exportgraphics(f2,'/results/logs/tora_relu_tanh_last.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'tora_relu_tanh.pdf', 'ContentType', 'vector');
        exportgraphics(f2,'tora_relu_tanh_last.pdf','ContentType', 'vector');
    end

end