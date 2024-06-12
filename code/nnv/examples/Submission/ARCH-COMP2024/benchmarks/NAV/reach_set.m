function rT = reach_set()

    %% Reachability analysis of NAV Benchmark
    
    %% Load Components 

    % Load the controller
    net = importNetworkFromONNX('networks/nn-nav-set.onnx', "InputDataFormats", "BC");
    net = matlab2nnv(net);
    % Load plant
    reachStep = 0.02;
    controlPeriod = 0.2;
    plant = NonLinearODE(4, 2, @dynamics, reachStep, controlPeriod, eye(4));
    plant.options.zonotopeOrder = 160;
    plant.options.intermediateOrder = 200;
    plant.options.taylorTerms = 3;
    plant.set_tensorOrder(3);


    %% Reachability analysis

    % Initial set
    lb = [2.999; 2.999; 0; 0];
    ub = [3.0; 3.0; 0; 0];
    init_set = Star(lb,ub);

    % Store all reachable sets
    reachAll = init_set;

    % Reachability options
    num_steps = 28;
    reachOptions.reachMethod = 'approx-star';

    % Begin computation
    t = tic;
    for i=1:num_steps

        % Compute controller output set
        input_set = net.reach(init_set,reachOptions);
        
        % input_set = input_set.getBox;
        % input_set = Star.get_hypercube_hull(input_set);
        % input_set = input_set.toStar;

        % Compute plant reachable set
        init_set = plant.stepReachStar(init_set, input_set,'lin-adaptive');
        % init_set = Star.get_hypercube_hull(init_set);
        % init_set = init_set.toStar;

        reachAll = [reachAll init_set];

    end

    rT = toc(t);
    
    R = reachAll;

    % Save results
    if is_codeocean
        save('/results/logs/nav_set.mat', 'R','rT','-v7.3');
    else
        save('nav_set.mat', 'R','rT','-v7.3');
    end


    %% Visualize results
    plant.get_point_sets;

    f = figure;
    rectangle('Position',[-0.5,-0.5,1,1],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth', 0.1); % goal region
    hold on;
    rectangle('Position',[1,1,1,1],'FaceColor',[0.7 0 0 0.8], 'EdgeColor','r', 'LineWidth', 0.1); % obstacle
    Star.plotBoxes_2D_noFill(plant.intermediate_pointSet,1,2,'b');
    grid;
    xlabel('x_1');
    ylabel('x_2');
    
    % f5 = figure;
    % % rectangle('Position',[-1,-1,2,2],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    % hold on;
    % Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
    % grid;
    % xlabel('x_3');
    % ylabel('x_4');

    % Save figure
    if is_codeocean
        exportgraphics(f,'/results/logs/nav-set.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'nav-set.pdf','ContentType', 'vector');
    end

    % end

end