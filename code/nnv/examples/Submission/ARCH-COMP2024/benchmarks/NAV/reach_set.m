function t = reach_set()

    %% Reachability analysis of NAV Benchmark
    
    %% Load Components 

    % Load the controller
    net = importNetworkFromONNX('networks/nn-nav-set.onnx', "InputDataFormats", "BC");
    net = matlab2nnv(net);
    % Load plant
    reachStep = 0.01;
    controlPeriod = 0.02;
    plant = NonLinearODE(4, 2, @dynamics, reachStep, controlPeriod, eye(4));
    plant.set_tensorOrder(2);


    %% Reachability analysis

    % Initial set
    lb = [2.9; 2.9; 0; 0];
    ub = [3.1; 3.1; 0; 0];
    init_set = Star(lb,ub);

    % Store all reachable sets
    reachAll = init_set;

    % Reachability options
    num_steps = 30;
    reachOptions.reachMethod = 'approx-star';

    % Begin computation
    t = tic;
    for i=1:num_steps

        % Compute controller output set
        input_set = net.reach(init_set,reachOptions);

        % Compute plant reachable set
        init_set = plant.stepReachStar(init_set, input_set,'lin');
        reachAll = [reachAll init_set];

    end

    t = toc(t);


    %% Visualize results
    plant.get_interval_sets;

    f1 = figure;
    % rectangle('Position',[-1,-1,2,2],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    hold on;
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
    grid;
    xlabel('x_1');
    ylabel('x_2');
    
    f2 = figure;
    % rectangle('Position',[-1,-1,2,2],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    hold on;
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
    grid;
    xlabel('x_3');
    ylabel('x_4');

    % Save figure
    % if is_codeocean
    % 
    % else
    % 
    % end

end