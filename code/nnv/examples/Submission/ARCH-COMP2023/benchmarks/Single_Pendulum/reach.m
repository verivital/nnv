function t = reach()

    %% Reachability analysis of Single Pendulum Benchmark
    
    %% Load Components 

    % Load the controller
    net = load_NN_from_mat('controller_single_pendulum.mat');
    % Load plant
    reachStep = 0.005;
    controlPeriod = 0.05;
    plant = NonLinearODE(3,1,@dynamics_sp, reachStep, controlPeriod, eye(3));
    plant.set_tensorOrder(2);
    

    %% Reachability analysis

    % Initial set (smaller than specification to prove not safe)
    % lb = [1.0; 0.0; 0];
    lb = [1.1999; 0.1999; 0];
    ub = [1.2; 0.2; 0];
    init_set = Star(lb,ub);
    % Store all reachable sets
    reachAll = init_set;
    num_steps = 11;
    reachOptions.reachMethod = 'approx-star';
    t = tic;
    for i=1:num_steps
        % Compute controller output set
        init_set_s = init_set.affineMap([1 0 0;0 1 0],[]);
        input_set = net.reach(init_set_s,reachOptions);
        % Compute plant reachable set
        init_set = plant.stepReachStar(init_set, input_set,'lin');
        reachAll = [reachAll init_set];
    end
    t = toc(t);
    
    %% Visualize results
    plant.get_interval_sets;
    
    f = figure;
    hold on;
    rectangle('Position',[0.5,1,1,1],'FaceColor',[1 0 0 0.5],'EdgeColor','r', 'LineWidth',0.1)
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,1,'b');
    % Plot only falsifying trace
%     plot(squeeze(sims(3,k,:)), squeeze(sims(1,k,:)), 'Color', [0 0 1 1]);
    grid;
    xlabel('Time (s)');
    ylabel('\theta');
    xlim([0 0.6])
    ylim([0.95 1.25])
    % Save figure
    if is_codeocean
        exportgraphics(f,'/results/logs/singlePendulum.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'singlePendulum.pdf','ContentType', 'vector');
    end

end