% function t = reach()

    %% Reachability analysis of Cartpole Benchmark
    
    %% Load Components 

    % Load the controller
    net = importNetworkFromONNX('model.onnx', "InputDataFormats", "BC");
    net = matlab2nnv(net);
    % Load plant
    reachStep = 0.002;
    controlPeriod = 0.02;
    plant = NonLinearODE(4,1,@dynamics, reachStep, controlPeriod, eye(4));
    plant.set_tensorOrder(2);


    %% Reachability analysis

    % Initial set
    lb = [-0.1; -0.05; -0.1; -0.05];
    ub = [0.1; 0.05; 0.1; 0.05];
    InitialSet = Box(lb,ub);
    init_sets = InitialSet.partition([1,2,3,4],[20,10,20,10]);
    for k=1:length(init_sets)
        init_set = init_sets(k);
        init_set = init_set.toStar;
        % Store all reachable sets
        reachAll = init_set;
        num_steps = 500;
        reachOptions.reachMethod = 'approx-star';
        t = tic;
        for i=1:num_steps
            disp(i);
            % Compute controller output set
            input_set = net.reach(init_set,reachOptions);
            % Compute plant reachable set
            init_set = plant.stepReachStar(init_set, input_set,'lin');
            reachAll = [reachAll init_set];
            toc(t);
        end
        t = toc(t);
    end
    
    %% Visualize results
    plant.get_interval_sets;
    
    f = figure;
    hold on;
    % rectangle('Position',[0.5,1,1,1],'FaceColor',[1 0 0 0.5],'EdgeColor','r', 'LineWidth',0.1)
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'r');
    % Plot only falsifying trace
%     plot(squeeze(sims(3,k,:)), squeeze(sims(1,k,:)), 'Color', [0 0 1 1]);
    grid;
    % xlabel('Time (s)');
    % ylabel('\theta');
    % xlim([0 0.6])
    % ylim([0.95 1.25])
    % Save figure
    % if is_codeocean
    %     exportgraphics(f,'/results/logs/cartpole.pdf', 'ContentType', 'vector');
    % else
    %     exportgraphics(f,'cartpole.pdf','ContentType', 'vector');
    % end

% end