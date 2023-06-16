function t = reach_more()

    %% Reachability analysis of Double Pendulum Benchmark
    
    %% Load components

    % load controller
    net = load_NN_from_mat('controller_double_pendulum_more_robust.mat');
    % Specify the reach step, has to be smaller than the control period
    reachStep = 0.005;
    % specify the control period as specified by the benchmark description
    controlPeriod = 0.02;
    % define the plant as specified by nnv
    plant = NonLinearODE(4,2,@dynamics_dp, reachStep, controlPeriod, eye(4));
    % plant.set_tensorOrder(3);
    plant.options.intermediateOrder = 20;
    plant.options.zonotopeOrder = 20;
    plant.options.taylorTerms  = 8;
    % plant.set_maxError(0.01*ones(4,1))
    
    %% Reachability analysis
    
    % Initial set
    lb = [1.0; 1.0; 1.0; 1.0];
    % ub = [1.3; 1.3;1.3;1.3];
    ub = [1.001; 1.001; 1.001; 1.001];
    init_set = Star(lb,ub);
    % Input set
    lb = [0;0];
    ub = [0;0];
    input_set = Star(lb,ub);
    % Store all reachable sets
    reachAll = init_set;
    % Execute reachabilty analysis
    num_steps = 6;
    reachOptions.reachMethod = 'approx-star';
    t = tic;
    for i=1:num_steps
        % Compute controller output set
        input_set = net.reach(init_set,reachOptions);
        % Compute plant reachable set
        init_set = plantReach(plant,init_set, input_set,'lin');
        reachAll = [reachAll init_set];
    end
    t = toc(t);
    
    % Save results
    if is_codeocean
        save('/results/logs/double_pendulum_more.mat', 'reachAll','t','-v7.3');
    else
        save('double_pendulum_more.mat','reachAll','t','-v7.3');
    end
    
    %% Visualize results
    plant.get_interval_sets;
    
    f = figure;
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
    hold on;
    rectangle('Position',[-0.5,-0.5,2,2],'FaceColor',[0 0.2 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    grid;
    xlabel('x1');
    ylabel('x2');
    
    f1 = figure;
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
    hold on;
    rectangle('Position',[-0.5,-0.5,2,2],'FaceColor',[0 0.2 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    grid;
    xlabel('x3');
    ylabel('x4');

    % Save figures
    if is_codeocean
        exportgraphics(f,'/results/logs/double_pendulum_more_1v2.pdf', 'ContentType', 'vector');
        exportgraphics(f1,'/results/logs/double_pendulum_more_3v4.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'double_pendulum_more_1v2.pdf','ContentType', 'vector');
        exportgraphics(f1,'double_pendulum_more_3v4.pdf','ContentType', 'vector');
    end

end

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
