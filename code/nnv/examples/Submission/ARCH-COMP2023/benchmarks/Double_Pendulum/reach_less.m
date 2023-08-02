function t = reach_less()

    %% Reachability analysis of Double Pendulum Benchmark
    

    %% Load components

    % load the controller
    net = load_NN_from_mat('controller_double_pendulum_less_robust.mat');
    % Specify the reach step, has to be smaller than the control period
    reachStep = 0.005;
    % specify the control period as specified by the benchmark description
    controlPeriod = 0.05;
    % define the plant as specified by nnv
    plant = NonLinearODE(4,2,@dynamics_dp, reachStep, controlPeriod, eye(4));
    
    %% Reachability analysis

    % Initial set
    % lb = [1.0; 1.0;1.0;1.0];
    lb = [1.29; 1.29;1.29;1.29];
    ub = [1.3; 1.3;1.3;1.3];
    % ub = [1.01; 1.01;1.01;1.01];
    init_set = Star(lb,ub);
    % Input set
    lb = [0;0];
    ub = [0;0];
    input_set = Star(lb,ub);
    % Store all reachable sets
    reachAll = init_set;
    % Execute reachabilty analysis
    % for i =1:steps
    num_steps = 4;
    reachOptions.reachMethod = 'approx-star';
    t = tic;
    for i=1:num_steps
        % Compute controller output set
        input_set = net.reach(init_set, reachOptions);
        % Compute plant reachable set
        init_set = plantReach(plant,init_set, input_set,'lin');
        reachAll = [reachAll init_set];
    end
    t = toc(t);
    
    % Save results
    if is_codeocean
        save('/results/logs/double_pendulum_less.mat', 'reachAll','t','-v7.3');
    else
        save('double_pendulum_less.mat','reachAll','t','-v7.3');
    end
    
    %% Visualize results
    plant.get_interval_sets;
    
    f = figure;
    hold on;
    rectangle('Position',[-1,-1,2.7,2.7],'FaceColor',[0 0.2 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
    grid;
    xlabel('x1');
    ylabel('x2');
    
    f1 = figure;
    hold on;
    rectangle('Position',[-1,-1,2.7,2.7],'FaceColor',[0 0.2 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
    grid;
    xlabel('x3');
    ylabel('x4');

    % Save figures
    if is_codeocean
        exportgraphics(f,'/results/logs/double_pendulum_less_1v2.pdf', 'ContentType', 'vector');
        exportgraphics(f1,'/results/logs/double_pendulum_less_3v4.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'double_pendulum_less_1v2.pdf','ContentType', 'vector');
        exportgraphics(f1,'double_pendulum_less_3v4.pdf','ContentType', 'vector');
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