function t = reach()

%% Reachability analysis of the aircarft quad benchmark


%% Load components
    net = load_NN_from_mat('model.mat');
    controlPeriod = 0.1;
    reachstep = 0.01;
    plant = NonLinearODE(12,3,@dynamics, reachstep, controlPeriod, eye(12));
    plant.set_taylorTerms(4);
    plant.set_zonotopeOrder(20);
    plant.set_tensorOrder(2);
%     plant.set_intermediateOrder(200);
    % plant.set_polytopeOrder(50);% error = 0.001;
    % error = 0.0005;
    % plant.options.maxError = [error; error; error; error];
    
    %% Reachability analysis

    % Initial set
    % lb = [-0.4; -0.4 ; -0.4; -0.4; -0.4; -0.4; 0; 0; 0; 0; 0; 0];
    % ub = [0.4; 0.4 ; 0.4; 0.4; 0.4; 0.4; 0; 0; 0; 0; 0; 0];
    lb = [-0.01; -0.01; -0.01; -0.01; -0.01; -0.01; 0; 0; 0; 0; 0; 0];
    ub = [0.01; 0.01 ; 0.01; 0.01; 0.01; 0.01; 0; 0; 0; 0; 0; 0];
    init_set = Star(lb,ub);
    % Store all reachable sets
    reachAll = init_set;
    % Execute reachabilty analysis
%     steps = 50;
    steps = 7;
    reachOptions.reachMethod = 'approx-star';
    t = tic;
    for i = 1:steps
        % Compute controller output set
        input_set = net.reach(init_set,reachOptions);
        % Compute plant reachable set
        init_set = plantReach(plant,init_set,input_set,'lin');
        reachAll = [reachAll init_set];
    end
    t = toc(t);
    
    % Save results
    if is_codeocean
        save('/results/logs/quad.mat', 'reachAll','t','-v7.3');
    else
        save('quad.mat', 'reachAll','t','-v7.3');
    end

    %% Visualize results

    plant.get_interval_sets;

    f = figure;
    hold on;
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
    grid;
    xlabel('x1');
    ylabel('x2');
%     saveas(f,[path_out, 'reach1v2.pdf']);
    
    f1 = figure;
    hold on;
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,4,'b');
    grid;
    xlabel('x1');
    ylabel('x2');
%     saveas(f1,[path_out, 'reach3v4.pdf']);

    % Save figures
    if is_codeocean
        exportgraphics(f,'/results/logs/quad_1v2.pdf', 'ContentType', 'vector');
        exportgraphics(f1,'/results/logs/quad_3v4.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'quad_1v2.pdf','ContentType', 'vector');
        exportgraphics(f,'quad_3v4.pdf','ContentType', 'vector');
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