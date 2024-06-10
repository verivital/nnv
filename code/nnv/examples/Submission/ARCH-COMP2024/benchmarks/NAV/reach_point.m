% function rT = reach_point()

    %% Reachability analysis of NAV Benchmark
    
    %% Load Components 

    % Load the controller
    net = importNetworkFromONNX('networks/nn-nav-point.onnx', "InputDataFormats", "BC");
    net = matlab2nnv(net);
    % net.InputSize = 4;
    
    % Load plant
    reachStep = 0.02;
    controlPeriod = 0.2;
    plant = NonLinearODE(4, 2, @dynamics, reachStep, controlPeriod, eye(4));
    plant.set_tensorOrder(2);

    % Create NNCS
    % nncs = NonlinearNNCS(net,plant);

    %% Reachability analysis

    % Initial set
    lb = [2.999; 2.999; 0; 0];
    ub = [3.0; 3.0; 0; 0];
    init_set = Star(lb,ub);

    % Store all reachable sets
    reachAll = init_set;

    % Reachability options
    num_steps = 25;
    reachOptions.reachMethod = 'approx-star';

    % Reachabilty analysis
    % reachPRM.ref_input = input_ref;
    % reachPRM.numSteps = 30;
    % reachPRM.init_set = init_set;
    % reachPRM.numCores = 1;
    % reachPRM.reachMethod = 'approx-star';
    % [R,rT] = nncs.reach(reachPRM);
    % Got an error from mpt when using Gurobi... 
    % Checking GUROBI license ... Did not find an environment variable GRB_LICENSE_FILE with GUROBI license.
    %
    % Error using optimset (line 177)
    % Invalid default value for property 'modules' in class 'mptopt':
    % No default options available for the function 'linprog'.
    % 
    % Error in mpt_solvers_options (line 601)
    %     options.linprog = optimset('linprog');
    % 
    % Error in mpt_subModules (line 28)
    %         options.(modules_list(i).name) = feval(f);
    % 
    % Error in Polyhedron (line 264)
	% 			    MPTOPTIONS = mptopt;
    % 
    % Error in NonlinearNNCS/nextInputSetStar (line 287)
    %                 I1 = Polyhedron('lb', lb, 'ub', ub);
    % 
    % Error in NonlinearNNCS/reach (line 183)
    %                 input_set = obj.nextInputSetStar(fb_I{1});
    % 
    % Error in reach_point (line 41)
    %     [R,rT] = nncs.reach(reachPRM);

    % Begin computation
    t = tic;
    for i=1:num_steps
        % disp(i);

        % Compute controller output set
        input_set = net.reach(init_set,reachOptions);

        % Compute plant reachable set
        init_set = plant.stepReachStar(init_set, input_set, 'poly');
        reachAll = [reachAll init_set];
        % toc(t);

    end

    rT = toc(t);

    R = reachAll;

    % Save results
    if is_codeocean
        save('/results/logs/nav_point.mat', 'R','rT','-v7.3');
    else
        save('nav_point.mat', 'R','rT','-v7.3');
    end


    %% Visualize results
    plant.get_interval_sets;

    f = figure;
    rectangle('Position',[-0.5,-0.5,1,1],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1); % goal region
    hold on;
    rectangle('Position',[1,1,1,1],'FaceColor',[0.7 0 0 0.8], 'EdgeColor','r', 'LineWidth', 0.1); % obstacle
    Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,1,2,'b');
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
        exportgraphics(f,'/results/logs/nav-point.pdf', 'ContentType', 'vector');
    else
        exportgraphics(f,'nav-point.pdf','ContentType', 'vector');
    end

    % end