function reach_acc_orig()

    % Reachability analysis of the "original" Adaptive Cruise Control presented in FM2019 

    % safety specification: relative distance > safe distance
    % dis = x_lead - x_ego  
    % safe distance between two cars, see here 
    % https://www.mathworks.com/help/mpc/examples/design-an-adaptive-cruise-control-system-using-model-predictive-control.html
    % dis_safe = D_default + t_gap * v_ego;
    
    %%  Load objects
    
    % Load NN 
    net = load_NN_from_mat('controller_5_20.mat'); % Load controller

    % Load plant
    reachStep = 0.01;
    controlPeriod = 0.1;
    output_mat = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % feedback: relative distance, relative velocity and ego-car velocity
    plant = NonLinearODE(6,1,@dynamicsACC, reachStep, controlPeriod, output_mat);
    plant.options.tensorOrder = 2;

    % Create NNCS object
    nncs = NonlinearNNCS(net,plant);


    %% Reachability analysis

    % Set reachability options
    tF = 5; % seconds
    input_ref = [30;1.4];
    % Initial set
    lb = [90; 32; 0; 10; 30; 0];
    ub = [110; 32.2; 0; 11; 30.2; 0];
    init_set = Star(lb,ub);
    reachPRM.ref_input = input_ref;
    reachPRM.numSteps = 50;
    reachPRM.init_set = init_set;
    reachPRM.numCores = 1;
    reachPRM.reachMethod = 'approx-star';

    % Execute reachabilty analysis
    [R,rT] = nncs.reach(reachPRM);

    % Save results
    save("results_original.mat","R","rT");
    

    %% Visualize results

    % Transform reach sets for visualization
    t_gap = 1.4;
    D_default = 10;
    outAll = [];
    safe_dis = [];
    for i=1:length(R)
        outAll = [outAll R(i).affineMap(output_mat,[])]; % distance between cars
        safe_dis = [safe_dis R(i).affineMap([0 0 0 0 t_gap 0], D_default)]; % safe distance
    end
    times = 0:controlPeriod:tF; % to plot in x-axis (time)

    % Create figure 
    f = figure;
    hold on;
    % Plot results
    pb = plot(0,85,'m');
    pr = plot(0,85,'r');
    Star.plotRanges_2D(outAll,2,times,'r'); % plot overapproximation of reach sets
    hold on;
    Star.plotRanges_2D(safe_dis,1,times,'m'); % plot overapproximation of reach sets

    % Enhance figure for paper
    ax = gca; % Get current axis
    ax.XAxis.FontSize = 15; % Set font size of axis
    ax.YAxis.FontSize = 15;
    xlabel('Time (s)');
    ylabel('Distance (m)')
    legend([pr,pb],{'rel dist','safe dist'},"Location","best",'FontSize',14);
    
    % Save figure
    if is_codeocean
        exportgraphics(f,'/results/logs/acc_orig.pdf','ContentType','vector');
    else
        exportgraphics(f,'acc_orig.pdf','ContentType','vector');
    end

end