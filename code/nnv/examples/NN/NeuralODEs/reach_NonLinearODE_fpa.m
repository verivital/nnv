function reach_NonLinearODE_fpa() 

    %% Reachability analysis of the CTRNN Fixed Point Attractor (FPA)
    
    reachstep = 0.01; % step size to compute reach sets
    final_time = 10; % Time horizon
    Initial_radius = 0.01; % Uncertainty in dynamics.
    model = NonLinearODE(5,1,@CTRNN_FPA, reachstep, final_time,eye(5));
    
    % Change default options
    model.options.timeStep = 0.05;
    model.options.taylorTerms = 4;
    model.options.zonotopeOrder = 20;
    model.options.alg = 'lin';
    model.options.tensorOrder = 2;
    
    % Initial states
    x0 = [0.21535, -0.58587, 0.8, 0.52323, 0.5]';
    lb = x0 - Initial_radius;
    ub = x0 + Initial_radius;
    init_set = Star(lb,ub);
    input_set = Star(0,0); % No inputs, but need to define it
    
    % Compute reachability analysis
    t = tic;
    R = model.stepReachStar(init_set,input_set);
    toc(t);

    % Get all sets
    model.get_interval_sets;
    Rall = model.intermediate_reachSet;

    % Visualize results
    figure;
    hold on;
    Star.plotBoxes_2D_noFill(Rall,1,2,'b'); % plot over-approx of reach sets
    grid;
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');

end