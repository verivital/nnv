function reach_neuralode_spiral2D_nonlinear()

    %% Reachability analysis of a spiral2D Neural ODE
    % Function defined in a different file for CORA
    C = eye(2); 
    reachstep = 0.01; % step size to compute reach sets
    final_time = 2.0; % Time horizon
    Initial_radius = 0.01; % Uncertainty in dynamics (initial states)
    % Contruct NeuralODE (odeblock)
    model = NonLinearODE(2,1,@spiral_non,reachstep,final_time,C); % Nonlinear ODE plant
    odelayer = ODEblockLayer(model, final_time, reachstep, true);
    neuralode = NN({odelayer});
    
    % Change default options
    model.options.timeStep = 0.01;
    model.options.taylorTerms = 2;
    model.options.zonotopeOrder = 20;
    model.options.alg = 'lin';
    model.options.tensorOrder = 2;
    
    % Initial states
    x0 = [2.0;0.0]; % This is like the initial input to the ODEblock (initial state)
    lb = x0 - Initial_radius;
    ub = x0 + Initial_radius;
    init_set = Star(lb,ub);
    
    % Compute reachability analysis
    t = tic;
    Rall = neuralode.reach(init_set);
    toc(t);
    
    % Visualize results
    figure;
    hold on;
    Star.plotBoxes_2D_noFill(Rall,1,2,'b'); % plot over-approx of reach sets
    grid;
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');


end