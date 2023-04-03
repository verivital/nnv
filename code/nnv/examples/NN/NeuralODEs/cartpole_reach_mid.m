function cartpole_reach_mid()
    
    %% Reachability analysis of the CTRNN Cartpole (12 dimensions)
    
    reachstep = 0.01; % step size to compute reach sets
    final_time = 1.0; % Time horizon
    Initial_radius = 1e-4; % Uncertainty in dynamics.
    model = NonLinearODE(12,1,@CartpoleCTRNN, reachstep, final_time,eye(12));
    
    % Change default options
    model.options.timeStep = 0.01;
    model.options.taylorTerms = 4;
    model.options.zonotopeOrder = 20;
    model.options.alg = 'lin';
    model.options.tensorOrder = 2;
    
    % Initial states
    x0 = [0, 0, 0.001, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
    lb = x0 - Initial_radius;
    ub = x0 + Initial_radius;
    init_set = Star(lb,ub);
    input_set = Star(0,0); % No inputs, but need to define it
    
    % Compute reachability analysis
    t = tic;
    R = model.stepReachStar(init_set,input_set);
    time = toc(t);
    Rall = model.intermediate_reachSet;
    % save results
    save('cartpole_reach_mid.mat','Rall','time')
end