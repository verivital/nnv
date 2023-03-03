function spiral2D_nonlinear_reach()

    %% Reachability analysis of a spiral2D Neural ODE
    % Function defined in a different file for CORA
    C = eye(2); 
    reachstep = 0.01; % step size to compute reach sets
    final_time = 10.0; % Time horizon
    Initial_radius = 0.01; % Uncertainty in dynamics (initial states)
    % Contruct NeuralODE (odeblock)
    model = NonLinearODE(2,1,@spiral_non,reachstep,final_time,C); % Nonlinear ODE plant 
    
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
    input_set = Star(0,0); % No inputs, but need to define it
    
    % Compute reachability analysis
    t = tic;
    R = model.stepReachStar(init_set,input_set);
    time = toc(t);
    Rall = model.intermediate_reachSet;
    save('spiral_nl_0.01.mat','Rall','time')
    
    
    %% Reachability scenario 2

    Initial_radius = 0.05; % Uncertainty in dynamics (initial states)
    % Contruct NeuralODE (odeblock)
    model = NonLinearODE(2,1,@spiral_non,reachstep,final_time,C); % Nonlinear ODE plant 
    
    % Change default options
    model.options.timeStep = 0.01;
    model.options.taylorTerms = 2;
    model.options.zonotopeOrder = 20;
    model.options.alg = 'lin';
    model.options.tensorOrder = 2;
    
    % Initial states
    lb = x0 - Initial_radius;
    ub = x0 + Initial_radius;
    init_set = Star(lb,ub);
    input_set = Star(0,0); % No inputs, but need to define it
    
    % Compute reachability analysis
    t = tic;
    R = model.stepReachStar(init_set,input_set);
    time = toc(t);
    Rall = model.intermediate_reachSet;
    save('spiral_nl_0.05.mat','Rall','time')

    
    %% Reachability scenario #3

    Initial_radius = 0.1; % Uncertainty in dynamics (initial states)
    % Contruct NeuralODE (odeblock)
    model = NonLinearODE(2,1,@spiral_non,reachstep,final_time,C); % Nonlinear ODE plant 
    
    % Change default options
    model.options.timeStep = 0.01;
    model.options.taylorTerms = 2;
    model.options.zonotopeOrder = 20;
    model.options.alg = 'lin';
    model.options.tensorOrder = 2;
    
    % Initial states
    lb = x0 - Initial_radius;
    ub = x0 + Initial_radius;
    init_set = Star(lb,ub);
    input_set = Star(0,0); % No inputs, but need to define it
    
    % Compute reachability analysis
    t = tic;
    R = model.stepReachStar(init_set,input_set);
    time = toc(t);
    Rall = model.intermediate_reachSet;
    save('spiral_nl_0.1.mat','Rall','time')

end