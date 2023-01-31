function reach_ilnode(nnpath,dynamics, dim, unc, tfl)
    % Reachability analysis of a Neural ODE

    %% Define layers and neural ODE
    % Function defined in a different file for CORA
    if tfl
        controlPeriod = 1.0;
        reachStep = 0.01;
    else
        controlPeriod = 30.0; % total seconds
        reachStep = 0.05; % reach time step
    end
    C = eye(2+dim); % Want to get both of the outputs from NeuralODE

    % Load neural network parameters
    load(nnpath);
    if dim == 3
        layer1 = FullyConnectedLayer(double(w1),double(b1)');
        layerOut = FullyConnectedLayer(double(w5),double(b5)');
    else
        layer1 = FullyConnectedLayer(double(Wb{1}),double(Wb{2})');
        layerOut = FullyConnectedLayer(double(Wb{9}),double(Wb{10}'));
    end
    odeblock = NonLinearODE(2+dim,1,dynamics,reachStep,controlPeriod,C); % Nonlinear ODE plant 
    odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
    odeblock.options.tensorOrder = 2;
    neuralode = NN({layer1, odelayer, layerOut});
    if tfl
        name = "DampedOsc_ilnodenl_true_"+string(dim);
    else
        name = "DampedOsc_ilnodenl_"+string(dim);
    end

    %% Reachability run #1 (Star + zonoF)
    % Setup
    x0 = [-1.4996;-0.4609]; % Initial state first trajectory
    lb = x0-unc;
    ub = x0+unc;
    R0 = Star(lb,ub);
    try
        t = tic;
        Rall = neuralode.reach(R0); % Reachability
        time = toc(t);
    catch ME
        warning('Reachability of ZonoF failed');
        Rall = ME;
        time = 'NA';
    end
    save("../../nnvresults/"+name+"_zonoF.mat",'time','Rall');


    %% Reachability run #2 (Star + zonoA)

    odeblock = NonLinearODE(2+dim,1,dynamics,reachStep,controlPeriod,C); % Nonlinear ODE plant
    % Change default options
    % odeblock.options.timeStep = 0.05;
    % odeblock.options.taylorTerms = 4;
    % odeblock.options.zonotopeOrder = 50;
    odeblock.options.alg = 'lin-adaptive';
    odeblock.options.tensorOrder = 2;
    odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
    neuralode = NN({layer1, odelayer, layerOut});

    try
        t = tic;
        Rall = neuralode.reach(R0); % Reachability
        time = toc(t);
    catch ME
        warning('Reachability of ZonoA failed');
        Rall = ME;
        time = 'NA';
    end
    save("../../nnvresults/"+name+"_zonoA.mat",'time','Rall');
    
    
    %% Reachability run #3 (Star + polyA)

    odeblock = NonLinearODE(2+dim,1,dynamics,reachStep,controlPeriod,C); % Nonlinear ODE plant
    % Change default options
    % odeblock.options.timeStep = 0.05;
    % odeblock.options.taylorTerms = 4;
    % odeblock.options.zonotopeOrder = 50;
    odeblock.options.alg = 'poly-adaptive';
    % odeblock.options.tensorOrder = 3;
    odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
    neuralode = NN({layer1, odelayer, layerOut});

    try
        t = tic;
        Rall = neuralode.reach(R0); % Reachability
        time = toc(t);
    catch ME
        warning('Reachability of PolyzonoA failed');
        Rall = ME;
        time = 'NA';
    end
    save("../../nnvresults/"+name+"_polyA.mat",'time','Rall');

    %% Reachability run #3 (Star + polyA)
    odeblock = NonLinearODE(2+dim,1,dynamics,reachStep,controlPeriod,C); % Nonlinear ODE plant
    % Change default options
    % odeblock.options.timeStep = 0.05;
    % odeblock.options.taylorTerms = 4;
    % odeblock.options.zonotopeOrder = 50;
    odeblock.options.alg = 'poly';
    odeblock.options.tensorOrder = 3;
    odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
    neuralode = NN({layer1, odelayer, layerOut});
    try
        t = tic;
        Rall = neuralode.reach(R0); % Reachability
        time = toc(t);
    catch ME
        warning('Reachability of PolyzonoF failed');
        Rall = ME;
        time = 'NA';
    end
    
    save("../../nnvresults/"+name+"_polyF.mat",'time','Rall');

end

