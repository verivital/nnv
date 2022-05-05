function reach_ilnode(nnpath,dynamics, dim, unc, tfl)
    % Reachability analysis of a Neural ODE

    %% Define layers and neural ODE
    % Function defined in a different file for CORA
    if tfl
        controlPeriod = 1.0;
        reachStep = 0.01;
    else
        controlPeriod = 30.0; % total seconds
%         controlPeriod = 0.05; % total seconds
        reachStep = 0.05; % smaller reachStep, more accurate, longer computation
    end
    C = eye(2+dim); % Want to get both of the outputs from NeuralODE
    % Load neural network parameters
    load(nnpath);
    if dim == 3
        layer1 = LayerS(double(w1),double(b1)','purelin');
        layerOut = LayerS(double(w5),double(b5)','purelin');
    else
        layer1 = LayerS(double(Wb{1}),double(Wb{2})','purelin');
        layerOut = LayerS(double(Wb{9}),double(Wb{10})','purelin');
    end
    odeblock = NonLinearODE(2+dim,1,dynamics,reachStep,controlPeriod,C); % Nonlinear ODE plant 
    odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
    neuralode = NeuralODE({layer1, odelayer, layerOut});
%     neuralode = NeuralODE({odelayer});
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
    if dim ~= 2
        try
            t = tic;
            Rall = neuralode.reach(R0); % Reachability
            time = toc(t);
        catch ME
            warning('Reachability of ZonoF failed');
            Rall = ME;
            time = 'NA';
        end
%     pointSets = neuralode.Layers{1,2}.odemodel.intermediate_pointSet;
        save("../../results/"+name+"_zonoF.mat",'time','Rall');
        disp('ZonoF has finished');
        disp("Time elapsed is "+string(time) + " seconds");
    end
    % Plot results
%     f = figure;
%     hold on;
%     Star.plotBoxes_2D_noFill(Rall,1,2,'b');
%     xlabel('x_1');
%     ylabel('x_2');
%     saveas(f,name+"_traj1.png");

    %% Reachability run #2 (Star + zonoA)
    odeblock = NonLinearODE(2+dim,1,dynamics,reachStep,controlPeriod,C); % Nonlinear ODE plant
    % Change default options
    % odeblock.options.timeStep = 0.05;
    % odeblock.options.taylorTerms = 4;
    % odeblock.options.zonotopeOrder = 50;
    odeblock.options.alg = 'lin-adaptive';
    % odeblock.options.tensorOrder = 3;
    odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
    neuralode = NeuralODE({layer1, odelayer, layerOut});
    if dim ~= 2
        try
            t = tic;
            Rall = neuralode.reach(R0); % Reachability
            time = toc(t);
        catch ME
            warning('Reachability of ZonoA failed');
            Rall = ME;
            time = 'NA';
        end
    
        disp('ZonoA has finished');
        disp("Time elapsed is "+string(time) + " seconds");
    %     pointSets = neuralode.Layers{1,2}.odemodel.intermediate_pointSet;
        save("../../results/"+name+"_zonoA.mat",'time','Rall');
    end
    %% Reachability run #3 (Star + polyA)
    odeblock = NonLinearODE(2+dim,1,dynamics,reachStep,controlPeriod,C); % Nonlinear ODE plant
    % Change default options
    % odeblock.options.timeStep = 0.05;
    % odeblock.options.taylorTerms = 4;
    % odeblock.options.zonotopeOrder = 50;
    odeblock.options.alg = 'poly-adaptive';
    % odeblock.options.tensorOrder = 3;
    odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
    neuralode = NeuralODE({layer1, odelayer, layerOut});

    if dim ~= 2
        try
            t = tic;
            Rall = neuralode.reach(R0); % Reachability
            time = toc(t);
        catch ME
            warning('Reachability of PolyzonoA failed');
            Rall = ME;
            time = 'NA';
        end
        disp('PolyA has finished');
        disp("Time elapsed is "+string(time) + " seconds");
    %     pointSets = neuralode.Layers{1,2}.odemodel.intermediate_pointSet;
        save("../../results/"+name+"_polyA.mat",'time','Rall');
    end

    %% Reachability run #3 (Star + polyA)
    odeblock = NonLinearODE(2+dim,1,dynamics,reachStep,controlPeriod,C); % Nonlinear ODE plant
    % Change default options
    % odeblock.options.timeStep = 0.05;
    % odeblock.options.taylorTerms = 4;
    % odeblock.options.zonotopeOrder = 50;
    odeblock.options.alg = 'poly';
    odeblock.options.tensorOrder = 3;
    odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
    neuralode = NeuralODE({layer1, odelayer, layerOut});

    try
        t = tic;
        Rall = neuralode.reach(R0); % Reachability
        time = toc(t);
    catch ME
        warning('Reachability of PolyzonoF failed');
        Rall = ME;
        time = 'NA';
    end
%     pointSets = neuralode.Layers{1,2}.odemodel.intermediate_pointSet;
    save("../../results/"+name+"_polyF.mat",'time','Rall');
    disp('PolyF has finished');
    disp("Time elapsed is "+string(time) + " seconds");

    %% Save results
    % May want to set equal axes so that the set representations are equally
    % visualized
%     save("reach"+string(dim)+".mat",'ta','tb','tc');
end

