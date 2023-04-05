function reach_neuralode_spiral2D_linear()

    % Reachability analysis of a spiral2D defned as a Neural ODE
    
    %% Define layers and neural ODE
    controlPeriod = 10; % total seconds
    reachStep = 0.01; % 1 second
    % Load parameters
    net_info = load('models/odeffnn_spiral.mat');
    Wb = net_info.Wb;
    % Contruct NeuralODE
    % ODEBlock only linear layers
    % Convert in form of a linear ODE model
    states = 2;
    outputs = 2;
    w1 = Wb{3};
    b1 = Wb{4}';
    w2 = Wb{5};
    b2 = Wb{6}';
    Aout = w2*w1;
    Bout = b2' + b1'*w2';
    Cout = eye(states);
    D = zeros(outputs,1);
    numSteps = controlPeriod/reachStep;
    odeblock = LinearODE(Aout,Bout',Cout,D,controlPeriod,numSteps);
    odelayer = ODEblockLayer(odeblock,controlPeriod,reachStep,true);
    neuralode = NN({odelayer});
    
    x0 = [2.0;0.0]; % This is like the initial input to the ODEblock (initial state)
    
    
    %% Reachability 

    % Setup
    Initial_radius = 0.01; % Uncertainty in dynamics (initial states)
    lb = x0 - Initial_radius;
    ub = x0 + Initial_radius;
    R0 = Star(lb,ub);
    % Reach
    t = tic;
    Rall = neuralode.reach(R0); % Reachability
    toc(t);
    
    % Visualize reach sets
    figure;
    hold on;
    Star.plotBoxes_2D_noFill(Rall,1,2,'b'); % plot over-approx of reach sets
    grid;
    % Add legend and labels
    xlabel('x_1');
    ylabel('x_2');
    

end